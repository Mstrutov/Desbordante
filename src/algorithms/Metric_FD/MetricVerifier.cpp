#include "MetricVerifier.h"

#include <chrono>
#include <easylogging++.h>

namespace algos {

unsigned long long MetricVerifier::Execute() {
    auto start_time = std::chrono::system_clock::now();

    std::shared_ptr<util::PLI const>
        pli = relation_->GetColumnData(lhs_indices_[0]).GetPliOwnership();

    for (size_t i = 1; i < lhs_indices_.size(); ++i) {
        pli = pli->Intersect(relation_->GetColumnData(lhs_indices_[i]).GetPositionListIndex());
    }

    model::TypedColumnData const& col = typed_relation_->GetColumnData(rhs_index_);

    metric_fd_holds = true;
    std::function<bool(util::PLI::Cluster const& cluster)> compareFunction;
    switch (metric_) {
        case Metric::integer:
        case Metric::floating_point:
            compareFunction = [this, &col](util::PLI::Cluster const& cluster) {
                return CompareNumericValues(cluster, col);
            };
            break;
        case Metric::levenshtein:
            compareFunction = [this, &col](util::PLI::Cluster const& cluster) {
                return CompareStringValues(cluster, col);
            };
            break;
    }

    for (util::PLI::Cluster const& cluster : pli->GetIndex()) {
        if (!compareFunction(cluster)) {
            metric_fd_holds = false;
            break;
        }
    }

    if (metric_fd_holds) {
        LOG(INFO) << "Metric fd holds.";
    } else {
        LOG(INFO) << "Metric fd does not hold.";
    }
    auto elapsed_milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );
    return elapsed_milliseconds.count();
}

bool MetricVerifier::CompareNumericValues
    (util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const {
    auto const& type = dynamic_cast<model::INumericType const&>(col.GetType());
    std::vector<std::byte const*> const& data = col.GetData();
    std::byte const* max_value = nullptr;
    std::byte const* min_value = nullptr;
    for (int row_index : cluster) {
        if (col.IsNull(row_index) || col.IsEmpty(row_index)) {
            continue;
        }
        if (max_value == nullptr) {
            max_value = data[row_index];
            min_value = data[row_index];
            continue;
        }
        if (type.Compare(data[row_index], max_value) == model::CompareResult::kGreater) {
            max_value = data[row_index];
        } else if (type.Compare(data[row_index], min_value) == model::CompareResult::kLess) {
            min_value = data[row_index];
        }
        if (type.Dist(max_value, min_value) > parameter_)
            return false;
    }
    return true;
}

bool MetricVerifier::CompareStringValues
    (const util::PLI::Cluster& cluster, const model::TypedColumnData& col) const {
    auto const& type = dynamic_cast<model::StringType const&>(col.GetType());
    std::vector<std::byte const*> const& data = col.GetData();
    for (size_t i = 0; i < cluster.size() - 1; ++i) {
        if (col.IsNull(cluster[i]) || col.IsEmpty(cluster[i]))
            continue;
        double max_dist = 0;
        for (size_t j = i + 1; j < cluster.size(); ++j) {
            if (col.IsNull(cluster[j]) || col.IsEmpty(cluster[j]))
                continue;
            max_dist = fmax(max_dist, type.Dist(data[cluster[i]], data[cluster[j]]));
            if (max_dist > parameter_) return false;
        }
        if (max_dist * 2 <= parameter_) return true;
    }
    return true;
}

} // namespace algos