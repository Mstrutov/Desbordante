#include "MetricVerifier.h"

#include <chrono>
#include <functional>
#include <easylogging++.h>

unsigned long long MetricVerifier::Execute() {
    LOG(INFO) << "parameter: " << parameter_ << ", RHS index: " << rhs_index_ << ", LHS indices:";
    for (unsigned int& index : lhs_indices_) {
        LOG(INFO) << index;
    }
    auto start_time = std::chrono::system_clock::now();

    std::shared_ptr<util::PLI const>
        pli = relation_->GetColumnData(lhs_indices_[0]).GetPliOwnership();

    for (size_t i = 1; i < lhs_indices_.size(); ++i) {
        pli = pli->Intersect(relation_->GetColumnData(lhs_indices_[i]).GetPositionListIndex());
    }

    model::TypedColumnData const& col = typed_relation_->GetColumnData(rhs_index_);
    LOG(INFO) << "Column Type: " << col.GetTypeId();

    metric_fd_holds = true;
    std::function<bool(util::PLI::Cluster const& cluster)> compareFunction =
        [](util::PLI::Cluster const& cluster) { return false; };
    if (col.IsNumeric()) {
        compareFunction = [this, &col](util::PLI::Cluster const& cluster) {
            return CompareNumericValues(cluster, col);
        };
    } else {
        metric_fd_holds = false;
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
(util::PLI::Cluster const & cluster, model::TypedColumnData const & col) const {
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
    }
    std::unique_ptr<std::byte[]> res(type.Allocate());
    type.Sub(max_value, min_value, res.get());
    switch (type.GetTypeId()) {
        case model::TypeId::kInt:
            return std::abs(model::INumericType::GetValue<model::Int>(res.get())) <= parameter_;
        case model::TypeId::kDouble:
            return PositiveDoubleLessOrEq(std::fabs(model::INumericType::GetValue<model::Double>(res.get())),
                                          parameter_);
        default:return false;
    }
}

bool MetricVerifier::PositiveDoubleLessOrEq(long double l, long double r) {
    return l - r <= std::numeric_limits<double>::epsilon();
}
