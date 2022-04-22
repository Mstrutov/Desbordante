#include "MetricVerifier.h"

#include <chrono>

#include "easylogging++.h"

namespace algos {

MetricVerifier::MetricVerifier(Config const& config)
    : Primitive(config.data, config.separator, config.has_header, {}),
      metric_(Metric::_from_string(config.metric.c_str())),
      lhs_indices_(config.lhs_indices),
      rhs_index_(config.rhs_index),
      parameter_(config.parameter),
      q_(config.q),
      dist_to_null_infinity_(config.dist_to_null_infinity) {
    relation_ =
        ColumnLayoutRelationData::CreateFrom(input_generator_, config.is_null_equal_null);
    input_generator_.Reset();
    typed_relation_ =
        model::ColumnLayoutTypedRelationData::CreateFrom(input_generator_,
                                                         config.is_null_equal_null);
}

MetricVerifier::MetricVerifier(Config const& config,
                               std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation,
                               std::shared_ptr<ColumnLayoutRelationData> relation)
    : Primitive(config.data, config.separator, config.has_header, {}),
      metric_(Metric::_from_string(config.metric.c_str())),
      lhs_indices_(config.lhs_indices),
      rhs_index_(config.rhs_index),
      parameter_(config.parameter),
      q_(config.q),
      dist_to_null_infinity_(config.dist_to_null_infinity),
      typed_relation_(std::move(typed_relation)),
      relation_(std::move(relation)) {}

unsigned long long MetricVerifier::Execute() {
    auto start_time = std::chrono::system_clock::now();

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: metric FD verifying is meaningless.");
    }

    size_t cols_count = relation_->GetSchema()->GetNumColumns();
    std::sort(lhs_indices_.begin(), lhs_indices_.end());
    lhs_indices_.erase(std::unique(lhs_indices_.begin(), lhs_indices_.end()), lhs_indices_.end());

    bool indices_out_of_range = rhs_index_ >= cols_count
        || std::any_of(lhs_indices_.begin(), lhs_indices_.end(),
                       [cols_count](unsigned int i) { return i >= cols_count; });
    if (indices_out_of_range) {
        throw std::runtime_error(
            "Column index should be less than the number of columns in the dataset.");
    }

    model::TypedColumnData const& col = typed_relation_->GetColumnData(rhs_index_);
    auto type_id = col.GetTypeId();

    if (type_id == +model::TypeId::kMixed) {
        throw std::runtime_error("Column contains values of different types.");
    }
    metric_fd_holds_ = type_id == +model::TypeId::kUndefined
                       ? !dist_to_null_infinity_
                       : VerifyMetricFD(col);

    if (metric_fd_holds_) {
        LOG(INFO) << "Metric fd holds.";
    } else {
        LOG(INFO) << "Metric fd does not hold.";
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

bool MetricVerifier::VerifyMetricFD(model::TypedColumnData const& col) const {
    std::shared_ptr<util::PLI const>
        pli = relation_->GetColumnData(lhs_indices_[0]).GetPliOwnership();

    for (size_t i = 1; i < lhs_indices_.size(); ++i) {
        pli = pli->Intersect(relation_->GetColumnData(lhs_indices_[i]).GetPositionListIndex());
    }

    std::function<bool(util::PLI::Cluster const& cluster)> compare_function;
    switch (metric_) {
    case Metric::euclidian:
        if (!col.IsNumeric()) {
            throw std::runtime_error("\"euclidian\" metric does not match RHS column type");
        }
        compare_function = [this, &col](util::PLI::Cluster const& cluster) {
            return CompareNumericValues(cluster, col);
        };
        break;
    case Metric::levenshtein:
        if (col.GetTypeId() != +model::TypeId::kString) {
            throw std::runtime_error("\"levenshtein\" metric does not match RHS column type");
        }
        compare_function = [this, &col](util::PLI::Cluster const& cluster) {
            auto const& type = dynamic_cast<model::StringType const&>(col.GetType());
            return CompareStringValues(
                cluster, col, [&type](std::byte const* a, std::byte const* b) {
                    return type.Dist(a, b);
                });
        };
        break;
    case Metric::cosine:
        if (col.GetTypeId() != +model::TypeId::kString) {
            throw std::runtime_error("\"cosine\" metric does not match RHS column type");
        }
        compare_function = [this, &col](util::PLI::Cluster const& cluster) {
            auto const& type = dynamic_cast<model::StringType const&>(col.GetType());
            std::unordered_map<std::string, util::QGramVector> q_gram_map;
            return CompareStringValues(cluster, col, GetCosineDistFunction(type, q_gram_map));
        };
        break;
    }

    return std::all_of(pli->GetIndex().cbegin(), pli->GetIndex().cend(), compare_function);
}

std::function<double(std::byte const*, std::byte const*)> MetricVerifier::GetCosineDistFunction(
    model::StringType const& type,
    std::unordered_map<std::string, util::QGramVector>& q_gram_map) const {
    return [this, &type, &q_gram_map](std::byte const* a, std::byte const* b) -> double {
        std::string str1 = type.ValueToString(a);
        std::string str2 = type.ValueToString(b);
        if (str1.length() < q_ || str2.length() < q_) {
            throw std::runtime_error("q-gram length should not exceed the minimum string length "
                                     "in the dataset");
        }
        util::QGramVector& v1 = q_gram_map.try_emplace(str1, str1, q_).first->second;
        util::QGramVector& v2 = q_gram_map.try_emplace(str2, str2, q_).first->second;
        return v1.CosineDistance(v2);
    };
}

bool MetricVerifier::CompareNumericValues(
    util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const {
    auto const& type = dynamic_cast<model::INumericType const&>(col.GetType());
    std::vector<std::byte const*> const& data = col.GetData();
    std::byte const* max_value = nullptr;
    std::byte const* min_value = nullptr;
    for (int row_index : cluster) {
        if (col.IsNull(row_index) || col.IsEmpty(row_index)) {
            if (dist_to_null_infinity_) return false;
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

bool MetricVerifier::CompareStringValues(
    util::PLI::Cluster const& cluster,
    model::TypedColumnData const& col,
    std::function<double(std::byte const*, std::byte const*)> const& distance_function) const {
    std::vector<std::byte const*> const& data = col.GetData();
    for (size_t i = 0; i < cluster.size() - 1; ++i) {
        if (col.IsNull(cluster[i]) || col.IsEmpty(cluster[i])) {
            if (dist_to_null_infinity_) return false;
            continue;
        }
        for (size_t j = i + 1; j < cluster.size(); ++j) {
            if (col.IsNull(cluster[j]) || col.IsEmpty(cluster[j])) {
                if (dist_to_null_infinity_) return false;
                continue;
            }
            if (distance_function(data[cluster[i]], data[cluster[j]]) > parameter_) return false;
        }
    }
    return true;
}

}  // namespace algos
