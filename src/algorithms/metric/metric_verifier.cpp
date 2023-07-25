#include "algorithms/metric/metric_verifier.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <easylogging++.h>

#include "util/config/equal_nulls/option.h"
#include "util/config/indices/option.h"
#include "util/config/names_and_descriptions.h"
#include "util/config/option_using.h"
#include "util/config/tabular_data/input_table/option.h"

namespace algos::metric {

MetricVerifier::MetricVerifier() : Algorithm({}) {
    RegisterOptions();
}

void MetricVerifier::ValidateRhs(util::config::IndicesType const& rhs_indices) {
    assert(!rhs_indices.empty());
    if (rhs_indices.size() == 1) {
        util::config::IndexType column_index = rhs_indices[0];
        model::TypedColumnData const& column = typed_relation_->GetColumnData(column_index);
        model::TypeId type_id = column.GetTypeId();
        if (type_id == +model::TypeId::kUndefined) {
            throw std::invalid_argument("Column with index \"" + std::to_string(column_index)
                                        + "\" type undefined.");
        }
        if (type_id == +model::TypeId::kMixed) {
            throw std::invalid_argument("Column with index \"" + std::to_string(column_index)
                                        + "\" contains values of different types.");
        }

        if (metric_ == +Metric::euclidean) {
            if (!column.IsNumeric()) {
                throw std::invalid_argument("\"Euclidean\" metric is only available for numeric "
                                            "columns.");
            }
            return;
        }
        if (type_id == +model::TypeId::kString) return;
        throw std::invalid_argument("The chosen metric is available only for string columns.");
    }
    if (metric_ == +Metric::euclidean) {
        for (util::config::IndexType column_index : rhs_indices) {
            model::TypedColumnData const& column = typed_relation_->GetColumnData(column_index);
            model::TypeId type_id = column.GetTypeId();
            if (type_id == +model::TypeId::kUndefined) {
                throw std::invalid_argument("Column with index \"" + std::to_string(column_index)
                                            + "\" type undefined.");
            }
            if (type_id == +model::TypeId::kMixed) {
                throw std::invalid_argument("Column with index \"" + std::to_string(column_index)
                                            + "\" contains values of different types.");
            }

            if (!column.IsNumeric()) {
                throw std::invalid_argument("\"Euclidean\" metric is only available for numeric "
                                            "columns, column with index "
                                            + std::to_string(column_index) + " is not numeric");
            }
        }
        return;
    }
    throw std::invalid_argument("Multidimensional RHS is not available for the chosen metric");
}

void MetricVerifier::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    auto check_parameter = [](long double parameter) {
        if (parameter < 0) throw std::invalid_argument("Parameter out of range");
    };
    auto get_schema_columns = [this]() { return relation_->GetSchema()->GetNumColumns(); };
    auto check_rhs = [this](util::config::IndicesType const& rhs_indices) {
        ValidateRhs(rhs_indices);
    };
    auto need_algo_and_q = [this]([[maybe_unused]] util::config::IndicesType const& _) {
        return metric_ == +Metric::cosine;
    };
    auto need_algo_only = [this](util::config::IndicesType const& rhs_indices) {
        assert(metric_ == +Metric::levenshtein || metric_ == +Metric::euclidean);
        return metric_ == +Metric::levenshtein || rhs_indices.size() != 1;
    };
    auto algo_check = [this](MetricAlgo metric_algo) {
        assert(!(metric_ == +Metric::euclidean && rhs_indices_.size() == 1));
        if (metric_algo == +MetricAlgo::calipers) {
            if (!(metric_ == +Metric::euclidean && rhs_indices_.size() == 2))
                throw std::invalid_argument("\"calipers\" algorithm is only available for "
                                            "2-dimensional RHS and \"euclidean\" metric.");
        }
    };
    auto q_check = [](unsigned int q) {
        if (q <= 0) throw std::invalid_argument("Q-gram length should be greater than zero.");
    };

    RegisterInitialLoadOption(util::config::TableOpt(&input_table_));
    RegisterInitialLoadOption(util::config::EqualNullsOpt(&is_null_equal_null_));
    RegisterInitialExecOption(util::config::LhsIndicesOpt(&lhs_indices_, get_schema_columns));
    RegisterInitialExecOption(Option{&dist_from_null_is_infinity_, kDistFromNullIsInfinity,
                                     kDDistFromNullIsInfinity, false});
    RegisterInitialExecOption(
            Option{&parameter_, kParameter, kDParameter}.SetValueCheck(check_parameter));
    RegisterInitialExecOption(Option{&metric_, kMetric, kDMetric}.SetConditionalOpts(
            {{{}, {util::config::RhsIndicesOpt.GetName()}}}));
    RegisterOption(Option{&algo_, kMetricAlgorithm, kDMetricAlgorithm}.SetValueCheck(algo_check));
    RegisterOption(Option{&q_, kQGramLength, kDQGramLength, 2u}.SetValueCheck(q_check));
    RegisterOption(util::config::RhsIndicesOpt(&rhs_indices_, get_schema_columns, check_rhs)
                           .SetConditionalOpts({{need_algo_and_q, {kMetricAlgorithm, kQGramLength}},
                                                {need_algo_only, {kMetricAlgorithm}}}));
}

void MetricVerifier::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);
    input_table_->Reset();
    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: metric FD verifying is meaningless.");
    }
    typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, is_null_equal_null_);
}

void MetricVerifier::ResetState() {
    metric_fd_holds_ = false;
}

unsigned long long MetricVerifier::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();

    points_calculator_ = std::make_unique<PointsCalculator>(dist_from_null_is_infinity_,
                                                            typed_relation_, rhs_indices_);
    highlight_calculator_ = std::make_unique<HighlightCalculator>(typed_relation_, rhs_indices_);
    assert(points_calculator_.get() != nullptr || highlight_calculator_.get() != nullptr);

    VerifyMetricFD();

    if (metric_fd_holds_) {
        LOG(DEBUG) << "Metric fd holds.";
    } else {
        LOG(DEBUG) << "Metric fd does not hold.";
    }

    SortHighlightsByDistanceDescending();

    VisualizeHighlights();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

std::string MetricVerifier::GetStringValue(util::config::IndicesType const& index_vec,
                                           ClusterIndex row_index) const {
    model::TypedColumnData const& col = typed_relation_->GetColumnData(index_vec[0]);
    if (col.IsNull(row_index)) {
        return "NULL";
    }
    if (col.IsEmpty(row_index)) {
        return "EMPTY";
    }
    if (index_vec.size() == 1) {
        return col.GetType().ValueToString(col.GetData()[row_index]);
    }
    std::string value("(");
    for (size_t j = 0; j < index_vec.size(); ++j) {
        model::TypedColumnData const& coord_col = typed_relation_->GetColumnData(index_vec[j]);
        value += coord_col.GetType().ValueToString(coord_col.GetData()[row_index]);
        if (j == index_vec.size() - 1) {
            break;
        }
        value += ", ";
    }
    value.push_back(')');
    return value;
}

void MetricVerifier::VisualizeHighlights() const {
    if (highlight_calculator_->GetHighlights().empty()) {
        return;
    }
    for (auto const& cluster_highlight : highlight_calculator_->GetHighlights()) {
        LOG(DEBUG) << "----------------------------------------- LHS value: "
                   << GetStringValue(lhs_indices_, cluster_highlight[0].data_index);
        for (auto const& highlight : cluster_highlight) {
            bool is_null =
                    typed_relation_->GetColumnData(rhs_indices_[0]).IsNull(highlight.data_index);
            bool is_empty =
                    typed_relation_->GetColumnData(rhs_indices_[0]).IsEmpty(highlight.data_index);
            std::string value = GetStringValue(rhs_indices_, highlight.data_index);
            std::string begin_desc;
            if (!is_empty) {
                begin_desc = std::string("[") + (highlight.max_distance <= parameter_ ? "✓" : "X") +
                             "] | max dist: " + std::to_string(highlight.max_distance) + "\t| ";
            }
            std::string end_desc;
            if (!is_null && !is_empty) {
                end_desc = "\t| furthest point index: " +
                           std::to_string(highlight.furthest_data_index) +
                           "\t| furthest point value: " +
                           GetStringValue(rhs_indices_, highlight.furthest_data_index);
            }
            LOG(DEBUG) << begin_desc << "index: " << highlight.data_index << "\t| value: " << value
                       << end_desc;
        }
    }
}

void MetricVerifier::VerifyMetricFD() {
    std::shared_ptr<util::PLI const> pli =
            relation_->GetColumnData(lhs_indices_[0]).GetPliOwnership();

    for (size_t i = 1; i < lhs_indices_.size(); ++i) {
        pli = pli->Intersect(relation_->GetColumnData(lhs_indices_[i]).GetPositionListIndex());
    }

    metric_fd_holds_ = true;
    auto cluster_func = GetClusterFunction();
    for (auto const& cluster : pli->GetIndex()) {
        if (!cluster_func(cluster)) {
            metric_fd_holds_ = false;
            if (algo_ == +MetricAlgo::approx) {
                break;
            }
        }
    }
}

ClusterFunction MetricVerifier::GetClusterFunctionForOneDimension() {
    model::TypedColumnData const& col = typed_relation_->GetColumnData(rhs_indices_[0]);

    if (metric_ == +Metric::euclidean) {
        assert(col.IsNumeric());
        return CalculateClusterFunction<IndexedOneDimensionalPoint>(
                [this](auto const& cluster) {
                    return points_calculator_->CalculateIndexedPoints(cluster);
                },
                [this](auto const& points) { return CompareNumericValues(points); },
                [this](auto const& points, std::vector<Highlight>&& cluster_highlights) {
                    return highlight_calculator_->CalculateOneDimensionalHighlights(
                            points, std::move(cluster_highlights));
                });
    }

    assert(col.GetTypeId() == +model::TypeId::kString);
    auto const& type = static_cast<model::StringType const&>(col.GetType());

    std::function<ClusterFunction(DistanceFunction<std::byte const*>)> verify_func;
    if (algo_ == +MetricAlgo::brute) {
        verify_func = [this](auto dist_func) {
            return CalculateClusterFunction<IndexedOneDimensionalPoint>(
                    [this](auto const& cluster) {
                        return points_calculator_->CalculateIndexedPoints(cluster);
                    },
                    [this, dist_func](auto const& points) {
                        return this->BruteVerifyCluster(points, dist_func);
                    },
                    [this, dist_func](auto const& points,
                                      std::vector<Highlight>&& cluster_highlights) {
                        return highlight_calculator_->CalculateHighlightsForStrings(
                                points, std::move(cluster_highlights), dist_func);
                    });
        };
    } else {
        verify_func = [this](auto const& dist_func) {
            return CalculateApproxClusterFunction<std::byte const*>(
                    [this](auto const& cluster) {
                        return points_calculator_->CalculatePoints(cluster);
                    },
                    dist_func);
        };
    }

    if (metric_ == +Metric::levenshtein) {
        return verify_func(
                [&type](std::byte const* l, std::byte const* r) { return type.Dist(l, r); });
    }

    return [this, &type, verify_func](util::PLI::Cluster const& cluster) {
        std::unordered_map<std::string, util::QGramVector> q_gram_map;
        return verify_func(GetCosineDistFunction(type, q_gram_map))(cluster);
    };
}

ClusterFunction MetricVerifier::GetClusterFunctionForSeveralDimensions() {
    if (algo_ == +MetricAlgo::calipers) {
        return [this](util::PLI::Cluster const& cluster) {
            auto result = points_calculator_->CalculateMultidimensionalPointsForCalipers(cluster);
            if (!CheckMFDFailIfHasNulls(result.has_nulls) &&
                CalipersCompareNumericValues(result.points)) {
                return true;
            }

            auto result_indexed =
                    points_calculator_->CalculateMultidimensionalIndexedPoints(cluster);
            highlight_calculator_->CalculateMultidimensionalHighlights(
                    result_indexed.points, std::move(result_indexed.cluster_highlights));
            return false;
        };
    }
    if (algo_ == +MetricAlgo::brute) {
        return CalculateClusterFunction<IndexedVector>(
                [this](auto const& cluster) {
                    return points_calculator_->CalculateMultidimensionalIndexedPoints(cluster);
                },
                [this](auto const& points) {
                    return BruteVerifyCluster<std::vector<long double>>(points,
                                                                        util::EuclideanDistance);
                },
                [this](auto const& points, std::vector<Highlight>&& cluster_highlights) {
                    return highlight_calculator_->CalculateMultidimensionalHighlights(
                            points, std::move(cluster_highlights));
                });
    }
    auto points_func = [this](auto const& cluster) {
        return points_calculator_->CalculateMultidimensionalPointsForApprox(cluster);
    };
    return CalculateApproxClusterFunction<std::vector<long double>>(points_func,
                                                                    util::EuclideanDistance);
}

ClusterFunction MetricVerifier::GetClusterFunction() {
    if (rhs_indices_.size() == 1) {
        return GetClusterFunctionForOneDimension();
    }
    return GetClusterFunctionForSeveralDimensions();
}

DistanceFunction<std::byte const*> MetricVerifier::GetCosineDistFunction(
        model::StringType const& type,
        std::unordered_map<std::string, util::QGramVector>& q_gram_map) const {
    return [this, &type, &q_gram_map](std::byte const* a, std::byte const* b) -> long double {
        std::string str1 = type.ValueToString(a);
        std::string str2 = type.ValueToString(b);
        if (str1.length() < q_ || str2.length() < q_) {
            throw std::runtime_error(
                    "q-gram length should not exceed the minimum string length "
                    "in the dataset.");
        }
        util::QGramVector const& v1 = q_gram_map.try_emplace(str1, str1, q_).first->second;
        util::QGramVector const& v2 = q_gram_map.try_emplace(str2, str2, q_).first->second;
        return v1.CosineDistance(v2);
    };
}

template <typename T>
ClusterFunction MetricVerifier::CalculateClusterFunction(
        IndexedPointsFunction<T> points_func, CompareFunction<T> compare_func,
        HighlightFunction<T> highlight_func) const {
    return [this, points_func, compare_func, highlight_func](util::PLI::Cluster const& cluster) {
        auto result = points_func(cluster);
        if (!CheckMFDFailIfHasNulls(result.has_nulls) && compare_func(result.points)) {
            return true;
        }
        highlight_func(result.points, std::move(result.cluster_highlights));
        return false;
    };
}

template <typename T>
ClusterFunction MetricVerifier::CalculateApproxClusterFunction(
        PointsFunction<T> points_func, DistanceFunction<T> dist_func) const {
    return [points_func, dist_func, this](util::PLI::Cluster const& cluster) {
        auto result = points_func(cluster);
        return !CheckMFDFailIfHasNulls(result.has_nulls) &&
               ApproxVerifyCluster(result.points, dist_func);
    };
}

bool MetricVerifier::CompareNumericValues(
        std::vector<IndexedOneDimensionalPoint> const& points) const {
    if (points.size() < 2) {
        return true;
    }
    model::TypedColumnData const& col = typed_relation_->GetColumnData(rhs_indices_[0]);
    auto const& type = static_cast<model::INumericType const&>(col.GetType());

    std::byte const* max_value = points[0].point;
    std::byte const* min_value = points[0].point;
    for (size_t i = 1; i < points.size(); ++i) {
        if (type.Compare(points[i].point, max_value) == model::CompareResult::kGreater) {
            max_value = points[i].point;
        } else if (type.Compare(points[i].point, min_value) == model::CompareResult::kLess) {
            min_value = points[i].point;
        }
        if (type.Dist(max_value, min_value) > parameter_) {
            return false;
        }
    }
    return true;
}

template <typename T>
bool MetricVerifier::ApproxVerifyCluster(std::vector<T> const& points,
                                         DistanceFunction<T> const& dist_func) const {
    if (points.size() < 2) {
        return true;
    }
    return std::all_of(std::next(points.cbegin()), points.cend(),
                       [this, &points, &dist_func](auto const& p) {
                           return dist_func(points[0], p) * 2 <= parameter_;
                       });
}

template <typename T>
bool MetricVerifier::BruteVerifyCluster(std::vector<IndexedPoint<T>> const& points,
                                        DistanceFunction<T> const& dist_func) const {
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            if (dist_func(points[i].point, points[j].point) > parameter_) {
                return false;
            }
        }
    }
    return true;
}

bool MetricVerifier::CalipersCompareNumericValues(std::vector<util::Point>& points) const {
    auto pairs = util::GetAntipodalPairs(util::CalculateConvexHull(points));
    return std::all_of(pairs.cbegin(), pairs.cend(), [this](auto const& pair) {
        return util::Point::EuclideanDistance(pair.first, pair.second) <= parameter_;
    });
}

}  // namespace algos::metric
