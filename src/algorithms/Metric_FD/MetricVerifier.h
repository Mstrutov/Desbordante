#pragma once

#include "Primitive.h"
#include "FDAlgorithm.h"
#include "ColumnLayoutTypedRelationData.h"
#include "ColumnLayoutRelationData.h"

class MetricVerifier : public algos::Primitive {
private:
    using Config = FDAlgorithm::Config;
    Config config_;

    std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;
    std::shared_ptr<ColumnLayoutRelationData> relation_; //temporarily parsing 2 times

    bool CompareNumericValues
        (util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const;
    static bool PositiveDoubleLessOrEq(long double l, long double r);
public:
    std::vector<unsigned int> lhs_indices_;
    unsigned int rhs_index_;
    long double parameter_; //long double temporarily

    bool metric_fd_holds = false;

    unsigned long long Execute() override;

    explicit MetricVerifier(Config const& config)
        : Primitive(config.data, config.separator, config.has_header, {}),
          lhs_indices_(config.GetSpecialParam<std::vector<unsigned int>>("lhs_indices")),
          rhs_index_(config.GetSpecialParam<unsigned int>("rhs_index")),
          parameter_(config.GetSpecialParam<long double>("parameter")),
          config_(config) {
        relation_ =
            ColumnLayoutRelationData::CreateFrom(input_generator_, config_.is_null_equal_null);
        CSVParser generator2(config_.data, config_.separator, config_.has_header);
        typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(generator2,
                                                                           config_.is_null_equal_null);
    }
};

