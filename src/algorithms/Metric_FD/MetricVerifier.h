#pragma once

#include "Primitive.h"
#include "FDAlgorithm.h"
#include "ColumnLayoutTypedRelationData.h"
#include "ColumnLayoutRelationData.h"

namespace algos {

BETTER_ENUM(Metric, char,
            integer = 0,
            floating_point,
            levenshtein
);

class MetricVerifier : public algos::Primitive {
private:
    std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;
    std::shared_ptr<ColumnLayoutRelationData> relation_; //temporarily parsing 2 times

    bool CompareNumericValues
        (util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const;
    bool CompareStringValues
        (util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const;
public:
    struct Config {
        std::filesystem::path data{}; /* Path to input file */
        char separator = ',';         /* Separator for csv */
        bool has_header = true;       /* Indicates if input file has header */
        bool is_null_equal_null = true; /* Is NULL value equals another NULL value */
        std::string metric;
        long double parameter;
        std::vector<unsigned int> lhs_indices;
        unsigned int rhs_index;
    };

    Metric metric_;
    std::vector<unsigned int> lhs_indices_;
    unsigned int rhs_index_;
    double parameter_;

    bool metric_fd_holds = false;

    unsigned long long Execute() override;

    explicit MetricVerifier(Config const& config)
        : Primitive(config.data, config.separator, config.has_header, {}),
          metric_(Metric::_from_string(config.metric.c_str())),
          lhs_indices_(config.lhs_indices),
          rhs_index_(config.rhs_index),
          parameter_(config.parameter) {
        relation_ =
            ColumnLayoutRelationData::CreateFrom(input_generator_, config.is_null_equal_null);
        CSVParser generator2(config.data, config.separator, config.has_header);
        typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(generator2, config.is_null_equal_null);
    }
};

} // namespace algos