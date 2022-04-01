#pragma once

#include "Primitive.h"
#include "FDAlgorithm.h"
#include "ColumnLayoutTypedRelationData.h"
#include "ColumnLayoutRelationData.h"

namespace algos {

BETTER_ENUM(Metric, char,
            euclidian = 0,
            levenshtein
)

class MetricVerifier : public algos::Primitive {
private:
    Metric metric_;
    std::vector<unsigned int> lhs_indices_;
    unsigned int rhs_index_;
    double parameter_;
    bool dist_to_null_infinity_;

    bool metric_fd_holds_ = false;

    std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;
    std::shared_ptr<ColumnLayoutRelationData> relation_; // temporarily parsing twice

    bool CompareNumericValues(
        util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const;
    bool CompareStringValues(
        util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const;
    bool VerifyMetricFD(model::TypedColumnData const& col) const;
public:
    struct Config {
        std::filesystem::path data{}; /* Path to input file */
        char separator = ',';         /* Separator for csv */
        bool has_header = true;       /* Indicates if input file has header */
        bool is_null_equal_null = true; /* Is NULL value equals another NULL value */
        std::string metric;
        double parameter;
        std::vector<unsigned int> lhs_indices;
        unsigned int rhs_index;
        bool dist_to_null_infinity;
    };

    unsigned long long Execute() override;

    bool GetResult() const {
        return metric_fd_holds_;
    }

    void SetParameter(double parameter) {
        parameter_ = parameter;
    }

    explicit MetricVerifier(Config const& config);

    explicit MetricVerifier(Config const& config,
                            std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation,
                            std::shared_ptr<ColumnLayoutRelationData> relation);
};

}  // namespace algos
