#pragma once

#include "Primitive.h"
#include "FDAlgorithm.h"
#include "ColumnLayoutTypedRelationData.h"
#include "ColumnLayoutRelationData.h"
#include "QGramVector.h"

namespace algos {

BETTER_ENUM(Metric, char,
            euclidian = 0,  /* Standard metric for calculating the distance between numeric
                             * values */
            levenshtein,    /* Levenshtein distance between strings */
            cosine          /* Represent strings as q-gram vectors and calculate cosine distance
                             * between these vectors */
)

class MetricVerifier : public algos::Primitive {
private:
    Metric metric_;
    std::vector<unsigned int> lhs_indices_;
    unsigned int rhs_index_;
    double parameter_;
    unsigned int q_;
    bool dist_to_null_infinity_;

    bool metric_fd_holds_ = false;

    std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;
    std::shared_ptr<ColumnLayoutRelationData> relation_; // temporarily parsing twice

    bool CompareNumericValues(
        util::PLI::Cluster const& cluster, model::TypedColumnData const& col) const;
    bool CompareStringValues(
        util::PLI::Cluster const& cluster,
        model::TypedColumnData const& col,
        std::function<double(std::byte const*, std::byte const*)> const& distance_function) const;
    std::function<double(std::byte const*, std::byte const*)> GetCosineDistFunction(
        model::StringType const& type,
        std::unordered_map<std::string, util::QGramVector>& q_gram_map) const;
    bool VerifyMetricFD(model::TypedColumnData const& col) const;

public:
    struct Config {
        std::filesystem::path data{};   /* Path to input file */
        char separator = ',';           /* Separator for csv */
        bool has_header = true;         /* Indicates if input file has header */
        bool is_null_equal_null = true; /* Is NULL value equals another NULL value */
        std::string metric;             /* Metric to verify metric FD */
        double parameter;               /* The max possible distance between 2 values in cluster
                                         * at which metric FD holds */
        std::vector<unsigned int> lhs_indices; /* Indices of LHS columns */
        unsigned int rhs_index;         /* Index of RHS column */
        bool dist_to_null_infinity;     /* Determines whether distance to NULL value is infinity
                                         * or zero */
        unsigned int q;                 /* Q-gram length for cosine metric */
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
