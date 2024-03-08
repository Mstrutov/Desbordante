#pragma once

#include <memory>

#include "algorithms/md/hymd/preprocessing/similarity_measure/similarity_measure.h"
#include "algorithms/md/hymd/similarity_measure_creator.h"

namespace algos::hymd::preprocessing::similarity_measure {

class LevenshteinSimilarityMeasure final : public SimilarityMeasure {
    bool const is_null_equal_null_;
    model::md::DecisionBoundary const min_sim_;
    std::size_t const size_limit_;
    static constexpr auto kName = "levenshtein_similarity";

public:
    class Creator final : public SimilarityMeasureCreator {
        model::md::DecisionBoundary const min_sim_;
        bool const is_null_equal_null_;
        std::size_t const size_limit_;

    public:
        Creator(model::md::DecisionBoundary min_sim = 0.7, bool is_null_equal_null = true,
                std::size_t size_limit = 0);

        std::unique_ptr<SimilarityMeasure> MakeMeasure() const final {
            return std::make_unique<LevenshteinSimilarityMeasure>(min_sim_, is_null_equal_null_,
                                                                  size_limit_);
        }
    };

    [[nodiscard]] indexes::ColumnMatchSimilarityInfo MakeIndexes(
            std::shared_ptr<DataInfo const> data_info_left,
            std::shared_ptr<DataInfo const> data_info_right,
            std::vector<indexes::PliCluster> const& clusters_right,
            util::WorkerThreadPool& thread_pool) const final;

    LevenshteinSimilarityMeasure(model::md::DecisionBoundary min_sim, bool is_null_equal_null,
                                 std::size_t size_limit);
};

}  // namespace algos::hymd::preprocessing::similarity_measure
