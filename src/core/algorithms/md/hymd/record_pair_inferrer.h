#pragma once

#include <list>
#include <unordered_set>

#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/pair_comparison_result.h"
#include "algorithms/md/hymd/recommendation.h"
#include "algorithms/md/hymd/similarity_data.h"
#include "util/ratio.h"

namespace algos::hymd {

class RecordPairInferrer {
public:
    struct PairStatistics;

    struct PhaseSwitchHeuristicParameters {
        static constexpr std::size_t kPairsRequiredForPhaseSwitch = 5;

        // Represents the ratio in 5.3.1 of the "Efficient Discovery of Matching Dependencies"
        // article, except instead of refined MDs, removed MDs are used in the numerator.
        util::Ratio<std::size_t> final_lattice_ratio;
        static constexpr util::Ratio<std::size_t> kStaleRatio{1, 21};
        static constexpr util::Ratio<std::size_t> kFinalLatticeMult{1, 2};

        // v1 * ratio >= v2 ?
        static bool IsGe(std::size_t v1, util::Ratio<std::size_t> ratio, std::size_t v2) noexcept {
            auto const& [numer, denom] = ratio;
            return v1 * numer >= v2 * denom;
        }
    };

private:
    SimilarityData* const similarity_data_;
    lattice::MdLattice* const lattice_;

    // Metanome uses a linked list for some reason.
    std::unordered_set<PairComparisonResult> comparisons_to_process_;
    std::unordered_set<PairComparisonResult> processed_comparisons_;

    util::WorkerThreadPool* pool_;

    RecordIdentifier next_left_record_ = 0;

    PhaseSwitchHeuristicParameters heuristic_parameters{{1, 1}};

    bool const avoid_same_comparison_processing_ = true;

    PairStatistics ProcessPairComparison(PairComparisonResult const& pair_comparison_result);
    template <typename StatisticsType>
    bool ShouldStopInferring(StatisticsType const& statistics) const noexcept;

public:
    RecordPairInferrer(SimilarityData* similarity_data, lattice::MdLattice* lattice,
                       util::WorkerThreadPool* pool) noexcept
        : similarity_data_(similarity_data), lattice_(lattice), pool_(pool) {}

    bool InferFromRecordPairs(Recommendations recommendations);
};

}  // namespace algos::hymd
