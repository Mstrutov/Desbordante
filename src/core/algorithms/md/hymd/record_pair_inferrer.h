#pragma once

#include <list>
#include <unordered_set>

#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/pair_comparison_result.h"
#include "algorithms/md/hymd/recommendation.h"
#include "algorithms/md/hymd/similarity_data.h"

namespace algos::hymd {

class RecordPairInferrer {
private:
    struct Statistics;
    struct PairStatistics;

    SimilarityData* const similarity_data_;
    lattice::MdLattice* const lattice_;

    // Metanome uses a linked list for some reason.
    std::unordered_set<PairComparisonResult> comparisons_to_process_;
    std::unordered_set<PairComparisonResult> processed_comparisons_;

    util::WorkerThreadPool* pool_;

    RecordIdentifier next_left_record_ = 0;

    static constexpr std::size_t kPairsRequiredForPhaseSwitch = 5;

    std::size_t final_lattice_numerator_ = 1;
    // Represents the reciprocal of the ratio in 5.3.1 of the "Efficient Discovery of Matching
    // Dependencies" article, except instead of refined MDs, removed MDs are used in the numerator
    // (denominator in the reciprocal).
    std::size_t final_lattice_denominator_ = 1;
    static constexpr std::size_t kFinalLatticeHeuristicGrowthNumerator = 1;
    static constexpr std::size_t kFinalLatticeHeuristicGrowthDenominator = 2;

    static constexpr std::size_t kStaleHeuristicNumerator = 1;
    static constexpr std::size_t kStaleHeuristicDenominator = 21;

    bool const avoid_same_comparison_processing_ = true;

    PairStatistics ProcessPairComparison(PairComparisonResult const& pair_comparison_result);
    bool ShouldStopInferring(Statistics const& statistics) const noexcept;

public:
    RecordPairInferrer(SimilarityData* similarity_data, lattice::MdLattice* lattice,
                       util::WorkerThreadPool* pool) noexcept
        : similarity_data_(similarity_data), lattice_(lattice), pool_(pool) {}

    bool InferFromRecordPairs(Recommendations recommendations);
};

}  // namespace algos::hymd
