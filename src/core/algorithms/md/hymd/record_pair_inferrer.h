#pragma once

#include <list>
#include <unordered_set>

#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/pair_comparison_result.h"
#include "algorithms/md/hymd/recommendation.h"
#include "algorithms/md/hymd/similarity_data.h"

namespace algos::hymd {

class RecordPairInferrer {
public:
    struct PairStatistics;

    struct PhaseSwitchHeuristicParameters {
        static constexpr std::size_t kPairsRequiredForPhaseSwitch = 5;

        // Represents the ratio in 5.3.1 of the "Efficient Discovery of Matching Dependencies"
        // article, except instead of refined MDs, removed MDs are used in the numerator.
        double final_lattice_coefficient;
        static constexpr double kStaleCoefficient = 1 / 21.;
        static constexpr double kFinalLatticeGrowth = 1 / 2.;
    };

private:
    enum class InferenceStatus {
        KeepGoing,
        LatticeIsAlmostFinal,
        PairsAreStale,
    };

    SimilarityData* const similarity_data_;
    lattice::MdLattice* const lattice_;

    // Metanome uses a linked list for some reason.
    std::unordered_set<PairComparisonResult> comparisons_to_process_;
    std::unordered_set<PairComparisonResult> processed_comparisons_;

    util::WorkerThreadPool* pool_;

    RecordIdentifier next_left_record_ = 0;

    PhaseSwitchHeuristicParameters heuristic_parameters{1 / 100.};

    bool const avoid_same_comparison_processing_ = true;

    PairStatistics ProcessPairComparison(PairComparisonResult const& pair_comparison_result);
    template <typename StatisticsType>
    InferenceStatus Evaluate(StatisticsType const& statistics) const noexcept;

public:
    RecordPairInferrer(SimilarityData* similarity_data, lattice::MdLattice* lattice,
                       util::WorkerThreadPool* pool) noexcept
        : similarity_data_(similarity_data), lattice_(lattice), pool_(pool) {}

    bool InferFromRecordPairs(Recommendations recommendations);
};

}  // namespace algos::hymd
