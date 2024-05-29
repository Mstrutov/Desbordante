#include "algorithms/md/hymd/record_pair_inferrer.h"

#include <cstddef>

#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "model/index.h"

namespace algos::hymd {

struct RecordPairInferrer::PairStatistics {
    std::size_t rhss_removed = 0;
    std::size_t all_rhss_removed = 0;
    std::size_t invalidated_number = 0;
};

struct RecordPairInferrer::Statistics {
    std::size_t samplings_started = 0;

    std::size_t pairs_processed = 0;
    std::size_t pairs_inspected = 0;

    std::size_t mds_removed = 0;
    std::size_t all_rhss_removed = 0;

    std::size_t invalidated_number = 0;

    void AddPairStatistics(PairStatistics const& pair_statistics) noexcept {
        mds_removed += pair_statistics.rhss_removed;
        all_rhss_removed += pair_statistics.all_rhss_removed;
        invalidated_number += pair_statistics.invalidated_number;
    }
};

bool RecordPairInferrer::ShouldStopInferring(Statistics const& statistics) const noexcept {
    // NOTE: this is the condition from the original implementation. I believe it severely reduces
    // the benefits of the hybrid approach on datasets where inference from record pairs is
    // efficient.
    // return statistics.samplings_started >= 2 || statistics.pairs_processed > 100;

    // Without this, the algorithm will switch phase after processing the first record pair in the
    // case when only one table is inspected.
    bool const not_enough_data = statistics.pairs_inspected < kPairsRequiredForPhaseSwitch;
    if (not_enough_data) return false;
    // Modified phase switch heuristic described in "Efficient Discovery of Matching Dependencies".
    bool const lattice_is_almost_final = statistics.pairs_inspected * final_lattice_numerator_ >=
                                         final_lattice_denominator_ * statistics.mds_removed;
    // New phase switch heuristic: if there are too many pairs that share similarity classifier
    // boundaries with already processed pairs, switch phase.
    bool const pairs_are_stale = statistics.pairs_inspected * kStaleHeuristicNumerator >=
                                 kStaleHeuristicDenominator * statistics.pairs_processed;
    return lattice_is_almost_final || pairs_are_stale;
}

auto RecordPairInferrer::ProcessPairComparison(PairComparisonResult const& pair_comparison_result)
        -> PairStatistics {
    using MdRefiner = lattice::MdLattice::MdRefiner;
    std::size_t mds_removed = 0;
    std::size_t all_rhss_removed = 0;
    std::size_t total_invalidated = 0;
    std::vector<MdRefiner> refiners = lattice_->CollectRefinersForViolated(pair_comparison_result);
    for (MdRefiner& refiner : refiners) {
        std::size_t const rhss_removed = refiner.Refine();
        mds_removed += rhss_removed;
        std::size_t const invalidated = refiner.InvalidatedNumber();
        total_invalidated += invalidated;
        all_rhss_removed += (rhss_removed == invalidated);
    }
    return {mds_removed, all_rhss_removed, total_invalidated};
}

bool RecordPairInferrer::InferFromRecordPairs(Recommendations recommendations) {
    Statistics statistics;

    auto process_collection = [&](auto& collection, auto get_sim_vec) {
        while (!collection.empty()) {
            if (ShouldStopInferring(statistics)) {
                final_lattice_denominator_ *= kFinalLatticeHeuristicGrowthDenominator;
                final_lattice_numerator_ *= kFinalLatticeHeuristicGrowthNumerator;
                return true;
            }
            PairComparisonResult const& pair_comparison_result =
                    get_sim_vec(collection.extract(collection.begin()).value());
            ++statistics.pairs_inspected;
            if (avoid_same_comparison_processing_) {
                bool const not_seen_before =
                        processed_comparisons_.insert(pair_comparison_result).second;
                if (!not_seen_before) continue;
            }
            statistics.AddPairStatistics(ProcessPairComparison(pair_comparison_result));
            ++statistics.pairs_processed;
        }
        return false;
    };
    if (process_collection(recommendations, [&](Recommendation& rec) {
            // TODO: parallelize similarity vector calculation
            return similarity_data_->CompareRecords(*rec.left_record, *rec.right_record);
        })) {
        return false;
    }
    auto move_out = [&](PairComparisonResult& pair_comp_res) { return std::move(pair_comp_res); };
    if (process_collection(comparisons_to_process_, move_out)) {
        return false;
    }
    std::size_t const left_size = similarity_data_->GetLeftSize();
    while (next_left_record_ < left_size) {
        ++statistics.samplings_started;
        comparisons_to_process_ = similarity_data_->CompareAllWith(next_left_record_);
        ++next_left_record_;
        if (process_collection(comparisons_to_process_, move_out)) {
            return false;
        }
    }
    return true;
}

}  // namespace algos::hymd
