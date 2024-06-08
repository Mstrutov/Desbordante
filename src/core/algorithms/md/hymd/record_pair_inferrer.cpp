#include "algorithms/md/hymd/record_pair_inferrer.h"

#include <cstddef>
#include <optional>

#include <boost/circular_buffer.hpp>

#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "model/index.h"

namespace algos::hymd {

constexpr std::size_t kSmallestWindowSize = 1000;
static_assert(kSmallestWindowSize > 0);

struct RecordPairInferrer::PairStatistics {
    std::size_t rhss_removed = 0;
    std::size_t all_rhss_removed = 0;
    std::size_t invalidated_number = 0;
};

struct Statistics {
    std::size_t samplings_started = 0;

    std::size_t pairs_processed = 0;
    std::size_t pairs_inspected = 0;

    std::size_t mds_removed = 0;
    std::size_t all_rhss_removed = 0;

    std::size_t invalidated_number = 0;

    void AddPairStats(RecordPairInferrer::PairStatistics pair_statistics) noexcept {
        mds_removed += pair_statistics.rhss_removed;
        all_rhss_removed += pair_statistics.all_rhss_removed;
        invalidated_number += pair_statistics.invalidated_number;
    }
};

struct WindowStatistics : public Statistics {
    using OptionalStats = std::optional<RecordPairInferrer::PairStatistics>;

    boost::circular_buffer<OptionalStats> pairs;

    WindowStatistics(std::size_t window_size) : pairs(window_size) {}

    void AddPairStatistics(RecordPairInferrer::PairStatistics&& pair_statistics) {
        if (pairs.full()) {
            if (OptionalStats const& first_el = *pairs.begin(); first_el.has_value()) {
                mds_removed -= first_el->rhss_removed;
                all_rhss_removed -= first_el->all_rhss_removed;
                invalidated_number -= first_el->invalidated_number;
            } else {
                ++pairs_processed;
            }
        } else {
            ++pairs_processed;
            ++pairs_inspected;
        }
        AddPairStats(pair_statistics);
        pairs.push_back(std::move(pair_statistics));
    }

    void AddSkipped() {
        if (pairs.full()) {
            if (OptionalStats const& first_el = *pairs.begin(); first_el.has_value()) {
                --pairs_processed;
                mds_removed -= first_el->rhss_removed;
                all_rhss_removed -= first_el->all_rhss_removed;
                invalidated_number -= first_el->invalidated_number;
            }
        } else {
            ++pairs_inspected;
        }
        pairs.push_back();
    }
};

struct TotalStatistics : public Statistics {
    void AddPairStatistics(RecordPairInferrer::PairStatistics&& pair_statistics) noexcept {
        ++pairs_inspected;
        ++pairs_processed;
        AddPairStats(pair_statistics);
    }

    void AddSkipped() noexcept {
        ++pairs_inspected;
    }
};

template <typename StatisticsType>
auto RecordPairInferrer::Evaluate(StatisticsType const& statistics) const noexcept
        -> InferenceStatus {
    // NOTE: this is the condition from the original implementation. I believe it severely reduces
    // the benefits of the hybrid approach on datasets where inference from record pairs is
    // efficient.
    // return statistics.samplings_started >= 2 || statistics.pairs_processed > 100;

    // Without this, the algorithm will switch phase after processing the first record pair in the
    // case when only one table is inspected.
    bool const not_enough_data = statistics.pairs_inspected <
                                 PhaseSwitchHeuristicParameters::kPairsRequiredForPhaseSwitch;
    if (not_enough_data) return InferenceStatus::KeepGoing;
    // Modified phase switch heuristic described in "Efficient Discovery of Matching Dependencies".
    bool const lattice_is_almost_final =
            statistics.pairs_inspected * heuristic_parameters.final_lattice_coefficient >=
            statistics.mds_removed;
    if (lattice_is_almost_final) return InferenceStatus::LatticeIsAlmostFinal;

    // New phase switch heuristic: if there are too many pairs that share similarity classifier
    // boundaries with already processed pairs, switch phase.
    bool const pairs_are_stale =
            statistics.pairs_inspected * PhaseSwitchHeuristicParameters::kStaleCoefficient >=
            statistics.pairs_processed;
    if (pairs_are_stale) return InferenceStatus::PairsAreStale;
    return InferenceStatus::KeepGoing;
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
    TotalStatistics statistics{};

    auto process_collection = [&](auto& collection, auto get_sim_vec) {
        while (!collection.empty()) {
            switch (Evaluate(statistics)) {
                case InferenceStatus::KeepGoing:
                    break;
                case InferenceStatus::LatticeIsAlmostFinal:
                    heuristic_parameters.final_lattice_coefficient *=
                            PhaseSwitchHeuristicParameters::kFinalLatticeGrowth;
                    [[fallthrough]];
                case InferenceStatus::PairsAreStale:
                    return true;
                default:
                    assert(false);
                    __builtin_unreachable();
            }
            PairComparisonResult const& pair_comparison_result =
                    get_sim_vec(collection.extract(collection.begin()).value());
            if (avoid_same_comparison_processing_) {
                bool const not_seen_before =
                        processed_comparisons_.insert(pair_comparison_result).second;
                if (!not_seen_before) {
                    statistics.AddSkipped();
                    continue;
                }
            }
            statistics.AddPairStatistics(ProcessPairComparison(pair_comparison_result));
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
