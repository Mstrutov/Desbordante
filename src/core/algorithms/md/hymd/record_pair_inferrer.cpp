#include "algorithms/md/hymd/record_pair_inferrer.h"

#include <cstddef>

#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "model/index.h"

namespace algos::hymd {

struct RecordPairInferrer::Statistics {
    std::size_t samplings_started = 0;
    std::size_t sim_vecs_processed = 0;
};

bool RecordPairInferrer::ShouldStopInferring(Statistics const& statistics) const noexcept {
    return statistics.samplings_started >= 2 || statistics.sim_vecs_processed > 100;
    /*
    return records_checked >= 5 &&
           (mds_refined == 0 || records_checked / mds_refined >= efficiency_reciprocal_);
    */
}

void RecordPairInferrer::ProcessSimVec(SimilarityVector const& sim) {
    using MdRefiner = lattice::MdLattice::MdRefiner;
    std::vector<MdRefiner> refiners = lattice_->CollectRefinersForViolated(sim);
    for (MdRefiner& refiner : refiners) {
        refiner.Refine();
    }
}

bool RecordPairInferrer::InferFromRecordPairs(Recommendations recommendations) {
    Statistics statistics;

    auto process_collection = [&](auto& collection, auto get_sim_vec) {
        while (!collection.empty()) {
            if (ShouldStopInferring(statistics)) {
                // efficiency_reciprocal_ *= 2;
                return true;
            }
            SimilarityVector const sim =
                    get_sim_vec(collection.extract(collection.begin()).value());
            if (avoid_same_sim_vec_processing_) {
                bool const not_seen_before = checked_sim_vecs_.insert(sim).second;
                if (!not_seen_before) continue;
            }
            ProcessSimVec(sim);
            ++statistics.sim_vecs_processed;
        }
        return false;
    };
    if (process_collection(recommendations, [&](Recommendation& rec) {
            // TODO: parallelize similarity vector calculation
            return similarity_data_->GetSimilarityVector(*rec.left_record, *rec.right_record);
        })) {
        return false;
    }
    auto move_out = [&](SimilarityVector& v) { return std::move(v); };
    if (process_collection(sim_vecs_to_check_, move_out)) {
        return false;
    }
    std::size_t const left_size = similarity_data_->GetLeftSize();
    while (next_left_record_ < left_size) {
        ++statistics.samplings_started;
        sim_vecs_to_check_ = similarity_data_->GetSimVecs(next_left_record_);
        ++next_left_record_;
        if (process_collection(sim_vecs_to_check_, move_out)) {
            return false;
        }
    }
    return true;
}

}  // namespace algos::hymd
