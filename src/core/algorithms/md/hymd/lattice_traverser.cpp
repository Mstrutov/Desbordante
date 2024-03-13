#include "algorithms/md/hymd/lattice_traverser.h"

#include "algorithms/md/hymd/invalidated_rhs.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/utility/set_for_scope.h"
#include "model/index.h"

namespace algos::hymd {

void LatticeTraverser::LowerAndSpecialize(Validator::Result& validation_result,
                                          lattice::ValidationInfo& validation_info) {
    DecisionBoundaryVector& lhs_bounds = validation_info.node_info->lhs_bounds;
    DecisionBoundaryVector& rhs_bounds = *validation_info.node_info->rhs_bounds;

    bool const is_unsupported = validation_result.is_unsupported;
    // TODO: move the below to another class.
    if (is_unsupported) {
        // Does practically nothing, just repeating what is in the article.
        std::fill(rhs_bounds.begin(), rhs_bounds.end(), kLowestBound);
        // TODO: specializations can be removed from the MD lattice. If not worth it, removing just
        // this node and its children should be cheap. Though, destructors also take time.
        lattice_->MarkUnsupported(lhs_bounds);
        return;
    }

    InvalidatedRhss const& invalidated = validation_result.invalidated;
    Rhss rhss;
    rhss.reserve(invalidated.size());
    for (auto const& [index, old_bound, actual_bound] : invalidated) {
        rhs_bounds[index] = actual_bound;
        rhss.emplace_back(index, old_bound);
    }
    lattice_->Specialize(lhs_bounds, lhs_bounds, rhss);
}

bool LatticeTraverser::TraverseLattice(bool const traverse_all) {
    while (level_getter_->AreLevelsLeft()) {
        std::vector<lattice::ValidationInfo> mds = level_getter_->GetCurrentMds();
        if (mds.empty()) {
            continue;
        }

        std::size_t const mds_size = mds.size();
        std::vector<Validator::Result> results(mds_size);
        auto validate_at_index = [&](model::Index i) { results[i] = validator_.Validate(mds[i]); };
        pool_->ExecIndex(validate_at_index, mds_size);
        pool_->WorkUntilComplete();
        auto viol_func = [this, &results]() {
            for (Validator::Result& result : results) {
                for (std::vector<Recommendation> const& rhs_violations : result.recommendations) {
                    recommendations_.insert(rhs_violations.begin(), rhs_violations.end());
                };
            }
        };
        pool_->ExecSingle(viol_func);
        for (model::Index i = 0; i < mds_size; ++i) {
            LowerAndSpecialize(results[i], mds[i]);
        }
        pool_->WorkUntilComplete();
        if (!traverse_all) return false;
    }
    return true;
}

}  // namespace algos::hymd
