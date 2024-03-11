#include "algorithms/md/hymd/specializer.h"

#include <functional>

#include "algorithms/md/hymd/utility/set_for_scope.h"

namespace algos::hymd {

std::optional<model::md::DecisionBoundary> Specializer::SpecializeOneLhs(
        model::Index col_match_index, model::md::DecisionBoundary lhs_bound) const {
    std::vector<model::md::DecisionBoundary> const& decision_bounds =
            (*column_matches_info_)[col_match_index].similarity_info.lhs_bounds;
    auto end_bounds = decision_bounds.end();
    auto upper = std::upper_bound(decision_bounds.begin(), end_bounds, lhs_bound);
    if (upper == end_bounds) return std::nullopt;
    return *upper;
}

void Specializer::Specialize(DecisionBoundaryVector& lhs_bounds,
                             DecisionBoundaryVector const& specialize_past, Rhss const& rhss) {
    using model::md::DecisionBoundary, model::Index;
    auto specialize_all_lhs = [this, col_matches_num = lhs_bounds.size(), &lhs_bounds,
                               it_begin = rhss.begin(), it_end = rhss.end(),
                               &specialize_past](auto handle_same_lhs_as_rhs) {
        for (Index lhs_spec_index = 0; lhs_spec_index < col_matches_num; ++lhs_spec_index) {
            std::optional<DecisionBoundary> const specialized_lhs_bound =
                    SpecializeOneLhs(lhs_spec_index, specialize_past[lhs_spec_index]);
            if (!specialized_lhs_bound.has_value()) continue;
            auto context = utility::SetForScope(lhs_bounds[lhs_spec_index], *specialized_lhs_bound);
            if (lattice_->IsUnsupported(lhs_bounds)) continue;

            for (auto it = it_begin; it != it_end; ++it) {
                auto const& [rhs_index, old_rhs_bound] = *it;
                if (rhs_index == lhs_spec_index) {
                    handle_same_lhs_as_rhs(old_rhs_bound, *specialized_lhs_bound, rhs_index);
                    for (++it; it != it_end; ++it) {
                        auto const& [rhs_index, old_rhs_bound] = *it;
                        lattice_->AddIfMinimal(lhs_bounds, old_rhs_bound, rhs_index);
                    }
                    break;
                }
                lattice_->AddIfMinimal(lhs_bounds, old_rhs_bound, rhs_index);
            }
        }
    };
    if (prune_nondisjoint_) {
        specialize_all_lhs([](...) {});
    } else {
        specialize_all_lhs([this, &lhs_bounds](DecisionBoundary old_rhs_bound,
                                               DecisionBoundary specialized_lhs_bound,
                                               Index rhs_index) {
            if (old_rhs_bound > specialized_lhs_bound) {
                lattice_->AddIfMinimal(lhs_bounds, old_rhs_bound, rhs_index);
            }
        });
    }
}

}  // namespace algos::hymd
