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

void Specializer::SpecializeFor(SimilarityVector const& sim, DecisionBoundaryVector& lhs_bounds,
                                model::Index rhs_index, model::md::DecisionBoundary old_rhs_bound) {
    using model::md::DecisionBoundary, model::Index;
    auto const add_lhs_specialized_md = [this, &lhs_bounds, old_rhs_bound, rhs_index, &sim](
                                                Index i, auto should_add) {
        DecisionBoundary& lhs_sim = lhs_bounds[i];
        std::optional<DecisionBoundary> new_lhs_value = SpecializeOneLhs(i, sim[i]);
        if (!new_lhs_value.has_value()) return;
        auto context = utility::SetForScope(lhs_sim, *new_lhs_value);
        if (should_add(old_rhs_bound, lhs_sim)) {
            lattice_->AddIfMinimalAndNotUnsupported(lhs_bounds, old_rhs_bound, rhs_index);
        }
    };
    auto true_func = [](...) { return true; };
    for (Index i = 0; i < rhs_index; ++i) {
        add_lhs_specialized_md(i, true_func);
    }
    if (!prune_nondisjoint_) {
        add_lhs_specialized_md(rhs_index, std::greater<>{});
    }
    std::size_t const col_match_number = lhs_bounds.size();
    for (Index i = rhs_index + 1; i < col_match_number; ++i) {
        add_lhs_specialized_md(i, true_func);
    }
}

void Specializer::SpecializeInvalidated(DecisionBoundaryVector& lhs_bounds,
                                        InvalidatedRhss const& invalidated) {
    using model::md::DecisionBoundary, model::Index;
    auto specialize_all_lhs = [this, col_matches_num = lhs_bounds.size(), &lhs_bounds,
                               it_begin = invalidated.begin(),
                               it_end = invalidated.end()](auto handle_same_lhs_as_rhs) {
        for (Index lhs_spec_index = 0; lhs_spec_index < col_matches_num; ++lhs_spec_index) {
            DecisionBoundary& lhs_bound = lhs_bounds[lhs_spec_index];
            std::optional<DecisionBoundary> const specialized_lhs_bound =
                    SpecializeOneLhs(lhs_spec_index, lhs_bound);
            if (!specialized_lhs_bound.has_value()) continue;
            auto context = utility::SetForScope(lhs_bound, *specialized_lhs_bound);
            if (lattice_->IsUnsupported(lhs_bounds)) continue;

            for (auto it = it_begin; it != it_end; ++it) {
                auto const& [rhs_index, old_rhs_bound, _] = *it;
                if (rhs_index == lhs_spec_index) {
                    handle_same_lhs_as_rhs(old_rhs_bound, *specialized_lhs_bound, rhs_index);
                    for (++it; it != it_end; ++it) {
                        auto const& [rhs_index, old_rhs_bound, _] = *it;
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
