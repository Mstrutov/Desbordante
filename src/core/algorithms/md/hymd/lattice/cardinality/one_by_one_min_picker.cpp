#include "algorithms/md/hymd/lattice/cardinality/one_by_one_min_picker.h"

#include <algorithm>
#include <cassert>

namespace algos::hymd::lattice::cardinality {

void OneByOnePicker::NewBatch(std::size_t elements) {
    currently_picked_.clear();
    currently_picked_.reserve(elements);
}

void OneByOnePicker::AddGeneralizations(MdLatticeNodeInfo& node_info,
                                        boost::dynamic_bitset<>& considered_indices) {
    DecisionBoundaryVector const& lhs_bounds_cur = node_info.lhs_bounds;
    auto cur_begin = lhs_bounds_cur.begin();
    auto cur_end = lhs_bounds_cur.end();
    for (ValidationInfo& prev_info : currently_picked_) {
        DecisionBoundaryVector const& lhs_bounds_prev = prev_info.node_info->lhs_bounds;
        boost::dynamic_bitset<>& indices_prev = prev_info.rhs_indices;
        auto cur_it = cur_begin;
        auto prev_it = lhs_bounds_prev.begin();
        while (cur_it != cur_end) {
            model::md::DecisionBoundary const cur_bound = *cur_it;
            model::md::DecisionBoundary const prev_bound = *prev_it;
            if (cur_bound < prev_bound) {
                for (++cur_it, ++prev_it; cur_it != cur_end; ++cur_it, ++prev_it) {
                    model::md::DecisionBoundary const cur_bound = *cur_it;
                    model::md::DecisionBoundary const prev_bound = *prev_it;
                    if (cur_bound > prev_bound) {
                        goto incomparable;
                    }
                }
                indices_prev -= considered_indices;
                break;
            } else if (cur_bound > prev_bound) {
                for (++cur_it, ++prev_it; cur_it != cur_end; ++cur_it, ++prev_it) {
                    model::md::DecisionBoundary const cur_bound = *cur_it;
                    model::md::DecisionBoundary const prev_bound = *prev_it;
                    if (cur_bound < prev_bound) {
                        goto incomparable;
                    }
                }
                considered_indices -= indices_prev;
                if (considered_indices.none()) return;
                break;
            } else {
                ++cur_it;
                ++prev_it;
            }
        }
    incomparable:;
    }
    assert(!considered_indices.none());
    currently_picked_.emplace_back(&node_info, std::move(considered_indices));
}

std::vector<ValidationInfo> OneByOnePicker::GetAll() noexcept {
    static_assert(kNeedsEmptyRemoval,
                  "This picker needs post-processing to remove candidates with empty RHS indices");
    return std::move(currently_picked_);
}

}  // namespace algos::hymd::lattice::cardinality
