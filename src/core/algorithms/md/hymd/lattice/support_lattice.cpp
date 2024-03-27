#include "algorithms/md/hymd/lattice/support_lattice.h"

#include <cassert>

#include "algorithms/md/hymd/utility/get_first_non_zero_index.h"

namespace {
using algos::hymd::utility::GetFirstNonZeroIndex;
using model::Index, model::md::DecisionBoundary;
}  // namespace

namespace algos::hymd::lattice {

bool SupportLattice::IsUnsupported(Node const& cur_node, DecisionBoundaryVector const& lhs_bounds,
                                   Index this_node_index) const {
    if (cur_node.is_unsupported) return true;
    NodeChildren const& children = cur_node.children;
    std::size_t const child_array_size = children.size();
    for (Index child_array_index = FindFirstNonEmptyIndex(children, 0);
         child_array_index != child_array_size;
         child_array_index = FindFirstNonEmptyIndex(children, child_array_index + 1)) {
        Index const next_node_index = this_node_index + child_array_index;
        DecisionBoundary const generalization_boundary_limit = lhs_bounds[next_node_index];
        for (auto const& [generalization_boundary, node] : *children[child_array_index]) {
            if (generalization_boundary > generalization_boundary_limit) break;
            if (IsUnsupported(node, lhs_bounds, next_node_index + 1)) return true;
        }
    }
    return false;
}

void SupportLattice::MarkNewLhs(Node& cur_node, DecisionBoundaryVector const& lhs_bounds,
                                Index cur_node_index) {
    assert(IsEmpty(cur_node.children));
    std::size_t const col_match_number = lhs_bounds.size();
    Node* cur_node_ptr = &cur_node;
    for (Index next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != col_match_number; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        std::size_t const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = col_match_number - next_node_index;
        cur_node_ptr = &cur_node_ptr->children[child_array_index]
                                .emplace()
                                .try_emplace(lhs_bounds[next_node_index], next_child_array_size)
                                .first->second;
    }
    cur_node_ptr->is_unsupported = true;
}

bool SupportLattice::IsUnsupported(DecisionBoundaryVector const& lhs_bounds) const {
    return IsUnsupported(root_, lhs_bounds, 0);
}

void SupportLattice::MarkUnsupported(DecisionBoundaryVector const& lhs_bounds) {
    std::size_t const col_match_number = lhs_bounds.size();
    Node* cur_node_ptr = &root_;
    for (Index cur_node_index = 0,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != col_match_number; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        DecisionBoundary const next_bound = lhs_bounds[next_node_index];
        Index const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = col_match_number - next_node_index;
        auto [boundary_map, is_first_arr] =
                TryEmplaceChild(cur_node_ptr->children, child_array_index);
        if (is_first_arr) {
            Node& new_node =
                    boundary_map.try_emplace(next_bound, next_child_array_size).first->second;
            MarkNewLhs(new_node, lhs_bounds, next_node_index + 1);
            return;
        }
        auto [it_map, is_first_map] = boundary_map.try_emplace(next_bound, next_child_array_size);
        Node& next_node = it_map->second;
        if (is_first_map) {
            MarkNewLhs(next_node, lhs_bounds, next_node_index + 1);
            return;
        }
        cur_node_ptr = &next_node;
    }
    // Can only happen if the root is unsupported.
    cur_node_ptr->is_unsupported = true;
}

}  // namespace algos::hymd::lattice
