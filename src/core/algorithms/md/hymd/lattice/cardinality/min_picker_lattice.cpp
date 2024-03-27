#include "algorithms/md/hymd/lattice/cardinality/min_picker_lattice.h"

#include <cassert>

#include "algorithms/md/hymd/utility/get_first_non_zero_index.h"

namespace {
using namespace algos::hymd;
using model::Index;
using utility::GetFirstNonZeroIndex;
}  // namespace

namespace algos::hymd::lattice::cardinality {

void MinPickerLattice::AddNewLhs(Node& cur_node, ValidationInfo* validation_info,
                                 Index cur_node_index) {
    assert(IsEmpty(cur_node.children));
    DecisionBoundaryVector const& lhs_bounds = validation_info->node_info->lhs_bounds;
    size_t const col_match_number = lhs_bounds.size();
    Node* cur_node_ptr = &cur_node;
    for (Index next_node_index = utility::GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != col_match_number; cur_node_index = next_node_index + 1,
               next_node_index = utility::GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        std::size_t const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = col_match_number - next_node_index;
        cur_node_ptr = &cur_node_ptr->children[child_array_index]
                                .emplace()
                                .try_emplace(lhs_bounds[next_node_index], next_child_array_size)
                                .first->second;
    }
    cur_node_ptr->task_info = validation_info;
}

void MinPickerLattice::ExcludeGeneralizationRhs(Node const& cur_node,
                                                MdLatticeNodeInfo const& lattice_node_info,
                                                model::Index cur_node_index,
                                                boost::dynamic_bitset<>& considered_indices) {
    if (cur_node.task_info != nullptr) {
        boost::dynamic_bitset<> const& cur_node_indices = cur_node.task_info->rhs_indices;
        considered_indices -= cur_node_indices;
        if (considered_indices.none()) return;
    }
    DecisionBoundaryVector const& lhs_bounds = lattice_node_info.lhs_bounds;
    Index const next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
    Index const child_array_index = next_node_index - cur_node_index;
    OptionalChild const& optional_child = cur_node.children[child_array_index];
    if (!optional_child.has_value()) return;
    assert(next_node_index < lhs_bounds.size());
    BoundaryMap<Node> const& threshold_mapping = *optional_child;
    model::md::DecisionBoundary const next_lhs_bound = lhs_bounds[next_node_index];
    for (auto const& [threshold, node] : threshold_mapping) {
        assert(threshold > kLowestBound);
        if (threshold > next_lhs_bound) break;
        ExcludeGeneralizationRhs(node, lattice_node_info, next_node_index + 1, considered_indices);
        if (considered_indices.none()) return;
    }
}

void MinPickerLattice::RemoveSpecializations(Node& cur_node,
                                             MdLatticeNodeInfo const& lattice_node_info,
                                             model::Index cur_node_index,
                                             boost::dynamic_bitset<> const& picked_indices) {
    // All MDs in the tree are of the same cardinality.
    DecisionBoundaryVector const& lhs_bounds = lattice_node_info.lhs_bounds;
    model::Index const next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
    NodeChildren& children = cur_node.children;
    ValidationInfo*& task_info = cur_node.task_info;
    if (next_node_index == lhs_bounds.size()) {
        assert(children.empty());
        if (task_info != nullptr) {
            boost::dynamic_bitset<>& this_node_rhs_indices = task_info->rhs_indices;
            this_node_rhs_indices -= picked_indices;
            if (this_node_rhs_indices.none()) task_info = nullptr;
        }
        return;
    }
    model::Index const child_array_index = next_node_index - cur_node_index;
    OptionalChild& optional_child = children[child_array_index];
    if (!optional_child.has_value()) return;
    model::md::DecisionBoundary const next_node_bound = lhs_bounds[next_node_index];
    BoundaryMap<Node>& threshold_mapping = *optional_child;
    auto mapping_end = threshold_mapping.end();
    for (auto it_map = threshold_mapping.lower_bound(next_node_bound); it_map != mapping_end;
         ++it_map) {
        auto& node = it_map->second;
        RemoveSpecializations(node, lattice_node_info, next_node_index + 1, picked_indices);
    }
}

void MinPickerLattice::GetAll(Node& cur_node, std::vector<ValidationInfo>& collected,
                              model::Index cur_node_index) {
    ValidationInfo*& task_info = cur_node.task_info;
    if (task_info != nullptr) {
        collected.push_back(std::move(*task_info));
    }
    NodeChildren& children = cur_node.children;
    std::size_t const child_array_size = children.size();
    for (model::Index child_array_index = FindFirstNonEmptyIndex(children, 0);
         child_array_index != child_array_size;
         child_array_index = FindFirstNonEmptyIndex(children, child_array_index + 1)) {
        model::Index const next_node_index = cur_node_index + child_array_index;
        for (auto& [boundary, node] : *children[child_array_index]) {
            GetAll(node, collected, next_node_index + 1);
        }
    }
}

void MinPickerLattice::Add(ValidationInfo* validation_info) {
    DecisionBoundaryVector const& lhs_bounds = validation_info->node_info->lhs_bounds;
    std::size_t const col_match_number = lhs_bounds.size();

    Node* cur_node_ptr = &root_;
    for (Index cur_node_index = 0,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index);
         next_node_index != col_match_number; cur_node_index = next_node_index + 1,
               next_node_index = GetFirstNonZeroIndex(lhs_bounds, cur_node_index)) {
        model::md::DecisionBoundary const next_lhs_bound = lhs_bounds[next_node_index];
        model::Index const child_array_index = next_node_index - cur_node_index;
        std::size_t const next_child_array_size = col_match_number - next_node_index;
        auto [threshold_map, is_first_arr] =
                TryEmplaceChild(cur_node_ptr->children, child_array_index);
        if (is_first_arr) {
            Node& new_node =
                    threshold_map.try_emplace(next_lhs_bound, next_child_array_size).first->second;
            AddNewLhs(new_node, validation_info, next_node_index + 1);
            return;
        }
        auto [it_map, is_first_map] =
                threshold_map.try_emplace(next_lhs_bound, next_child_array_size);
        Node& next_node = it_map->second;
        if (is_first_map) {
            AddNewLhs(next_node, validation_info, next_node_index + 1);
            return;
        }
        cur_node_ptr = &next_node;
    }
    ValidationInfo*& task_info = cur_node_ptr->task_info;
    assert(task_info == nullptr);
    task_info = validation_info;
}

void MinPickerLattice::NewBatch(std::size_t max_elements) {
    info_.clear();
    std::size_t const vec_capacity = info_.capacity();
    if (max_elements > vec_capacity) {
        info_.reserve(std::max(vec_capacity * 2, max_elements));
    }
    root_ = Node{root_.children.size()};
}

void MinPickerLattice::AddGeneralizations(MdLatticeNodeInfo& lattice_node_info,
                                          boost::dynamic_bitset<>& considered_indices) {
    ExcludeGeneralizationRhs(root_, lattice_node_info, 0, considered_indices);
    if (considered_indices.none()) return;
    RemoveSpecializations(root_, lattice_node_info, 0, considered_indices);
    ValidationInfo& added_ref =
            info_.emplace_back(&lattice_node_info, std::move(considered_indices));
    Add(&added_ref);
}

std::vector<ValidationInfo> MinPickerLattice::GetAll() noexcept(kNeedsEmptyRemoval) {
    if constexpr (kNeedsEmptyRemoval) {
        // TODO: investigate different orders.
        return std::move(info_);

    } else {
        std::vector<ValidationInfo> collected;
        collected.reserve(info_.size());
        GetAll(root_, collected, 0);
        return collected;
    }
}

}  // namespace algos::hymd::lattice::cardinality
