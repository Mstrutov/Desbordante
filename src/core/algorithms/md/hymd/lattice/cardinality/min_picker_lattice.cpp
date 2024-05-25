#include "algorithms/md/hymd/lattice/cardinality/min_picker_lattice.h"

#include <cassert>

#include "algorithms/md/hymd/lowest_bound.h"

namespace {
using namespace algos::hymd;
using model::Index;
}  // namespace

namespace algos::hymd::lattice::cardinality {

void MinPickerLattice::AddNewLhs(Node& cur_node, ValidationInfo* validation_info,
                                 MdLhs::iterator cur_lhs_iter) {
    auto set_info = SetInfoAction(validation_info);
    AddUnchecked(&cur_node, validation_info->messenger->GetLhs(), cur_lhs_iter, set_info);
}

void MinPickerLattice::ExcludeGeneralizationRhs(Node const& cur_node,
                                                MdLattice::MdVerificationMessenger const& messenger,
                                                MdLhs::iterator cur_lhs_iter,
                                                boost::dynamic_bitset<>& considered_indices) {
    if (cur_node.task_info != nullptr) {
        boost::dynamic_bitset<> const& cur_node_indices = cur_node.task_info->rhs_indices;
        considered_indices -= cur_node_indices;
        return;
    }
    auto const& [child_array_index, next_lhs_bound] = *cur_lhs_iter;
    ++cur_lhs_iter;
    OptionalChild const& optional_child = cur_node.children[child_array_index];
    if (!optional_child.has_value()) return;
    Node::BoundMap const& bound_map = *optional_child;
    for (auto const& [bound, node] : bound_map) {
        if (bound > next_lhs_bound) break;
        ExcludeGeneralizationRhs(node, messenger, cur_lhs_iter, considered_indices);
        if (considered_indices.none()) return;
    }
}

void MinPickerLattice::RemoveSpecializations(Node& cur_node,
                                             MdLattice::MdVerificationMessenger const& messenger,
                                             MdLhs::iterator cur_lhs_iter,
                                             boost::dynamic_bitset<> const& picked_indices) {
    // All MDs in the tree are of the same cardinality.
    MdLhs const& lhs = messenger.GetLhs();
    ValidationInfo*& task_info = cur_node.task_info;
    if (cur_lhs_iter == lhs.end()) {
        assert(cur_node.IsEmpty());
        if (task_info != nullptr) {
            boost::dynamic_bitset<>& this_node_rhs_indices = task_info->rhs_indices;
            this_node_rhs_indices -= picked_indices;
            if (this_node_rhs_indices.none()) task_info = nullptr;
        }
        return;
    }
    auto const& [child_array_index, next_node_bound] = *cur_lhs_iter;
    ++cur_lhs_iter;
    OptionalChild& optional_child = cur_node.children[child_array_index];
    if (!optional_child.has_value()) return;
    Node::BoundMap& bound_map = *optional_child;
    auto mapping_end = bound_map.end();
    for (auto it_map = bound_map.lower_bound(next_node_bound); it_map != mapping_end; ++it_map) {
        auto& node = it_map->second;
        RemoveSpecializations(node, messenger, cur_lhs_iter, picked_indices);
    }
}

void MinPickerLattice::GetAll(Node& cur_node, std::vector<ValidationInfo>& collected) {
    ValidationInfo*& task_info = cur_node.task_info;
    if (task_info != nullptr) {
        collected.push_back(std::move(*task_info));
    }
    auto collect = [&](BoundMap& bound_map, model::Index) {
        for (auto& [boundary, node] : bound_map) {
            GetAll(node, collected);
        }
    };
    cur_node.ForEachNonEmpty(collect);
}

void MinPickerLattice::Add(ValidationInfo* validation_info) {
    auto add_new = [this](auto&&... args) { AddNewLhs(std::forward<decltype(args)>(args)...); };
    auto set_info = SetInfoAction(validation_info);
    CheckedAdd(&root_, validation_info->messenger->GetLhs(), validation_info, add_new, set_info);
}

void MinPickerLattice::NewBatch(std::size_t max_elements) {
    info_.clear();
    std::size_t const vec_capacity = info_.capacity();
    if (max_elements > vec_capacity) {
        info_.reserve(std::max(vec_capacity * 2, max_elements));
    }
    root_ = {root_.children.size()};
}

void MinPickerLattice::AddGeneralizations(MdLattice::MdVerificationMessenger& messenger,
                                          boost::dynamic_bitset<>& considered_indices) {
    MdLhs const& lhs = messenger.GetLhs();
    if (!lhs.IsEmpty()) [[likely]] {
        ExcludeGeneralizationRhs(root_, messenger, lhs.begin(), considered_indices);
        if (considered_indices.none()) return;
    }
    RemoveSpecializations(root_, messenger, lhs.begin(), considered_indices);
    ValidationInfo& added_ref = info_.emplace_back(&messenger, std::move(considered_indices));
    Add(&added_ref);
}

std::vector<ValidationInfo> MinPickerLattice::GetAll() noexcept(kNeedsEmptyRemoval) {
    if constexpr (kNeedsEmptyRemoval) {
        // TODO: investigate different orders.
        return std::move(info_);
    } else {
        std::vector<ValidationInfo> collected;
        collected.reserve(info_.size());
        GetAll(root_, collected);
        return collected;
    }
}

}  // namespace algos::hymd::lattice::cardinality
