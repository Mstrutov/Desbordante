#pragma once

#include <cassert>
#include <cstddef>
#include <map>
#include <optional>
#include <vector>

#include "algorithms/md/hymd/md_lhs.h"
#include "model/index.h"

namespace algos::hymd::lattice {
template <typename NodeType>
class NodeBase {
public:
    using BoundMap = std::map<model::md::DecisionBoundary, NodeType>;
    using OptionalChild = std::optional<BoundMap>;
    using Children = std::vector<OptionalChild>;

    Children children;

    NodeBase(std::size_t children_number) : children(children_number) {}

    std::size_t GetChildArraySize(model::Index child_array_index) const noexcept {
        return children.size() - (child_array_index + 1);
    }

    void ForEachNonEmpty(auto action) {
        ForEachNonEmpty(*this, action);
    }

    void ForEachNonEmpty(auto action) const {
        ForEachNonEmpty(*this, action);
    }

    template <typename Self>
    static void ForEachNonEmpty(Self& self, auto action) {
        model::Index index = 0;
        for (std::size_t const array_size = self.children.size(); index != array_size; ++index) {
            auto& b_map = self.children[index];
            if (b_map.has_value()) action(*b_map, index);
        };
    }

    template <typename... Args>
    std::pair<BoundMap&, bool> TryEmplaceChild(model::Index index, Args&&... args) {
        OptionalChild& optional_child = children[index];
        if (optional_child.has_value()) {
            return {*optional_child, false};
        }
        return {optional_child.emplace(std::forward<Args>(args)...), true};
    }

    bool IsEmpty() const noexcept {
        return std::none_of(children.begin(), children.end(),
                            std::mem_fn(&OptionalChild::has_value));
    }

protected:
    template <typename... Args>
    NodeType* AddOneUncheckedBase(model::Index child_array_index, model::md::DecisionBoundary bound,
                                  Args... args) {
        return &children[child_array_index]
                        .emplace()
                        .try_emplace(bound, args..., GetChildArraySize(child_array_index))
                        .first->second;
    }
};

template <typename NodeType>
void AddUnchecked(NodeType* cur_node_ptr, MdLhs const& lhs, model::Index cur_node_index,
                  auto final_node_action) {
    assert(cur_node_ptr->IsEmpty());
    for (MdElement element = lhs.FindNextNonZero(cur_node_index); lhs.IsNotEnd(element);
         element = lhs.FindNextNonZero(cur_node_index)) {
        auto const& [next_node_index, next_bound] = element;
        std::size_t const child_array_index = next_node_index - cur_node_index;
        cur_node_ptr = cur_node_ptr->AddOneUnchecked(child_array_index, next_bound);
        cur_node_index = next_node_index + 1;
    }
    final_node_action(cur_node_ptr);
}

template <typename NodeType>
void CheckedAdd(NodeType* cur_node_ptr, MdLhs const& lhs, auto const& info, auto unchecked_add,
                auto final_node_action) {
    model::Index cur_node_index = 0;
    for (MdElement element = lhs.FindNextNonZero(cur_node_index); lhs.IsNotEnd(element);
         element = lhs.FindNextNonZero(cur_node_index)) {
        auto const& [next_node_index, next_bound] = element;
        model::Index const child_array_index = next_node_index - cur_node_index;
        cur_node_index = next_node_index + 1;
        std::size_t const next_child_array_size =
                cur_node_ptr->GetChildArraySize(child_array_index);
        auto [boundary_map, is_first_arr] = cur_node_ptr->TryEmplaceChild(child_array_index);
        if (is_first_arr) {
            NodeType& new_node =
                    boundary_map.try_emplace(next_bound, next_child_array_size).first->second;
            unchecked_add(new_node, info, cur_node_index);
            return;
        }
        auto [it_map, is_first_map] = boundary_map.try_emplace(next_bound, next_child_array_size);
        NodeType& next_node = it_map->second;
        if (is_first_map) {
            unchecked_add(next_node, info, cur_node_index);
            return;
        }
        cur_node_ptr = &next_node;
    };
    final_node_action(cur_node_ptr);
}
}  // namespace algos::hymd::lattice
