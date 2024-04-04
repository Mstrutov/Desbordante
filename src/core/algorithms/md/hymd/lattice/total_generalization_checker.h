#pragma once

#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/md_lhs.h"
#include "model/index.h"

namespace algos::hymd::lattice {
template <typename NodeType>
class TotalGeneralizationChecker {
    using Unspecialized = NodeType::Unspecialized;
    using OptionalChild = NodeType::OptionalChild;
    Unspecialized unspecialized_;

public:
    TotalGeneralizationChecker(Unspecialized unspecialized) noexcept
        : unspecialized_(unspecialized) {}

    bool HasGeneralizationInChildren(NodeType const& node, model::Index const node_index,
                                     model::Index const start_index) const {
        MdLhs const& lhs = NodeType::GetLhs(unspecialized_);
        for (MdElement element = lhs.FindNextNonZero(start_index); lhs.IsNotEnd(element);
             element = lhs.FindNextNonZero(element.index + 1)) {
            auto const& [next_node_index, generalization_bound_limit] = element;
            model::Index const child_array_index = next_node_index - node_index;
            OptionalChild const& optional_child = node.children[child_array_index];
            if (!optional_child.has_value()) continue;
            for (auto const& [generalization_bound, node] : *optional_child) {
                if (generalization_bound > generalization_bound_limit) break;
                if (HasGeneralization(node, next_node_index + 1)) return true;
            }
        }
        return false;
    }

    bool HasGeneralization(NodeType const& node, model::Index node_index) const {
        return node.ContainsGeneralizationOf(unspecialized_) ||
               HasGeneralizationInChildren(node, node_index, node_index);
    }

    bool HasGeneralization(NodeType const& node) const {
        return HasGeneralization(node, 0);
    }
};
}  // namespace algos::hymd::lattice
