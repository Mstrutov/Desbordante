#pragma once

#include "algorithms/md/hymd/md_element.h"
#include "algorithms/md/hymd/md_lhs.h"
#include "model/index.h"

namespace algos::hymd::lattice {
template <typename NodeType>
class TotalGeneralizationChecker {
    using Unspecialized = NodeType::Specialization::Unspecialized;
    using OptionalChild = NodeType::OptionalChild;
    Unspecialized unspecialized_;

public:
    TotalGeneralizationChecker(Unspecialized unspecialized) noexcept
        : unspecialized_(unspecialized) {}

    bool HasGeneralizationInChildren(NodeType const& node, MdLhs::iterator next_iter,
                                     model::Index child_array_index = 0) const {
        MdLhs const& lhs = NodeType::GetLhs(unspecialized_);
        for (MdLhs::iterator end_iter = lhs.end(); next_iter != end_iter; ++child_array_index) {
            auto const& [index_delta, generalization_bound_limit] = *next_iter;
            child_array_index += index_delta;
            ++next_iter;
            OptionalChild const& optional_child = node.children[child_array_index];
            if (!optional_child.has_value()) continue;
            for (auto const& [generalization_bound, node] : *optional_child) {
                if (generalization_bound > generalization_bound_limit) break;
                if (HasGeneralization(node, next_iter)) return true;
            }
        }
        return false;
    }

    bool HasGeneralization(NodeType const& node, MdLhs::iterator next_iter,
                           model::Index child_array_index = 0) const {
        return node.ContainsGeneralizationOf(unspecialized_) ||
               HasGeneralizationInChildren(node, next_iter, child_array_index);
    }

    bool HasGeneralization(NodeType const& node) const {
        return HasGeneralization(node, NodeType::GetLhs(unspecialized_).begin());
    }
};
}  // namespace algos::hymd::lattice
