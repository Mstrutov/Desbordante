#pragma once

#include "algorithms/md/hymd/lattice/lhs_specialization.h"
#include "algorithms/md/hymd/lattice/total_generalization_checker.h"

namespace algos::hymd::lattice {
template <typename NodeType>
class SpecGeneralizationChecker {
    using Specialization = NodeType::Specialization;
    Specialization const& specialization_;
    TotalGeneralizationChecker<NodeType> total_checker_{specialization_.ToUnspecialized()};

    bool HasGeneralizationTotal(NodeType const& node, model::Index node_index) const {
        return total_checker_.HasGeneralization(node, node_index);
    }

    bool HasChildGenSpec(NodeType const& node, model::Index node_index,
                         model::Index next_node_index, model::md::DecisionBoundary bound_limit,
                         auto gen_method, auto get_b_map_iter) const {
        model::Index const child_array_index = next_node_index - node_index;
        typename NodeType::OptionalChild const& optional_child = node.children[child_array_index];
        if (!optional_child.has_value()) return false;
        typename NodeType::BoundMap const& b_map = *optional_child;
        for (auto spec_iter = get_b_map_iter(b_map), end_iter = b_map.end(); spec_iter != end_iter;
             ++spec_iter) {
            auto const& [generalization_bound, node] = *spec_iter;
            if (generalization_bound > bound_limit) break;
            model::Index const fol_index = next_node_index + 1;
            if ((this->*gen_method)(node, fol_index)) return true;
        }
        return false;
    }

public:
    SpecGeneralizationChecker(Specialization const& specialization)
        : specialization_(specialization) {}

    bool HasGeneralizationInChildren(NodeType const& node, model::Index node_index,
                                     model::Index start_index) const {
        LhsSpecialization const& lhs_specialization = specialization_.GetLhsSpecialization();
        auto const& [spec_index, spec_bound] = lhs_specialization.specialized;
        MdLhs const& old_lhs = lhs_specialization.old_lhs;
        using BoundMap = NodeType::BoundMap;
        for (MdElement result = old_lhs.FindNextNonZero(start_index); result.index < spec_index;
             result = old_lhs.FindNextNonZero(result.index + 1)) {
            auto const& [next_node_index, next_bound] = result;
            auto get_first = [](BoundMap const& b_map) { return b_map.begin(); };
            if (HasChildGenSpec(node, node_index, next_node_index, next_bound,
                                &SpecGeneralizationChecker::HasGeneralization, get_first))
                return true;
        }
        model::md::DecisionBoundary const old_bound = old_lhs[spec_index];
        auto get_higher = [&](BoundMap const& b_map) { return b_map.upper_bound(old_bound); };
        return HasChildGenSpec(node, node_index, spec_index, spec_bound,
                               &SpecGeneralizationChecker::HasGeneralizationTotal, get_higher);
    }

    bool HasGeneralization(NodeType const& node, model::Index node_index) const {
        return HasGeneralizationInChildren(node, node_index, node_index);
    }

    auto const& GetTotalChecker() const noexcept {
        return total_checker_;
    }
};
}  // namespace algos::hymd::lattice
