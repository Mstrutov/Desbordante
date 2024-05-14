#pragma once

#include "algorithms/md/hymd/lattice/lhs_specialization.h"
#include "algorithms/md/hymd/lattice/total_generalization_checker.h"

namespace algos::hymd::lattice {
template <typename NodeType>
class SpecGeneralizationChecker {
    using Specialization = NodeType::Specialization;
    using BoundMap = NodeType::BoundMap;
    using OptionalChild = NodeType::OptionalChild;

    Specialization const& specialization_;
    TotalGeneralizationChecker<NodeType> total_checker_{specialization_.ToUnspecialized()};

    bool HasGeneralizationTotal(NodeType const& node, MdLhs::iterator iter,
                                model::Index child_array_index) const {
        return total_checker_.HasGeneralization(node, iter, child_array_index);
    }

    bool HasChildGenSpec(NodeType const& node, model::Index child_array_index,
                         MdLhs::iterator fol_iter, model::md::DecisionBoundary bound_limit,
                         model::Index next_child_array_index, auto gen_method,
                         auto get_b_map_iter) const {
        OptionalChild const& optional_child = node.children[child_array_index];
        if (!optional_child.has_value()) return false;
        BoundMap const& b_map = *optional_child;
        for (auto spec_iter = get_b_map_iter(b_map), end_iter = b_map.end(); spec_iter != end_iter;
             ++spec_iter) {
            auto const& [generalization_bound, node] = *spec_iter;
            if (generalization_bound > bound_limit) break;
            if ((this->*gen_method)(node, fol_iter, next_child_array_index)) return true;
        }
        return false;
    }

public:
    SpecGeneralizationChecker(Specialization const& specialization)
        : specialization_(specialization) {}

    bool HasGeneralizationInChildren(NodeType const& node, MdLhs::iterator next_node_iter,
                                     model::Index child_array_index = 0) const {
        LhsSpecialization const& lhs_specialization = specialization_.GetLhsSpecialization();
        MdLhs const& old_lhs = lhs_specialization.old_lhs;
        auto spec_iter = lhs_specialization.specialization_data.spec_before;
        auto get_first = [](BoundMap const& b_map) { return b_map.begin(); };
        while (next_node_iter != spec_iter) {
            auto const& [delta, next_bound] = *next_node_iter;
            ++next_node_iter;
            child_array_index += delta;
            if (HasChildGenSpec(node, child_array_index, next_node_iter, next_bound, 0,
                                &SpecGeneralizationChecker::HasGeneralizationInChildren, get_first))
                return true;
            ++child_array_index;
        }
        auto const& [spec_delta, spec_bound] = lhs_specialization.specialization_data.new_child;
        child_array_index += spec_delta;
        if (spec_iter != old_lhs.end() && spec_iter->child_array_index == spec_delta) {
            model::md::DecisionBoundary const old_bound = spec_iter->decision_boundary;
            auto get_higher = [&](BoundMap const& b_map) { return b_map.upper_bound(old_bound); };
            return HasChildGenSpec(node, child_array_index, spec_iter + 1, spec_bound, 0,
                                   &SpecGeneralizationChecker::HasGeneralizationTotal, get_higher);
        } else {
            assert(spec_iter == old_lhs.end() || spec_iter->child_array_index > spec_delta);
            return HasChildGenSpec(node, child_array_index, spec_iter, spec_bound,
                                   -(spec_delta + 1),
                                   &SpecGeneralizationChecker::HasGeneralizationTotal, get_first);
        }
    }

    bool HasGeneralization(NodeType const& node, MdLhs::iterator iter) const {
        return HasGeneralizationInChildren(node, iter);
    }

    bool HasGeneralization(NodeType const& node) const {
        return HasGeneralization(node, specialization_.GetLhsSpecialization().old_lhs.begin());
    }

    auto const& GetTotalChecker() const noexcept {
        return total_checker_;
    }
};
}  // namespace algos::hymd::lattice
