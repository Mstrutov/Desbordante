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

    static auto GetFirstAction() noexcept {
        return [](BoundMap const& b_map) { return b_map.begin(); };
    }

    bool HasGeneralizationInChildren(NodeType const& node, MdLhs::iterator next_node_iter,
                                     model::Index child_array_index, auto spec_method,
                                     auto final_method) const {
        LhsSpecialization const& lhs_specialization = specialization_.GetLhsSpecialization();
        MdLhs::iterator spec_iter = lhs_specialization.specialization_data.spec_before;
        while (next_node_iter != spec_iter) {
            auto const& [delta, next_bound] = *next_node_iter;
            ++next_node_iter;
            child_array_index += delta;
            if (HasChildGenSpec(node, child_array_index, next_node_iter, next_bound, 0, spec_method,
                                GetFirstAction()))
                return true;
            ++child_array_index;
        }
        auto const& [spec_delta, spec_bound] = lhs_specialization.specialization_data.new_child;
        child_array_index += spec_delta;
        return (this->*final_method)(node, child_array_index, spec_iter, spec_bound, spec_delta);
    }

    bool ReplaceFinalCheck(NodeType const& node, model::Index child_array_index,
                           MdLhs::iterator spec_iter, model::md::DecisionBoundary spec_bound,
                           std::size_t spec_delta) const {
        assert(spec_iter != specialization_.GetLhsSpecialization().old_lhs.end() &&
               spec_iter->child_array_index == spec_delta);
        model::md::DecisionBoundary const old_bound = spec_iter->decision_boundary;
        auto get_higher = [&](BoundMap const& b_map) { return b_map.upper_bound(old_bound); };
        ++spec_iter;
        return HasChildGenSpec(node, child_array_index, spec_iter, spec_bound, 0,
                               &SpecGeneralizationChecker::HasGeneralizationTotal, get_higher);
    }

    bool NonReplaceFinalCheck(NodeType const& node, model::Index child_array_index,
                              MdLhs::iterator spec_iter, model::md::DecisionBoundary spec_bound,
                              std::size_t spec_delta) const {
        assert(spec_iter == specialization_.GetLhsSpecialization().old_lhs.end() ||
               spec_iter->child_array_index > spec_delta);
        return HasChildGenSpec(node, child_array_index, spec_iter, spec_bound, -(spec_delta + 1),
                               &SpecGeneralizationChecker::HasGeneralizationTotal,
                               GetFirstAction());
    }

public:
    SpecGeneralizationChecker(Specialization const& specialization)
        : specialization_(specialization) {}

    bool HasGeneralizationInChildrenReplace(NodeType const& node, MdLhs::iterator next_node_iter,
                                            model::Index child_array_index = 0) const {
        return HasGeneralizationInChildren(
                node, next_node_iter, child_array_index,
                &SpecGeneralizationChecker::HasGeneralizationInChildrenReplace,
                &SpecGeneralizationChecker::ReplaceFinalCheck);
    }

    bool HasGeneralizationInChildrenNonReplace(NodeType const& node, MdLhs::iterator next_node_iter,
                                               model::Index child_array_index = 0) const {
        return HasGeneralizationInChildren(
                node, next_node_iter, child_array_index,
                &SpecGeneralizationChecker::HasGeneralizationInChildrenNonReplace,
                &SpecGeneralizationChecker::NonReplaceFinalCheck);
    }

    bool HasGeneralizationReplace(NodeType const& node, MdLhs::iterator iter) const {
        return HasGeneralizationInChildrenReplace(node, iter);
    }

    bool HasGeneralizationNonReplace(NodeType const& node, MdLhs::iterator iter) const {
        return HasGeneralizationInChildrenNonReplace(node, iter);
    }

    bool HasGeneralizationReplace(NodeType const& node) const {
        return HasGeneralizationReplace(node,
                                        specialization_.GetLhsSpecialization().old_lhs.begin());
    }

    bool HasGeneralizationNonReplace(NodeType const& node) const {
        return HasGeneralizationNonReplace(node,
                                           specialization_.GetLhsSpecialization().old_lhs.begin());
    }

    auto const& GetTotalChecker() const noexcept {
        return total_checker_;
    }
};
}  // namespace algos::hymd::lattice
