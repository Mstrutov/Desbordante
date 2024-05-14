#pragma once

#include <cstddef>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "algorithms/md/hymd/lattice/md_lattice.h"
#include "algorithms/md/hymd/lattice/node_base.h"
#include "algorithms/md/hymd/lattice/validation_info.h"
#include "model/index.h"

namespace algos::hymd::lattice::cardinality {

class MinPickerLattice {
private:
    struct Node : NodeBase<Node> {
        ValidationInfo* task_info = nullptr;

        Node* AddOneUnchecked(model::Index child_array_index, model::md::DecisionBoundary bound) {
            return AddOneUncheckedBase(child_array_index, bound);
        }

        Node(std::size_t children_number) : NodeBase<Node>(children_number) {}
    };

    using BoundMap = Node::BoundMap;
    using OptionalChild = Node::OptionalChild;
    using NodeChildren = Node::Children;

    Node root_;
    std::vector<ValidationInfo> info_;

    static auto SetInfoAction(ValidationInfo* info) {
        return [info](Node* node) { node->task_info = info; };
    }
    void AddNewLhs(Node& cur_node, ValidationInfo* validation_info, MdLhs::iterator cur_lhs_iter);
    void ExcludeGeneralizationRhs(Node const& cur_node,
                                  MdLattice::MdVerificationMessenger const& messenger,
                                  MdLhs::iterator cur_lhs_iter,
                                  boost::dynamic_bitset<>& considered_indices);
    void RemoveSpecializations(Node& cur_node, MdLattice::MdVerificationMessenger const& messenger,
                               MdLhs::iterator cur_lhs_iter,
                               boost::dynamic_bitset<> const& picked_indices);
    void GetAll(Node& cur_node, std::vector<ValidationInfo>& collected);
    void Add(ValidationInfo* validation_info);

public:
    static constexpr bool kNeedsEmptyRemoval = false;

    void NewBatch(std::size_t max_elements);
    void AddGeneralizations(MdLattice::MdVerificationMessenger& messenger,
                            boost::dynamic_bitset<>& considered_indices);
    std::vector<ValidationInfo> GetAll() noexcept(kNeedsEmptyRemoval);

    MinPickerLattice(std::size_t col_matches_num) : root_(col_matches_num) {}
};

}  // namespace algos::hymd::lattice::cardinality
