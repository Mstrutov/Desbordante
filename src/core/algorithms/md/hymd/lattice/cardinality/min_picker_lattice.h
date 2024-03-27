#pragma once

#include <cstddef>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "algorithms/md/hymd/lattice/lattice_child_array.h"
#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "algorithms/md/hymd/lattice/validation_info.h"
#include "model/index.h"

namespace algos::hymd::lattice::cardinality {

class MinPickerLattice {
private:
    struct Node {
        LatticeChildArray<Node> children;
        ValidationInfo* task_info = nullptr;

        Node(std::size_t children_number) : children(children_number) {}
    };

    using BoundMap = BoundaryMap<Node>;
    using OptionalChild = std::optional<BoundMap>;
    using NodeChildren = LatticeChildArray<Node>;

    Node root_;
    std::vector<ValidationInfo> info_;

    void AddNewLhs(Node& cur_node, ValidationInfo* validation_info, model::Index cur_node_index);
    void ExcludeGeneralizationRhs(Node const& cur_node, MdLatticeNodeInfo const& lattice_node_info,
                                  model::Index cur_node_index,
                                  boost::dynamic_bitset<>& considered_indices);
    void RemoveSpecializations(Node& cur_node, MdLatticeNodeInfo const& lattice_node_info,
                               model::Index cur_node_index,
                               boost::dynamic_bitset<> const& picked_indices);
    void GetAll(Node& cur_node, std::vector<ValidationInfo>& collected,
                model::Index cur_node_index);
    void Add(ValidationInfo* validation_info);

public:
    static constexpr bool kNeedsEmptyRemoval = false;

    void NewBatch(std::size_t max_elements);
    void AddGeneralizations(MdLatticeNodeInfo& lattice_node_info,
                            boost::dynamic_bitset<>& considered_indices);
    std::vector<ValidationInfo> GetAll() noexcept(kNeedsEmptyRemoval);

    MinPickerLattice(std::size_t col_matches_num) : root_(col_matches_num) {}
};

}  // namespace algos::hymd::lattice::cardinality
