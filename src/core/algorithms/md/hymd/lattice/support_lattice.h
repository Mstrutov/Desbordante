#pragma once

#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lattice/lattice_child_array.h"
#include "model/index.h"

namespace algos::hymd::lattice {

class SupportLattice {
private:
    struct Node {
        LatticeChildArray<Node> children;
        bool is_unsupported = false;

        Node(std::size_t children_number) : children(children_number) {}
    };

    using NodeChildren = LatticeChildArray<Node>;

    Node root_;

    bool IsUnsupported(Node const& cur_node, DecisionBoundaryVector const& lhs_bounds,
                       model::Index cur_node_index) const;
    void MarkNewLhs(Node& cur_node, DecisionBoundaryVector const& lhs_bounds,
                    model::Index cur_node_index);

public:
    void MarkUnsupported(DecisionBoundaryVector const& lhs_bounds);
    bool IsUnsupported(DecisionBoundaryVector const& lhs_bounds) const;

    SupportLattice(std::size_t column_matches_size) : root_(column_matches_size) {}
};

}  // namespace algos::hymd::lattice
