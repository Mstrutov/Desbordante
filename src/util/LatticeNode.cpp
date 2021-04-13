//
// Created by alex on 13.04.2021.
//

#include "LatticeNode.h"

LatticeNode::LatticeNode(Column const& column)
    : Vertical(column) {
    category = NodeCategory::unvisited;
    partition = nullptr; //?
}

bool LatticeNode::isCandidate() {
    return category == NodeCategory::candidateMaximalNonDependency ||
           category == NodeCategory::candidateMinimalDependency;
}

bool LatticeNode::isMinimalDependency() {
    for ()
}
