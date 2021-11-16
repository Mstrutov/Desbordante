#include "LatticeObservations.h"

NodeCategory LatticeObservations::updateDependencyCategory(Vertical const& node) {
    NodeCategory newCategory;
    if (node.getArity() <= 1) {
        newCategory = NodeCategory::minimalDependency;
        (*this)[node] = newCategory;
        return newCategory;
    }

    auto columnIndices = node.getColumnIndicesRef(); //copy indices
    bool hasUncheckedSubset = false;

    for (size_t index = columnIndices.find_first();
         index < columnIndices.size();
         index = columnIndices.find_next(index)
    ) {
        columnIndices[index] = false; //remove one column
        auto const subsetNodeIter = this->find(Vertical(node.getSchema(), columnIndices));

        if (subsetNodeIter == this->end()) {
            //if we found unchecked subset of this node
            hasUncheckedSubset = true;
        } else {
            NodeCategory const &subsetVerticalCategory = subsetNodeIter->second;
            if (subsetVerticalCategory == NodeCategory::minimalDependency ||
                subsetVerticalCategory == NodeCategory::dependency ||
                subsetVerticalCategory == NodeCategory::candidateMinimalDependency
            ) {
                newCategory = NodeCategory::dependency;
                (*this)[node] = newCategory;
                return newCategory;
            }
        }

        columnIndices[index] = true; //restore removed column
    }
    newCategory = hasUncheckedSubset ? NodeCategory::candidateMinimalDependency : NodeCategory::minimalDependency;
    (*this)[node] = newCategory;
    return newCategory;
}

NodeCategory LatticeObservations::updateNonDependencyCategory(Vertical const& node, unsigned int rhsIndex) {
    auto columnIndices = node.getColumnIndicesRef();
    columnIndices[rhsIndex] = true;
    columnIndices.flip();

    NodeCategory newCategory;
    bool hasUncheckedSuperset = false;

    for (size_t index = columnIndices.find_first();
         index < columnIndices.size();
         index = columnIndices.find_next(index)
    ) {
        auto const supersetNodeIter = this->find(node.Union(*node.getSchema()->getColumn(index)));

        if (supersetNodeIter == this->end()) {
            //if we found unchecked superset of this node
            hasUncheckedSuperset = true;
        } else {
            NodeCategory const &supersetVerticalCategory = supersetNodeIter->second;
            if (supersetVerticalCategory == NodeCategory::maximalNonDependency ||
                supersetVerticalCategory == NodeCategory::nonDependency ||
                supersetVerticalCategory == NodeCategory::candidateMaximalNonDependency
            ) {
                newCategory = NodeCategory::nonDependency;
                (*this)[node] = newCategory;
                return newCategory;
            }
        }
    }
    newCategory = hasUncheckedSuperset ? NodeCategory::candidateMaximalNonDependency : NodeCategory::maximalNonDependency;
    (*this)[node] = newCategory;
    return newCategory;
}

bool LatticeObservations::isCandidate(Vertical const& node) const {
    auto nodeIter = this->find(node);
    if (nodeIter == this->end()) {
        return false;
    } else {
        return nodeIter->second == NodeCategory::candidateMaximalNonDependency ||
               nodeIter->second == NodeCategory::candidateMinimalDependency;
    }
}

std::unordered_set<Vertical>
LatticeObservations::getUncheckedSubsets(Vertical const& node, ColumnOrder const& columnOrder) const {
    auto indices = node.getColumnIndices();
    std::unordered_set<Vertical> uncheckedSubsets;

    for (int columnIndex : columnOrder.getOrderHighDistinctCount(node)) {
        indices[columnIndex] = false;
        Vertical subsetNode = Vertical(node.getSchema(), indices);
        if (this->find(subsetNode) == this->end()) {
            uncheckedSubsets.insert(std::move(subsetNode));
        }
        indices[columnIndex] = true;
    }

    return uncheckedSubsets;
}

std::unordered_set<Vertical>
LatticeObservations::getUncheckedSupersets(Vertical const& node, unsigned int rhsIndex, ColumnOrder const& columnOrder) const {
    auto flippedIndices = node.getColumnIndices().flip();
    std::unordered_set<Vertical> uncheckedSupersets;

    flippedIndices[rhsIndex] = false;

    for (int columnIndex : columnOrder.getOrderHighDistinctCount(Vertical(node.getSchema(), flippedIndices))) {
        auto indices = node.getColumnIndices();

        indices[columnIndex] = true;
        Vertical subsetNode = Vertical(node.getSchema(), indices);
        if (this->find(subsetNode) == this->end()) {
            uncheckedSupersets.insert(std::move(subsetNode));
        }
    }

    return uncheckedSupersets;
}
