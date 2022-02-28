#include "CandidateHashTree.h"

#include <algorithm>
#include "cassert"

void CandidateHashTree::appendRow(LeafRow row, HashTreeNode & subtreeRoot) {
    if (!subtreeRoot.children.empty()) {
        auto const hash = hashFunction(row, subtreeRoot.levelNumber);
        appendRow(std::move(row), subtreeRoot.children[hash]);
    } else {
        subtreeRoot.candidates.push_back(std::move(row));
        if (subtreeRoot.candidates.size() > minThreshold) {
            addLevel(subtreeRoot);
        }
    }
}

void CandidateHashTree::addLevel(HashTreeNode & leafNode) {
    unsigned const nextLevelNumber = leafNode.levelNumber + 1;

    //by this we make leaf node interior
    leafNode.children.reserve(branchingDegree);
    for (unsigned siblingNum = 0; siblingNum < branchingDegree; ++siblingNum) {
        //generate new leaf nodes
        leafNode.children.emplace_back(nextLevelNumber);
    }
    leafNode.children.shrink_to_fit();

    //distribute rows of an old leaf between new leaves
    for (auto & row : leafNode.candidates) {
        appendRow(std::move(row), leafNode);
    }

    leafNode.candidates.clear();
}

void CandidateHashTree::addCandidate(NodeIterator candidate, Node* parent) {
    appendRow(LeafRow(candidate, parent), root);
    ++totalRowCount;
}

unsigned CandidateHashTree::hashFunction(LeafRow const& nodeRow, unsigned const levelNum) const {
    auto const& nodeItems = nodeRow.candidateNode->items;
    assert(levelNum <= nodeItems.size());
    auto const currLevelElementID = nodeItems.at(levelNum - 1);
    return itemHash(currLevelElementID);
}

void CandidateHashTree::findAndVisitLeaves(HashTreeNode & subtreeRoot,
                                           std::vector<unsigned>::const_iterator const start,
                                           std::vector<unsigned> const& transactionItems,
                                           int transactionID) {
    unsigned const nextBranchNumber = itemHash(*start);
    auto & nextNode = subtreeRoot.children[nextBranchNumber];
    if (nextNode.children.empty()) {
        //if nextNode is a leaf itself, then we visit it and terminate the recursion
        visitLeaf(nextNode, transactionItems, transactionID);
    } else {
        for (auto newStart = std::next(start); newStart != transactionItems.end(); ++newStart) {
            findAndVisitLeaves(nextNode, newStart, transactionItems, transactionID);
        }
    }
}

void CandidateHashTree::visitLeaf(HashTreeNode & leaf, std::vector<unsigned> const& transactionItems, int tID) {
    if (leaf.lastVisitedTransactionID == tID) {
        return;
    }
    leaf.lastVisitedTransactionID = tID;

    for (auto & row : leaf.candidates) {
        auto const& candidateItems = row.candidateNode->items;
        if (std::includes(transactionItems.begin(), transactionItems.end(),
                          candidateItems.begin(), candidateItems.end())) {
            row.transactionCount++;
        }
    }
}

void CandidateHashTree::performCounting() {
    for (auto const& transaction : transactionalData->getTransactions()) { //TODO нужно ли по порядку или так норм?
        auto const& items = transaction.second.getItemsIDs();
        auto const transactionID = transaction.first;
        if (root.children.empty()) {
            //if the root is a leaf itself
            visitLeaf(root, items, transactionID);
        } else {
            for (auto start = items.begin(); start != items.end(); ++start) {
                findAndVisitLeaves(root, start, items, transactionID);
            }
        }
    }
}

void CandidateHashTree::prune(double minsup, HashTreeNode & subtreeRoot) {
    if (subtreeRoot.children.empty()) {
        for (auto & row : subtreeRoot.candidates) {
            double const support = static_cast<double>(row.transactionCount) / transactionalData->getNumTransactions();
            if (support < minsup) {
                candidates[row.parent].erase(row.candidateNode);
            } else {
                row.candidateNode->support = support;
            }
        }
    } else {
        for (auto & sibling : subtreeRoot.children) {
            prune(minsup, sibling);
        }
    }
}

void CandidateHashTree::pruneNodes(double minsup) {
    prune(minsup, root);
}

void CandidateHashTree::addCandidates() {
    for (auto& [node, candidateChildren] : candidates) {
        for (auto child = candidateChildren.begin(); child != candidateChildren.end(); ++child) {
            addCandidate(child, node);
        }
    }
}
