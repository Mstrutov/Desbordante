#include "CandidateHashTree.h"

#include <algorithm>

void CandidateHashTree::appendRow(LeafRow const row, HashTreeNode & subtreeRoot) {
    if (!subtreeRoot.children.empty()) {
        auto const hash = hashFunction(row, subtreeRoot.levelNumber);
        appendRow(row, subtreeRoot.children[hash]);
    } else {
        subtreeRoot.candidates.push_back(row);
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
        appendRow(row, leafNode); //TODO как копируем?
    }

    leafNode.candidates.clear();
}

void CandidateHashTree::addCandidate(NodeIterator nodeIter, std::list<Node> *const nodeContainer) {
    appendRow(LeafRow(nodeIter, nodeContainer), root);
    ++totalRowCount;
}

unsigned CandidateHashTree::hashFunction(LeafRow const& nodeRow, unsigned const levelNum) const {
    auto const& nodeItems = nodeRow.nodeIter->items;
    auto const currLevelElementID = nodeItems[levelNum - 1];
    return itemHash(currLevelElementID);
}

bool CandidateHashTree::findAndVisitLeaf(HashTreeNode & subtreeRoot,
                                         std::vector<unsigned>::const_iterator const start,
                                         std::vector<unsigned> const& transactionItems,
                                         int transactionID) {
    if (subtreeRoot.children.empty() || start == transactionItems.end()) { //found a leaf node or final. //TODO start = end?
        visitLeaf(subtreeRoot, transactionItems, transactionID);
        return true;
    } else {
        unsigned const nextBranchNumber = itemHash(*start);
        auto & nextNode = subtreeRoot.children[nextBranchNumber];

        /* Initially, the cycle started from the newStart = next(start) and passed newStart to findLeaf().
         * But if 'start' points to the latest element, next(start) will point to end(),
         * and cycle will even not enter its body. So in that case the function findAndVisitLeaf()
         * will not be executed for the just defined nextNode. The solution is to start
         * the cycle from 'start' and pass next(start) to the function.
         */
        for (auto newStart = start; newStart != transactionItems.end(); ++newStart) {
            if (findAndVisitLeaf(nextNode, std::next(newStart), transactionItems, transactionID)) {
                //this means nextNode is a leaf that we visited, so we can avoid traversing it in cycle again
                break;
            }
        }
        return false;
    }
}

void CandidateHashTree::visitLeaf(HashTreeNode & leaf, std::vector<unsigned> const& transactionItems, int tID) {
    if (leaf.lastVisitedTransactionID == tID) {
        return;
    }
    leaf.lastVisitedTransactionID = tID;

    for (auto & row : leaf.candidates) {
        auto const& candidateItems = row.nodeIter->items;
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

        for (auto start = items.begin(); start != items.end(); ++start) {
            if (findAndVisitLeaf(root, start, items, transactionID)) {
                break;
            };
        }
    }
}

void CandidateHashTree::prune(double minsup, HashTreeNode & subtreeRoot) {
    if (subtreeRoot.children.empty()) {
        for (auto & row : subtreeRoot.candidates) {
            double const support = static_cast<double>(row.transactionCount) / transactionalData->getNumTransactions();
            if (support < minsup) {
                row.nodeContainer->erase(row.nodeIter);
            } else {
                row.nodeIter->support = support;
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
