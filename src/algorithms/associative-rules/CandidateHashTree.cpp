#include "CandidateHashTree.h"

void CandidateHashTree::appendRow(LeafRow const row, HashTreeNode & subtreeRoot) {
    if (!subtreeRoot.siblings.empty()) {
        auto const hash = hashFunction(row, subtreeRoot.levelNumber);
        appendRow(row, subtreeRoot.siblings[hash]);
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
    leafNode.siblings.reserve(branchingDegree);
    for (unsigned siblingNum = 0; siblingNum < branchingDegree; ++siblingNum) {
        //generate new leaf nodes
        leafNode.siblings.emplace_back(nextLevelNumber);
    }
    leafNode.siblings.shrink_to_fit();

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
    //auto currLevelElementID = static_cast<double>(nodeItems[levelNum - 1]);
    //double branchNum = currLevelElementID / static_cast<unsigned>(universeSize) * static_cast<unsigned>(branchingDegree);
    //return static_cast<unsigned>(branchNum);
    auto const currLevelElementID = nodeItems[levelNum - 1];
    return itemHash(currLevelElementID);
}

bool CandidateHashTree::findLeaves(std::set<unsigned>::iterator const start,
                                   HashTreeNode & subtreeRoot, std::set<unsigned> const& transactionItems,
                                   int transactionID) {
    if (subtreeRoot.siblings.empty() || start == transactionItems.end()) { //found a leaf node or final. //TODO start = end?
        if (subtreeRoot.lastVisitedTransactionID != transactionID) {
            traverseAllLeaves(subtreeRoot, transactionItems);
            subtreeRoot.lastVisitedTransactionID = transactionID;
        }
        return true;
    } else {
        unsigned const nextBranchNumber = itemHash(*start);
        auto & nextNode = subtreeRoot.siblings[nextBranchNumber];
        for (auto newStart = std::next(start); newStart != transactionItems.end(); ++newStart) {
            if (findLeaves(newStart, nextNode, transactionItems, transactionID)) {
                return false;
            }
        }
        findLeaves(transactionItems.end(), nextNode, transactionItems, transactionID);
        return false;
        //findLeaves(std::next(start), nextNode, transactionItems);
    }
}

void CandidateHashTree::traverseAllLeaves(HashTreeNode & subtreeRoot, std::set<unsigned> const& transactionItems) {
    //if (subtreeRoot.siblings.empty()) {
        for (auto & row : subtreeRoot.candidates) {
            auto const& candidateItems = row.nodeIter->items;
            bool isSubset = true;

            for (unsigned candidateItem : candidateItems) {
                if (transactionItems.find(candidateItem) == transactionItems.end()) {
                    isSubset = false;
                    break;
                }
            }
            if (isSubset) {
                row.transactionCount++;
            }
        }
    //} else {     //TODO хз нуэно или нет, скорее всего нет, потому что мы рассматриваем подмножества меньшей мощности
        //for (auto & sibling : subtreeRoot.siblings) {
        //    traverseAllLeaves(sibling, transactionItems);
        //}
    //}
}

void CandidateHashTree::performCounting() {
    for (auto const& transaction : transactionalData->getTransactions()) { //TODO нужно ли по порядку или так норм?
        auto const& items = transaction.second.getItemsIDs();
        auto const transactionID = transaction.first;
        for (auto start = items.begin(); start != items.end(); ++start) {
            if (findLeaves(start, root, items, transactionID)) {
                break;
            };
        }
        findLeaves(items.end(), root, items, transactionID);//TODO хз надо или нет
        //findLeaves(items.begin(), root, items);
    }
}

void CandidateHashTree::prune(double minsup, HashTreeNode & subtreeRoot) {
    if (subtreeRoot.siblings.empty()) {
        for (auto & row : subtreeRoot.candidates) {
            if (row.transactionCount < minsup * transactionalData->getNumTransactions()) {
                row.nodeContainer->erase(row.nodeIter);
            }
        }
    } else {
        for (auto & sibling : subtreeRoot.siblings) {
            prune(minsup, sibling);
        }
    }
}

void CandidateHashTree::pruneNodes(double minsup) {
    prune(minsup, root);
}
