#pragma once

#include "TransactionalData.h"
#include "Node.h"

class TransactionalData;

class CandidateHashTree {
    using NodeIterator = std::list<Node>::iterator;
private:
    unsigned branchingDegree;
    unsigned minThreshold;
    unsigned totalRowCount = 0;

    TransactionalData const* const transactionalData = nullptr;

    struct LeafRow {
        NodeIterator nodeIter;
        std::list<Node>* const nodeContainer = nullptr;
        unsigned transactionCount = 0;

        LeafRow() = default;
        //TODO копирование?
        LeafRow(NodeIterator nodeIter, std::list<Node>* const nodeContainer)
            : nodeIter(nodeIter), nodeContainer(nodeContainer) {}
    };

    struct HashTreeNode {
        unsigned levelNumber;
        std::vector<HashTreeNode> siblings;
        std::list<LeafRow> candidates;

        explicit HashTreeNode(unsigned levelNumber)
            : levelNumber(levelNumber) {}
    };

    HashTreeNode root;
    unsigned hashFunction(LeafRow const& nodeRow, unsigned levelNum) const;
    unsigned itemHash(unsigned itemID) const noexcept { return itemID % branchingDegree; }

    void appendRow(LeafRow row, HashTreeNode & subtreeRoot);
    void addLevel(HashTreeNode & leafNode);
    void findLeaves(std::set<unsigned>::iterator start,
                    HashTreeNode & subtreeRoot, std::set<unsigned> const& transactionItems);
    static void traverseAllLeaves(HashTreeNode & subtreeRoot, std::set<unsigned> const& transactionItems);
    void prune(double minsup, HashTreeNode & subtreeRoot);
public:
    CandidateHashTree(TransactionalData const* transactionalData,
                      unsigned branchingDegree, unsigned minThreshold)
            : branchingDegree(branchingDegree), minThreshold(minThreshold),
              transactionalData(transactionalData), root(1) {}

    void addCandidate(NodeIterator nodeIter, std::list<Node>* nodeContainer); //maybe reference is possible
    unsigned size() const noexcept { return totalRowCount; };
    void performCounting();
    void pruneNodes(double minsup);
};
