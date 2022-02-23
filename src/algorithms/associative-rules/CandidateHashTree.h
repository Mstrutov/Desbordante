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
        int lastVisitedTransactionID = -1; //TODO ниче что не инициализируем?
        std::vector<HashTreeNode> children;
        std::list<LeafRow> candidates;

        HashTreeNode() = delete;
        explicit HashTreeNode(unsigned levelNumber)
            : levelNumber(levelNumber) {}
    };

    HashTreeNode root;
    unsigned hashFunction(LeafRow const& nodeRow, unsigned levelNum) const;
    unsigned itemHash(unsigned itemID) const noexcept { return itemID % branchingDegree; }

    void appendRow(LeafRow row, HashTreeNode & subtreeRoot);
    void addLevel(HashTreeNode & leafNode);
    void findAndVisitLeaves(HashTreeNode & subtreeRoot,
                            std::vector<unsigned int>::const_iterator start,
                            std::vector<unsigned> const& transactionItems,
                            int transactionID);
    static void visitLeaf(HashTreeNode & leaf, std::vector<unsigned> const& transactionItems, int tID);
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
