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
    std::unordered_map<Node*, std::list<Node>>& candidates;

    TransactionalData const* const transactionalData = nullptr;

    struct LeafRow {
        NodeIterator candidateNode;
        Node* const parent;
        unsigned transactionCount = 0;

        //LeafRow() = default;
        LeafRow(LeafRow&& other) = default;
        LeafRow(NodeIterator node, Node* parent)
            : candidateNode(node), parent(parent) {}

        LeafRow(LeafRow const& other) = delete;
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
    void addCandidates();
public:
    CandidateHashTree(TransactionalData const* transactionalData,
                      std::unordered_map<Node*, std::list<Node>>& candidates,
                      unsigned branchingDegree, unsigned minThreshold)
            : branchingDegree(branchingDegree), minThreshold(minThreshold),
              candidates(candidates), transactionalData(transactionalData), root(1) {
        addCandidates();
    }

    void addCandidate(NodeIterator candidate, Node* parent); //maybe reference is possible
    unsigned size() const noexcept { return totalRowCount; };
    void performCounting();
    void pruneNodes(double minsup);
};
