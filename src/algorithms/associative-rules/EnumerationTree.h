#pragma once

#include <vector>
#include <list>
#include <stack>
#include <queue>

#include "CandidateHashTree.h"
#include "Itemset.h"
#include "ARAlgorithm.h"
#include "Node.h"

class EnumerationTree : public ARAlgorithm {
private:
    std::unique_ptr<CandidateHashTree> candidateHashTree; //TODO где-то инициализировать

    Node root;
    //std::vector<unsigned> currentItemset;
    unsigned levelNumber = 1;

    bool generateNextCandidateLevel(); //or list?
    void performCounting();

    bool canBePruned(std::vector<unsigned> const& itemset);
    static void updatePath(std::stack<Node*> & path, std::list<Node> & vertices);
    void generateCandidates(Node* node);
    static void updatePath(std::stack<Node const*> & path, std::list<Node> const& vertices);
    static void updatePath(std::queue<Node const*> & path, std::list<Node> const& vertices);

    virtual unsigned long long generateAllRules() override;
    virtual unsigned long long findFrequent() override;
public:
    explicit EnumerationTree(double minsup, double minconf,
                             std::filesystem::path const& path,
                             char separator = ',',
                             bool hasHeader = true)
            : ARAlgorithm(minsup, minconf, path, separator, hasHeader) {}
};
