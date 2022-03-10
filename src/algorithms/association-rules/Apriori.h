#pragma once

#include <vector>
#include <list>
#include <stack>
#include <queue>

#include "CandidateHashTree.h"
#include "Itemset.h"
#include "ARAlgorithm.h"
#include "Node.h"

class Apriori : public ARAlgorithm {
private:
    std::unique_ptr<CandidateHashTree> candidate_hash_tree_; //TODO может убрать из полей, а создавать просто в методе?

    Node root_;
    std::unordered_map<Node*, std::list<Node>> candidates_;
    unsigned level_num_ = 1;

    bool GenerateNextCandidateLevel(); //or list?

    bool CanBePruned(std::vector<unsigned> const& itemset);
    static void UpdatePath(std::stack<Node*>& path, std::vector<Node>& vertices);
    void GenerateCandidates(std::vector<Node>& children);
    static void UpdatePath(std::stack<Node const*>& path, std::vector<Node> const& vertices);
    static void UpdatePath(std::queue<Node const*>& path, std::vector<Node> const& vertices);
    void CreateFirstLevelCandidates();
    void AppendToTree();

    double GetSupport(std::vector<unsigned> const& frequent_itemset) const override;
    unsigned long long GenerateAllRules() override;
    unsigned long long FindFrequent() override;
public:
    explicit Apriori(Config const& config)
        : ARAlgorithm(config, {"AR mining"}) {}

    std::list<std::set<std::string>> GetFrequentList() const override;
};
