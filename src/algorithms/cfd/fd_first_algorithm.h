#pragma once

#include "algorithms/cfd/cfd_discovery.h"
#include "algorithms/cfd/enums.h"
#include "algorithms/cfd/miner_node.h"
#include "algorithms/cfd/structures/partition_tidlist.h"
#include "algorithms/cfd/util/prefix_tree.h"

// see algorithms/cfd/LICENSE

namespace algos::cfd {

class FDFirstAlgorithm : public algos::cfd::CFDDiscovery {
    using PIdListMiners = std::vector<MinerNode<PartitionTIdList>>;
    using TIdListMiners = std::vector<MinerNode<SimpleTIdList>>;

private:
    unsigned min_supp_;
    unsigned max_cfd_size_;
    unsigned max_lhs_;
    double min_conf_;
    Substrategy substrategy_ = Substrategy::dfs;

    std::map<Itemset, PartitionTIdList> store_;
    PrefixTree<Itemset, Itemset> cand_store_;
    Itemset all_attrs_;
    std::map<std::pair<int, int>, std::vector<Itemset>> free_map_;
    std::set<Itemset> free_itemsets_;
    std::unordered_map<int, std::vector<Itemset>> rules_;

    void ResetStateCFD() final;
    void CheckForIncorrectInput() const;

    void FdsFirstDFS();
    void FdsFirstDFS(const Itemset &, const std::vector<MinerNode<PartitionTIdList>> &,
                     Substrategy = Substrategy::dfs);
    void MinePatternsBFS(const Itemset &lhs, int rhs, const PartitionTIdList &all_tids);
    void MinePatternsDFS(const Itemset &lhs, int rhs, const PartitionTIdList &all_tids);
    void MinePatternsDFS(const Itemset &, std::vector<MinerNode<SimpleTIdList>> &, const Itemset &,
                         int, RhsesPair2DList &, PartitionList &, std::vector<unsigned> &);
    std::vector<MinerNode<PartitionTIdList>> GetPartitionSingletons();

    bool Precedes(const Itemset &a, const Itemset &b);
    bool IsConstRule(const PartitionTIdList &items, int rhs_a);

    void MineFD(const MinerNode<PartitionTIdList> &inode, const Itemset &lhs, int rhs);
    std::pair<std::vector<const PartitionTIdList *>, std::vector<MinerNode<PartitionTIdList>>>
    ExpandMiningFd(const MinerNode<PartitionTIdList> &inode, int ix, const Itemset &iset,
                   const std::vector<MinerNode<PartitionTIdList>> &items) const;

    void FillMinePatternsVars(PartitionList &, RhsesPair2DList &, RuleIxs &, std::vector<int> &,
                              const Itemset &, int, const PartitionTIdList &) const;

    void AddCFDToCFDList(const std::vector<int> &sub, int out,
                         const MinerNode<SimpleTIdList> &inode, const PartitionList &partitions);

    void AnalyzeCFDFromPIdList(const std::pair<int, SimpleTIdList> &, const PartitionList &,
                               const std::vector<unsigned> &,
                               std::vector<MinerNode<SimpleTIdList>> &, const Itemset &);

    bool FillFreeMapAndItemsets(const PartitionList &partitions, const Itemset &lhs,
                                const Itemset &new_set, const SimpleTIdList &ij_tids, unsigned);

protected:
    void RegisterOptions();
    unsigned long long ExecuteInternal() final;

public:
    explicit FDFirstAlgorithm();
};
}  // namespace algos::cfd
