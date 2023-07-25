#include "algorithms/cfd/fd_first_algorithm.h"

#include <iterator>

#include <boost/unordered_map.hpp>
#include <easylogging++.h>

#include "algorithms/cfd/util/partition_tidlist_util.h"
#include "algorithms/cfd/util/partition_util.h"
#include "algorithms/cfd/util/set_util.h"
#include "algorithms/cfd/util/tidlist_util.h"
#include "util/config/equal_nulls/option.h"
#include "util/config/names_and_descriptions.h"
#include "util/config/option_using.h"

// see algorithms/cfd/LICENSE

namespace algos::cfd {

FDFirstAlgorithm::FDFirstAlgorithm() : CFDDiscovery({kDefaultPhaseName}) {
    RegisterOptions();
}

void FDFirstAlgorithm::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    RegisterInitialExecOption(Option{&min_supp_, kCfdMinimumSupport, kDCfdMinimumSupport, 0u});
    RegisterInitialExecOption(
            Option{&min_conf_, kCfdMinimumConfidence, kDCfdMinimumConfidence, 0.0});
    RegisterInitialExecOption(Option{&max_lhs_, kCfdMaximumLhs, kDCfdMaximumLhs, 0u});
    RegisterInitialExecOption(
            Option{&substrategy_, kCfdSubstrategy, kDCfdSubstrategy, +Substrategy::dfs});
}

void FDFirstAlgorithm::ResetStateCFD() {
    store_.clear();
    // cand_store_ does not have clear method
    all_attrs_.clear();
    free_map_.clear();
    free_itemsets_.clear();
    rules_.clear();
}

unsigned long long FDFirstAlgorithm::ExecuteInternal() {
    max_cfd_size_ = max_lhs_ + 1;
    CheckForIncorrectInput();
    auto start_time = std::chrono::system_clock::now();
    FdsFirstDFS();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    unsigned long long apriori_millis = elapsed_milliseconds.count();
    LOG(INFO) << "> CFD COUNT: " << cfd_list_.size();

    return apriori_millis;
}

void FDFirstAlgorithm::CheckForIncorrectInput() const {
    if (min_supp_ < 1) {
        throw std::invalid_argument("[ERROR] Illegal Support value: \"" +
                                    std::to_string(min_supp_) + "\"" + " is less than 1");
    }

    if (min_conf_ < 0 || min_conf_ > 1) {
        throw std::invalid_argument("[ERROR] Illegal Confidence value: \"" +
                                    std::to_string(min_conf_) + "\"" + " not in [0,1]");
    }

    if (max_cfd_size_ < 2) {
        throw std::invalid_argument("[ERROR] Illegal Max size value: \"" +
                                    std::to_string(max_cfd_size_) + "\"" + " is less than 1");
    }

    if (columns_number_ != 0 && tuples_number_ == 0) {
        throw std::invalid_argument(
                "[ERROR] Illegal columns_number and tuples_number values: columns_number is " +
                std::to_string(columns_number_) + " while tuples_number is 0");
    }

    if (tuples_number_ != 0 && columns_number_ == 0) {
        throw std::invalid_argument(
                "[ERROR] Illegal columns_number and tuples_number values: tuples_number is " +
                std::to_string(tuples_number_) + " while columnes_number is 0");
    }

    if (columns_number_ != 0 && tuples_number_ != 0 && min_supp_ > tuples_number_) {
        throw std::invalid_argument("[ERROR] Illegal Support value : " + std::to_string(min_supp_) +
                                    " is not in [1, " + std::to_string(tuples_number_) + "]");
    }
}

bool FDFirstAlgorithm::Precedes(const Itemset& a, const Itemset& b) {
    if (a.size() > b.size() || a == b) return false;
    Itemset pattern(relation_->GetAttrsNumber());
    for (int v : b) {
        if (v > 0) {
            pattern[relation_->GetAttrIndex(v)] = v;
        } else {
            pattern[-1 - v] = v;
        }
    }
    for (int i : a) {
        if (i > 0) {
            if (pattern[relation_->GetAttrIndex(i)] != i) return false;
        } else if (pattern[-1 - i] == 0) {
            return false;
        }
    }
    return true;
}

bool FDFirstAlgorithm::IsConstRule(const PartitionTIdList& items, int rhs_a) {
    int rhs_value;
    bool first = true;
    for (size_t pos_index = 0; pos_index <= items.tids.size(); pos_index++) {
        if (pos_index != items.tids.size() && items.tids[pos_index] != PartitionTIdList::SEP) {
            continue;
        }

        const auto tup = relation_->GetRow(items.tids[pos_index - 1]);
        if (first) {
            first = false;
            rhs_value = tup[rhs_a];
        } else if (tup[rhs_a] != rhs_value) {
            return false;
        }
    }
    return true;
}

void FDFirstAlgorithm::MineFD(const MinerNode<PartitionTIdList>& inode, const Itemset& lhs,
                              int rhs) {
    if (inode.tids.sets_number == 1 || IsConstRule(inode.tids, -1 - rhs)) return;
    const auto stored_sub = store_.find(lhs);
    if (stored_sub == store_.end()) {
        return;
    }
    bool lhs_gen = true;
    if (free_itemsets_.find(lhs) == free_itemsets_.end()) {
        lhs_gen = false;
    }
    if (rules_.find(rhs) != rules_.end()) {
        for (const auto& sub_rule : rules_[rhs]) {
            if (!std::any_of(sub_rule.begin(), sub_rule.end(),
                             [](int si) -> bool { return si < 0; }))
                continue;
            if (Precedes(sub_rule, lhs)) {
                lhs_gen = false;
            }
        }
    }
    if (lhs_gen) {
        // Here the confidence computing method from the paper is used
        double e = stored_sub->second.PartitionError(inode.tids);
        double conf = 1 - (e / TIdUtil::Support(stored_sub->second));
        if (conf >= min_conf_) {
            cfd_list_.emplace_back(lhs, rhs);
        }
        if (conf >= 1) {
            rules_[rhs].push_back(lhs);
        }
    }
}

// Initializing all objects that will be used in algorithm
void FDFirstAlgorithm::FdsFirstDFS() {
    all_attrs_ = Range(-static_cast<int>(relation_->GetAttrsNumber()), 0);
    PIdListMiners items = GetPartitionSingletons();
    for (auto& a : items) {
        a.candidates = all_attrs_;
        free_map_[std::make_pair(TIdUtil::Support(a.tids), a.tids.sets_number)].push_back(
                Itemset{a.item});
        free_itemsets_.insert(Itemset{a.item});
    }
    cand_store_ = PrefixTree<Itemset, Itemset>();
    store_[Itemset()] = PartitionTIdList(Iota(relation_->Size()));
    cand_store_.Insert(Itemset(), all_attrs_);
    FdsFirstDFS(Itemset(), items, substrategy_);
}

std::pair<std::vector<const PartitionTIdList*>, FDFirstAlgorithm::PIdListMiners>
FDFirstAlgorithm::ExpandMiningFd(const MinerNode<PartitionTIdList>& inode, int ix,
                                 const Itemset& iset, const PIdListMiners& items) const {
    const auto node_attrs = relation_->GetAttrVector(iset);
    std::vector<const PartitionTIdList*> expands;
    PIdListMiners tmp_suffix;
    for (int jx = static_cast<int>(items.size()) - 1; jx > ix; jx--) {
        const auto& jnode = items[jx];
        if (std::binary_search(node_attrs.begin(), node_attrs.end(),
                               relation_->GetAttrIndex(jnode.item)))
            continue;
        const Itemset newset = Join(iset, jnode.item);
        auto c = ConstructIntersection(inode.candidates, jnode.candidates);
        for (int zz : newset) {
            const auto zsub = ConstructSubset(newset, zz);
            const auto stored_sub = store_.find(zsub);
            if (stored_sub == store_.end()) {
                c.clear();
                break;
            }
            const Itemset& sub_cands = *cand_store_.Find(zsub);
            c = ConstructIntersection(c, sub_cands);
        }
        if (!c.empty()) {
            expands.emplace_back(&items[jx].tids);
            auto new_node = MinerNode<PartitionTIdList>(items[jx].item);
            new_node.candidates = c;
            new_node.prefix = ConstructSubset(newset, new_node.item);
            tmp_suffix.push_back(std::move(new_node));
        }
    }
    return {expands, tmp_suffix};
}

void FDFirstAlgorithm::FdsFirstDFS(const Itemset& prefix, const PIdListMiners& items,
                                   Substrategy ss) {
    for (int ix = static_cast<int>(items.size()) - 1; ix >= 0; ix--) {
        const MinerNode<PartitionTIdList>& inode = items[ix];
        const Itemset iset = Join(prefix, inode.item);
        const auto insect = ConstructIntersection(iset, inode.candidates);
        for (int out : insect) {
            const Itemset sub = ConstructSubset(iset, out);
            MineFD(inode, sub, out);

            if (ss == +Substrategy::dfs) {
                MinePatternsDFS(sub, out, inode.tids);
            } else if (ss == +Substrategy::bfs) {
                MinePatternsBFS(sub, out, inode.tids);
            }
        }

        if (inode.candidates.empty()) continue;
        if (iset.size() == max_cfd_size_) continue;
        store_[iset] = inode.tids;
        cand_store_.Insert(iset, inode.candidates);

        const auto [expands, tmp_suffix] = ExpandMiningFd(inode, ix, iset, items);

        const auto exps = PartitionTIdListUtil::ConstructIntersection(items[ix].tids, expands);
        PIdListMiners suffix;
        for (size_t e = 0; e < exps.size(); e++) {
            bool gen = true;
            const auto new_set = Join(tmp_suffix[e].prefix, tmp_suffix[e].item);
            const auto sp = std::make_pair(TIdUtil::Support(exps[e]), exps[e].sets_number);
            const auto free_map_elem = free_map_.find(sp);
            if (free_map_elem != free_map_.end()) {
                const auto free_cands = free_map_elem->second;
                for (const auto& sub_cand : free_cands) {
                    if (IsSubsetOf(sub_cand, new_set)) {
                        gen = false;
                        break;
                    }
                }
            }
            if (gen) {
                free_map_[sp].push_back(new_set);
                free_itemsets_.insert(new_set);
            }
            auto new_node = MinerNode<PartitionTIdList>(tmp_suffix[e].item, exps[e]);
            new_node.candidates = tmp_suffix[e].candidates;
            new_node.prefix = tmp_suffix[e].prefix;
            suffix.push_back(std::move(new_node));
        }
        if (!suffix.empty()) {
            std::sort(suffix.begin(), suffix.end(), [](const auto& a, const auto& b) {
                return a.tids.sets_number < b.tids.sets_number;
            });
            FdsFirstDFS(iset, suffix, ss);
        }
    }
}

void FDFirstAlgorithm::FillMinePatternsVars(PartitionList& partitions, RhsesPair2DList& pair_rhses,
                                            RuleIxs& rule_ixs, std::vector<int>& rhses,
                                            const Itemset& lhs, int rhs,
                                            const PartitionTIdList& all_tids) const {
    const auto lhs_attrs = relation_->GetAttrVector(lhs);
    unsigned count = 0;
    for (size_t pi = 0; pi <= all_tids.tids.size(); pi++) {
        if (pi != all_tids.tids.size() && all_tids.tids[pi] != PartitionTIdList::SEP) {
            count++;
            continue;
        }

        const Transaction& trans = relation_->GetRow(all_tids.tids[pi - 1]);
        const std::vector<int> lhs_constants = ConstructProjection(trans, lhs_attrs);
        const auto rule_i = rule_ixs.find(lhs_constants);
        if (rule_i == rule_ixs.end()) {
            rhses.push_back(trans[-1 - rhs]);
            auto rhses_pair_list = RhsesPairList();
            rhses_pair_list.emplace_back(trans[-1 - rhs], count);
            pair_rhses.push_back(std::move(rhses_pair_list));
            rule_ixs[lhs_constants] = partitions.size();
            partitions.emplace_back(lhs_constants, std::vector<unsigned>());
            partitions.back().second.push_back(count);
        } else {
            partitions[rule_i->second].second.push_back(count);
            pair_rhses[rule_i->second].emplace_back(trans[-1 - rhs], count);
        }
        count = 0;
    }
}

void FDFirstAlgorithm::AddCFDToCFDList(const std::vector<int>& sub, int out,
                                       const MinerNode<SimpleTIdList>& inode,
                                       const PartitionList& partitions) {
    bool lhs_gen = true;

    if (rules_.find(out) != rules_.end()) {
        for (const auto& sub_rule : rules_[out]) {
            if (out < 0 && !std::any_of(sub_rule.begin(), sub_rule.end(),
                                        [](int si) -> bool { return si < 0; }))
                continue;
            if (Precedes(sub_rule, sub)) {
                lhs_gen = false;
            }
        }
    }
    if (lhs_gen) {
        unsigned e = PartitionUtil::GetPartitionError(inode.tids, partitions);
        double conf = 1.0 - (static_cast<double>(e) / static_cast<double>(inode.node_supp));
        if (conf >= min_conf_) {
            cfd_list_.emplace_back(sub, out);
        }
        if (conf >= 1) {
            rules_[out].push_back(sub);
            if (out > 0) {
                rules_[-1 - relation_->GetAttrIndex(out)].push_back(sub);
            }
        }
    }
}

void FDFirstAlgorithm::AnalyzeCFDFromPIdList(const std::pair<int, SimpleTIdList>& item,
                                             const PartitionList& partitions,
                                             const std::vector<unsigned>& p_supps,
                                             TIdListMiners& items, const Itemset& lhs) {
    unsigned p_supp = PartitionUtil::GetPartitionSupport(item.second, p_supps);
    if (p_supp < min_supp_) {
        return;
    }
    const Itemset ns = Join(Itemset{item.first},
                            ConstructSubset(lhs, -1 - relation_->GetAttrIndex(item.first)));
    std::set<Itemset> nr_parts;
    for (int pid : item.second) {
        nr_parts.insert(partitions[pid].first);
    }
    bool gen = true;
    const auto sp = std::make_pair(p_supp, nr_parts.size());
    const auto free_map_pair = free_map_.find(sp);
    if (free_map_pair != free_map_.end()) {
        const auto free_cands = free_map_pair->second;
        for (const auto& sub_cand : free_cands) {
            if (IsSubsetOf(sub_cand, ns)) {
                gen = false;
                break;
            }
        }
    }
    if (gen) {
        free_map_[sp].push_back(ns);
        free_itemsets_.insert(ns);
    }
    items.emplace_back(item.first, item.second, p_supp);
}

bool FDFirstAlgorithm::FillFreeMapAndItemsets(const PartitionList& partitions, const Itemset& lhs,
                                              const Itemset& new_set, const SimpleTIdList& ij_tids,
                                              unsigned ij_supp) {
    if (ij_supp < min_supp_) {
        return false;
    }

    bool gen = true;
    const auto nas = relation_->GetAttrVectorItems(new_set);
    const Itemset ns = Join(new_set, SetDiff(lhs, nas));
    std::set<Itemset> nr_parts;
    for (int pid : ij_tids) {
        nr_parts.insert(partitions[pid].first);
    }
    const auto sp = std::make_pair(ij_supp, nr_parts.size());
    const auto free_map_pair = free_map_.find(sp);
    if (free_map_pair != free_map_.end()) {
        const auto free_cands = free_map_pair->second;
        for (const auto& sub_cand : free_cands) {
            if (IsSubsetOf(sub_cand, ns)) {
                gen = false;
                break;
            }
        }
    }
    if (gen) {
        free_map_[sp].push_back(ns);
        free_itemsets_.insert(ns);
    }
    return true;
}

void FDFirstAlgorithm::MinePatternsBFS(const Itemset& lhs, int rhs,
                                       const PartitionTIdList& all_tids) {
    std::map<int, SimpleTIdList> pid_lists;
    PartitionList partitions;
    RhsesPair2DList rhses_pairs;
    RuleIxs rule_ixs;
    std::vector<int> rhses;
    FillMinePatternsVars(partitions, rhses_pairs, rule_ixs, rhses, lhs, rhs, all_tids);
    TIdListMiners items;
    std::vector<unsigned> p_supps(partitions.size());
    int ri = 0;
    for (const auto& rule : partitions) {
        for (int i : rule.first) {
            pid_lists[i].push_back(ri);
        }
        unsigned s = std::accumulate(rule.second.begin(), rule.second.end(), 0u);
        p_supps[ri] = s;
        ri++;
    }

    for (const auto& item : pid_lists) {
        AnalyzeCFDFromPIdList(item, partitions, p_supps, items, lhs);
    }

    while (!items.empty()) {
        TIdListMiners suffix;
        for (size_t i = 0; i < items.size(); i++) {
            const auto& inode = items[i];
            Itemset iset = inode.prefix;
            iset.push_back(inode.item);
            const auto node_attrs = relation_->GetAttrVectorItems(iset);
            int out = (iset.size() == lhs.size()) ? GetMaxElem(rhses_pairs[inode.tids[0]]) : rhs;
            const auto sub = Join(iset, SetDiff(lhs, node_attrs));

            if (out > 0 || !PartitionUtil::IsConstRulePartition(inode.tids, rhses_pairs)) {
                AddCFDToCFDList(sub, out, inode, partitions);
            }
            for (size_t j = i + 1; j < items.size(); j++) {
                const auto& jnode = items[j];
                if (jnode.prefix != inode.prefix) continue;
                if (std::binary_search(node_attrs.begin(), node_attrs.end(),
                                       -1 - relation_->GetAttrIndex(jnode.item)))
                    continue;
                Itemset jset = jnode.prefix;
                jset.push_back(jnode.item);
                Itemset new_set = Join(jset, inode.item);
                SimpleTIdList ij_tids = ConstructIntersection(inode.tids, jnode.tids);
                unsigned ij_supp = PartitionUtil::GetPartitionSupport(ij_tids, p_supps);
                bool result = FillFreeMapAndItemsets(partitions, lhs, new_set, ij_tids, ij_supp);
                if (!result) continue;
                int jtem = new_set.back();
                new_set.pop_back();
                suffix.emplace_back(jtem, std::move(ij_tids), ij_supp, new_set);
            }
        }
        items.swap(suffix);
    }
}

void FDFirstAlgorithm::MinePatternsDFS(const Itemset& lhs, int rhs,
                                       const PartitionTIdList& all_tids) {
    std::map<int, SimpleTIdList> pid_lists;
    PartitionList partitions;
    RhsesPair2DList rhses_pairs;
    RuleIxs rule_ixs;
    std::vector<int> rhses;
    FillMinePatternsVars(partitions, rhses_pairs, rule_ixs, rhses, lhs, rhs, all_tids);
    TIdListMiners items;
    std::vector<unsigned> p_supps(partitions.size());
    int ri = 0;
    for (const auto& rule : partitions) {
        for (int i : rule.first) {
            pid_lists[i].push_back(ri);
        }
        unsigned s = std::accumulate(rule.second.begin(), rule.second.end(), 0u);
        p_supps[ri] = s;
        ri++;
    }

    for (const auto& item : pid_lists) {
        AnalyzeCFDFromPIdList(item, partitions, p_supps, items, lhs);
    }

    MinePatternsDFS(Itemset(), items, lhs, rhs, rhses_pairs, partitions, p_supps);
}

void FDFirstAlgorithm::MinePatternsDFS(const Itemset& prefix, TIdListMiners& items,
                                       const Itemset& lhs, int rhs, RhsesPair2DList& rhses_pair,
                                       PartitionList& partitions, std::vector<unsigned>& psupps) {
    for (int ix = static_cast<int>(items.size()) - 1; ix >= 0; ix--) {
        const auto& inode = items[ix];
        if (inode.tids.empty() && items[ix].tids.empty()) {
            LOG(INFO) << ix;
        }
        const Itemset iset = Join(prefix, inode.item);
        const auto node_attrs = relation_->GetAttrVectorItems(iset);
        int out = (iset.size() == lhs.size()) ? GetMaxElem(rhses_pair[inode.tids[0]]) : rhs;
        const auto sub = Join(iset, SetDiff(lhs, node_attrs));

        if (out > 0 || !PartitionUtil::IsConstRulePartition(inode.tids, rhses_pair)) {
            AddCFDToCFDList(sub, out, inode, partitions);
        }
        TIdListMiners suffix;
        for (size_t j = ix + 1; j < items.size(); j++) {
            const auto& jnode = items[j];
            if (std::binary_search(node_attrs.begin(), node_attrs.end(),
                                   -1 - relation_->GetAttrIndex(jnode.item)))
                continue;
            const Itemset new_set = Join(iset, jnode.item);
            SimpleTIdList ij_tids = ConstructIntersection(inode.tids, jnode.tids);
            unsigned ij_supp = PartitionUtil::GetPartitionSupport(ij_tids, psupps);

            bool result = FillFreeMapAndItemsets(partitions, lhs, new_set, ij_tids, ij_supp);
            if (!result) continue;
            suffix.emplace_back(jnode.item, ij_tids, ij_supp);
        }
        if (!suffix.empty()) {
            std::sort(suffix.begin(), suffix.end(), [](const auto& a, const auto& b) {
                return TIdUtil::Support(a.tids) < TIdUtil::Support(b.tids);
            });
            MinePatternsDFS(iset, suffix, lhs, rhs, rhses_pair, partitions, psupps);
        }
    }
}

// This method is used to get all partition singletons using MinerNode<PartitionTidList>
FDFirstAlgorithm::PIdListMiners FDFirstAlgorithm::GetPartitionSingletons() {
    std::vector<std::vector<std::vector<unsigned>>> partitions(relation_->GetAttrsNumber());
    std::unordered_map<int, std::pair<int, int>> attr_indices;

    for (size_t a = 0; a < relation_->GetAttrsNumber(); a++) {
        const auto& dom = relation_->GetDomain(a);
        partitions[a] = std::vector<std::vector<unsigned>>(dom.size());
        for (unsigned i = 0; i < dom.size(); i++) {
            partitions[a][i].reserve(relation_->Frequency(dom[i]));
            attr_indices[dom[i]] = std::make_pair(a, i);
        }
    }
    for (size_t row = 0; row < relation_->Size(); row++) {
        const auto& tup = relation_->GetRow(row);
        for (int item : tup) {
            const auto& attr_node_ix = attr_indices.at(item);
            partitions[attr_node_ix.first][attr_node_ix.second].push_back(row);
        }
    }

    PIdListMiners singletons;
    for (size_t a = 0; a < relation_->GetAttrsNumber(); a++) {
        int attr_item = -1 - static_cast<int>(a);
        auto new_node = MinerNode<PartitionTIdList>(attr_item);
        const auto& dom = relation_->GetDomain(a);
        new_node.tids.tids.reserve(relation_->Size() + dom.size() - 1);
        new_node.tids.sets_number = dom.size();
        for (unsigned i = 0; i < dom.size(); i++) {
            auto& ts = new_node.tids.tids;
            ts.insert(ts.end(), partitions[a][i].begin(), partitions[a][i].end());
            if (i != dom.size() - 1) {
                ts.push_back(PartitionTIdList::SEP);
            }
        }

        singletons.push_back(std::move(new_node));
    }
    return singletons;
}
}  // namespace algos::cfd
