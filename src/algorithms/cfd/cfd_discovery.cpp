#include "algorithms/cfd/cfd_discovery.h"

// see ./LICENSE

#include <iterator>
#include <thread>

#include <boost/unordered_map.hpp>
#include <easylogging++.h>

#include "algorithms/cfd/partition_table.h"
#include "algorithms/cfd/util/set_util.h"
#include "algorithms/cfd/util/cfd_output_util.h"
#include "algorithms/options/equal_nulls_opt.h"

namespace algos {

decltype(CFDDiscovery::MinSupportOpt) CFDDiscovery::MinSupportOpt{
        {config::names::kCfdMinimumSupport, config::descriptions::kDCfdMinimumSupport}, 0};

decltype(CFDDiscovery::ColumnsNumberOpt) CFDDiscovery::ColumnsNumberOpt{
        {config::names::kCfdColumnsNumber, config::descriptions::kDCfdColumnsNumber}, 0};

decltype(CFDDiscovery::TuplesNumberOpt) CFDDiscovery::TuplesNumberOpt{
        {config::names::kCfdTuplesNumber, config::descriptions::kDCfdTuplesNumber}, 0};

decltype(CFDDiscovery::MinConfidenceOpt) CFDDiscovery::MinConfidenceOpt{
        {config::names::kCfdMinimumConfidence, config::descriptions::kDCfdMinimumConfidence}, 0.0};

decltype(CFDDiscovery::MaxLhsSizeOpt) CFDDiscovery::MaxLhsSizeOpt{
        {config::names::kCfdMaximumLhs, config::descriptions::kDCfdMaximumLhs}, 0};

decltype(CFDDiscovery::AlgoOpt) CFDDiscovery::AlgoOpt{
        {config::names::kCfdAlgo, config::descriptions::kDCfdAlgo}};

CFDDiscovery::CFDDiscovery(std::vector<std::string_view> phase_names)
    : Primitive(std::move(phase_names)) {
    RegisterOptions();
    MakeOptionsAvailable(config::GetOptionNames(config::EqualNullsOpt, TuplesNumberOpt, ColumnsNumberOpt));
}

CFDDiscovery::CFDDiscovery() : CFDDiscovery({kDefaultPhaseName}) {}

void CFDDiscovery::FitInternal(model::IDatasetStream& data_stream) {
    if (relation_ == nullptr) {
        relation_ =
                CFDRelationData::CreateFrom(data_stream, is_null_equal_null_,
                                            columns_number_, tuples_number_);
    }

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: CFD mining is meaningless.");
    }
}

void CFDDiscovery::RegisterOptions() {
    RegisterOption(config::EqualNullsOpt.GetOption(&is_null_equal_null_));
    RegisterOption(MinSupportOpt.GetOption(&min_supp_));
    RegisterOption(MinConfidenceOpt.GetOption(&min_conf_));
    RegisterOption(TuplesNumberOpt.GetOption(&tuples_number_));
    RegisterOption(ColumnsNumberOpt.GetOption(&columns_number_));
    RegisterOption(MaxLhsSizeOpt.GetOption(&max_lhs_));
    RegisterOption(AlgoOpt.GetOption(&algo_));
}

void CFDDiscovery::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable(config::GetOptionNames(MinSupportOpt, MinConfidenceOpt, MaxLhsSizeOpt, AlgoOpt));
}

unsigned long long CFDDiscovery::ExecuteInternal() {
    max_cfd_size_ = max_lhs_ + 1;
    CheckForIncorrectInput();
    auto start_time = std::chrono::system_clock::now();

    if (algo_ == +CfdAlgo::fd_first_dfs_dfs) {
        FdsFirstDFS(min_supp_, max_cfd_size_, kDfs, min_conf_);
    }
    else if (algo_ == +CfdAlgo::fd_first_dfs_bfs) {
        FdsFirstDFS(min_supp_, max_cfd_size_, kBfs, min_conf_);
    }
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    long long apriori_millis = elapsed_milliseconds.count();
    LOG(INFO) << "> CFD COUNT: " << cfd_list_.size();

    return apriori_millis;
}

void CFDDiscovery::CheckForIncorrectInput() const {
    if (min_supp_ < 1) {
        throw std::invalid_argument(
                "[ERROR] Illegal Support value: \"" + std::to_string(min_supp_) +
                "\"" + " is less than 1");
    }

    if (min_conf_ < 0 || min_conf_ > 1) {
        throw std::invalid_argument(
                "[ERROR] Illegal Confidence value: \"" + std::to_string(min_conf_) +
                "\"" + " not in [0,1]");
    }

    if (max_cfd_size_ < 2) {
        throw std::invalid_argument(
                "[ERROR] Illegal Max size value: \"" + std::to_string(max_cfd_size_) +
                "\"" + " is less than 1");
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
        throw std::invalid_argument(
                "[ERROR] Illegal Support value : " + std::to_string(min_supp_) +
                " is not in [1, " + std::to_string(tuples_number_) + "]");
    }
}

[[maybe_unused]] int CFDDiscovery::NrCfds() const {
    return (int)cfd_list_.size();
}

CFDList CFDDiscovery::GetCfds() const {
    return cfd_list_;
}

std::string CFDDiscovery::GetCfdString(CFD const& cfd) const {
    return Output::CFDToString(cfd, relation_);
}

std::string CFDDiscovery::GetRelationString(char delim) const {
    return relation_->GetStringFormat(delim);
}

std::string CFDDiscovery::GetRelationString(const SimpleTidList& subset, char delim) const {
    return relation_->GetStringFormat(subset, delim);
}

bool CFDDiscovery::Precedes(const Itemset& a, const Itemset& b) {
    if (a.size() > b.size()) return false;
    if (a == b) return false;
    Itemset pattern(relation_->GetAttrsNumber());
    for (int v : b) {
        if (v > 0) pattern[relation_->GetAttrIndex(v)] = v;
        else pattern[-1-v] = v;
    }
    for (int i : a) {
        if (i > 0) {
            int bv = pattern[relation_->GetAttrIndex(i)];
            if (bv != i) {
                return false;
            }
        }
        else {
            if (pattern[-1-i] == 0) return false;
        }
    }
    return true;
}

bool CFDDiscovery::IsConstRule(const PartitionTidList& items, int rhs_a) {
    int rhs_value;
    bool first = true;
    for (size_t pos_index = 0; pos_index <= items.tids.size(); pos_index++) {
        if (pos_index == items.tids.size() || items.tids[pos_index] == PartitionTidList::SEP) {
            auto tup = relation_->GetRow(items.tids[pos_index - 1]);
            if (first) {
                first = false;
                rhs_value = tup[rhs_a];
            }
            else {
                if (tup[rhs_a] != rhs_value) return false;
            }
        }
    }
    return true;
}

// Initializing all objects that will be used in algorithm
void CFDDiscovery::FdsFirstDFS(int min_supp, unsigned max_lhs_size, SUBSTRATEGY ss, double min_conf) {
    min_supp_ = min_supp;
    min_conf_ = min_conf;
    max_cfd_size_ = max_lhs_size;
    all_attrs_ = Range(-(int)relation_->GetAttrsNumber(), 0);
    PartitionTable::database_row_number_ = (int)relation_->Size();
    std::vector<MinerNode<PartitionTidList> > items = GetPartitionSingletons();
    for (auto& a : items) {
        a.cands = all_attrs_;
        free_map_[std::make_pair(support(a.tids),a.tids.sets_number)].push_back(Itemset{a.item});
        free_itemsets_.insert(Itemset{a.item});
    }
    cand_store_ = PrefixTree<Itemset, Itemset>();
    store_[Itemset()] = Convert(Iota((int)relation_->Size()));
    cand_store_.Insert(Itemset(), all_attrs_);
    FdsFirstDFS(Itemset(), items, ss);
}

void CFDDiscovery::FdsFirstDFS(const Itemset &prefix, std::vector<MinerNode<PartitionTidList> > &items, SUBSTRATEGY ss) {
    for (int ix = items.size() - 1; ix >= 0; ix--) {
        MinerNode<PartitionTidList>& inode = items[ix];
        const Itemset iset = Join(prefix, inode.item);

        auto insect = Intersection(iset, inode.cands);
        for (int out : insect) {
            Itemset sub = Subset(iset, out);
            if (inode.tids.sets_number == 1 || IsConstRule(inode.tids, -1 - out)) continue;
            auto stored_sub = store_.find(sub);
            if (stored_sub == store_.end()) {
                continue;
            }
            bool lhs_gen = true;
            if (free_itemsets_.find(sub) == free_itemsets_.end()) {
                lhs_gen = false;
            }
            if (rules_.find(out) != rules_.end()) {
                for (const auto& sub_rule : rules_[out]) {
                    if (!Has(sub_rule, [](int si) -> bool { return si < 0; })) continue;
                    if (Precedes(sub_rule, sub)) {
                        lhs_gen = false;
                    }
                }
            }
            if (lhs_gen) {
                // Here the confidence computing method from the paper is used
                double e = PartitionTable::PartitionError(stored_sub->second, inode.tids);
                double conf = 1 - (e / support(stored_sub->second));
                if (conf >= min_conf_) {
                    cfd_list_.emplace_back(sub, out);
                }
                if (conf >= 1) {
                    rules_[out].push_back(sub);
                }
            }
            if (ss == SUBSTRATEGY::kDfs)
                MinePatternsDFS(sub, out, inode.tids);
            else if (ss == SUBSTRATEGY::kBfs)
                 MinePatternsBFS(sub, out, inode.tids);
        }
        if (inode.cands.empty()) continue;
        if (iset.size() == max_cfd_size_) continue;
        store_[iset] = inode.tids;
        cand_store_.Insert(iset, inode.cands);
        auto node_attrs = relation_->GetAttrVector(iset);
        std::vector<const PartitionTidList*> expands;
        std::vector<MinerNode<PartitionTidList> > tmp_suffix;
        for (int jx = (int)items.size() - 1; jx > ix; jx--) {
            const auto &jnode = items[jx];
            if (std::binary_search(node_attrs.begin(), node_attrs.end(),
                                   relation_->GetAttrIndex(jnode.item))) continue;
            Itemset newset = Join(iset, jnode.item);
            auto c = Intersection(inode.cands, jnode.cands);
            for (int zz : newset) {
                auto zsub = Subset(newset, zz);
                auto stored_sub = store_.find(zsub);
                if (stored_sub == store_.end()) {
                    c.clear();
                    break;
                }
                const Itemset & sub_cands = *cand_store_.Find(zsub);
                c = Intersection(c, sub_cands);
            }
            if (c.size()) {
                expands.push_back(&items[jx].tids);
                tmp_suffix.emplace_back(items[jx].item);
                tmp_suffix.back().cands = c;
                tmp_suffix.back().prefix = Subset(newset, tmp_suffix.back().item);
            }
        }
        const auto exps = PartitionTable::Intersection(items[ix].tids, expands);
        std::vector<MinerNode<PartitionTidList> > suffix;

        for (size_t e = 0; e < exps.size(); e++) {
            bool gen = true;
            auto new_set = Join(tmp_suffix[e].prefix, tmp_suffix[e].item);
            auto sp = std::make_pair(support(exps[e]), exps[e].sets_number);
            auto free_map_elem = free_map_.find(sp);
            if (free_map_elem != free_map_.end()) {
                auto free_cands = free_map_elem->second;
                for (const auto &sub_cand : free_cands) {
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
            suffix.emplace_back(tmp_suffix[e].item, exps[e]);
            suffix.back().cands = tmp_suffix[e].cands;
            suffix.back().prefix = tmp_suffix[e].prefix;
        }
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [](const MinerNode<PartitionTidList>& a, const MinerNode<PartitionTidList>& b) {
                return a.tids.sets_number < b.tids.sets_number;
            });
            FdsFirstDFS(iset, suffix, ss);
        }
    }
}



int CFDDiscovery::GetPartitionSupport(const SimpleTidList& pids, const std::vector<int>& psupps) {
    int res = 0;
    for (int p : pids) {
        res += psupps[p];
    }
    return res;
}

int CFDDiscovery::GetPartitionError(const SimpleTidList& pids, const std::vector<std::pair<Itemset, std::vector<int> > >& partitions) {
    int res = 0;
    for (int p : pids) {
        int max = 0;
        int total = 0;
        for (const auto& rhs : partitions[p].second) {
            total += rhs;
            if (rhs > max) {
                max = rhs;
            }
        }
        res += total - max;
    }
    return res;
}

bool CFDDiscovery::IsConstRulePartition(const SimpleTidList &items,
                                        const std::vector<std::vector<std::pair<int, int> > > &rhses) {
    int rhs_value;
    bool first = true;
    for (int pi : items) {
        if (first) {
            first = false;
            rhs_value = rhses[pi][0].first;
        }
        for (auto rp : rhses[pi]) {
            int r = rp.first;
            if (r != rhs_value) return false;
        }
    }
    return true;
}

void CFDDiscovery::MinePatternsBFS(const Itemset &lhs, int rhs, const PartitionTidList & all_tids) {
    std::map<int, SimpleTidList> pid_lists;
    auto lhs_attrs = relation_->GetAttrVector(lhs);
    std::vector<std::pair<Itemset, std::vector<int> > > partitions;
    std::vector<std::vector<std::pair<int,int> > > rhses_pairs;
    boost::unordered_map<Itemset, int> rule_ixs;
    std::vector<int> rhses;
    int count = 0;
    for (size_t pi = 0; pi <= all_tids.tids.size(); pi++) {
        if (pi == all_tids.tids.size() || all_tids.tids[pi] == PartitionTidList::SEP ) {
            const Transaction& trans = relation_->GetRow(all_tids.tids[pi - 1]);
            Itemset lhs_constants = Projection(trans, lhs_attrs);
            auto rule_i = rule_ixs.find(lhs_constants);
            if (rule_i == rule_ixs.end()) {
                rhses.push_back(trans[-1-rhs]);
                rhses_pairs.emplace_back();
                rhses_pairs.back().emplace_back(trans[-1 - rhs],count);
                rule_ixs[lhs_constants] = partitions.size();
                partitions.emplace_back(lhs_constants, std::vector<int>());
                partitions.back().second.push_back(count);
            }
            else {
                partitions[rule_i->second].second.push_back(count);
                rhses_pairs[rule_i->second].emplace_back(trans[-1 - rhs], count);
            }
            count = 0;
        }
        else {
            count++;
        }
    }
    std::vector<MinerNode<SimpleTidList> > items;
    std::vector<int> p_supps(partitions.size());
    int ri = 0;
    for (const auto& rule : partitions) {
        for (int i : rule.first) {
            pid_lists[i].push_back(ri);
        }
        int s = 0;
        for (const auto& rhs_int : rule.second) {
            s += rhs_int;
        }
        p_supps[ri] = s;
        ri++;
    }
    for (const auto& item : pid_lists) {
        int p_supp = GetPartitionSupport(item.second, p_supps);
        if (p_supp >= (int)min_supp_) {
            Itemset ns = Join(Itemset{item.first},
                              Subset(lhs, -1 - relation_->GetAttrIndex(item.first)));
            std::set<Itemset> nr_parts;
            for (int pid : item.second) {
                nr_parts.insert(partitions[pid].first);
            }
            bool gen = true;
            auto sp = std::make_pair(p_supp, nr_parts.size());
            auto free_map_pair = free_map_.find(sp);
            if (free_map_pair != free_map_.end()) {
                auto free_cands = free_map_pair->second;
                for (const auto & sub_cand : free_cands) {
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
    }
    while (!items.empty()) {
        std::vector<MinerNode<SimpleTidList> > suffix;
        for (size_t i = 0; i < items.size(); i++) {
            const auto &inode = items[i];
            Itemset iset = inode.prefix;
            iset.push_back(inode.item);
            auto node_attrs = relation_->GetAttrVectorItems(iset);
            bool lhs_gen = true;
            int out = (iset.size() == lhs.size()) ? GetMaxElem(rhses_pairs[inode.tids[0]]) : rhs;
            auto sub = Join(iset, SetDiff(lhs, node_attrs));

            if (out > 0 || !IsConstRulePartition(inode.tids, rhses_pairs)) {
                if (rules_.find(out) != rules_.end()) {
                    for (const auto & sub_rule : rules_[out]) {
                        if (out < 0 && !Has(sub_rule, [](int si) -> bool { return si < 0; })) continue;
                        if (Precedes(sub_rule, sub)) {
                            lhs_gen = false;
                        }
                    }
                }
                if (lhs_gen) {
                    int e = GetPartitionError(inode.tids, partitions);
                    double conf = 1.0 - ((double) e / (double) inode.node_supp);
                    if (conf >= min_conf_) {
                        cfd_list_.emplace_back(sub, out);
                    }
                    if (conf >= 1) {
                        rules_[out].push_back(sub);
                        if (out > 0) rules_[-1 - relation_->GetAttrIndex(out)].push_back(sub);
                    }
                }
            }
            for (size_t j = i + 1; j < items.size(); j++) {
                const auto& jnode = items[j];
                if (jnode.prefix != inode.prefix) continue;
                if (std::binary_search(node_attrs.begin(), node_attrs.end(), -1 - relation_->GetAttrIndex(jnode.item))) continue;
                Itemset jset = jnode.prefix;
                jset.push_back(jnode.item);
                Itemset new_set = Join(jset, inode.item);
                SimpleTidList ij_tids = Intersection(inode.tids, jnode.tids);
                int ij_supp = GetPartitionSupport(ij_tids, p_supps);
                if (ij_supp >= (int)min_supp_){
                    bool gen = true;
                    auto nas = relation_->GetAttrVectorItems(new_set);
                    Itemset ns = Join(new_set, SetDiff(lhs, nas));
                    std::set<Itemset> nr_parts;
                    for (int pid : ij_tids) {
                        nr_parts.insert(partitions[pid].first);
                    }
                    auto sp = std::make_pair(ij_supp, nr_parts.size());
                    auto free_map_pair = free_map_.find(sp);
                    if (free_map_pair != free_map_.end()) {
                        auto free_cands = free_map_pair->second;
                        for (const auto & sub_cand : free_cands) {
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
                    int jtem = new_set.back();
                    new_set.pop_back();
                    suffix.emplace_back(jtem, std::move(ij_tids), ij_supp, new_set);
                }
            }
        }
        items.swap(suffix);
    }
}

void CFDDiscovery::MinePatternsDFS(const Itemset &lhs, int rhs, const PartitionTidList & all_tids) {
    std::map<int, SimpleTidList> pid_lists;
    auto lhs_attrs = relation_->GetAttrVector(lhs);
    std::vector<std::pair<Itemset, std::vector<int> > > partitions;
    std::vector<std::vector<std::pair<int,int> > > pair_rhses;
    boost::unordered_map<Itemset, int> rule_ixs;
    std::vector<int> rhses;
    int count = 0;
    for (size_t pi = 0; pi <= all_tids.tids.size(); pi++) {
        if (pi == all_tids.tids.size() || all_tids.tids[pi] == PartitionTidList::SEP ) {
            const Transaction& trans = relation_->GetRow(all_tids.tids[pi - 1]);
            std::vector<int> lhs_constants = Projection(trans, lhs_attrs);
            auto rule_i = rule_ixs.find(lhs_constants);
            if (rule_i == rule_ixs.end()) {
                rhses.push_back(trans[-1-rhs]);
                pair_rhses.emplace_back();
                pair_rhses.back().emplace_back(trans[-1-rhs],count);
                rule_ixs[lhs_constants] = (int)partitions.size();
                partitions.emplace_back(lhs_constants, std::vector<int>());
                partitions.back().second.push_back(count);
            }
            else {
                partitions[rule_i->second].second.push_back(count);
                pair_rhses[rule_i->second].emplace_back(trans[-1-rhs], count);
            }
            count = 0;
        }
        else {
            count++;
        }
    }
    std::vector<MinerNode<SimpleTidList> > items;
    std::vector<int> psupps(partitions.size());
    int ri = 0;
    for (const auto& rule : partitions) {
        for (int i : rule.first) {
            pid_lists[i].push_back(ri);
        }
        int s = 0;
        for (const auto& rhs_items : rule.second) {
            s += rhs_items;
        }
        psupps[ri] = s;
        ri++;
    }
    for (const auto& item : pid_lists) {
        int psupp = GetPartitionSupport(item.second, psupps);
        if (psupp >= (int)min_supp_) {
            Itemset ns = Join(Itemset{item.first},
                              Subset(lhs, -1 - relation_->GetAttrIndex(item.first)));
            std::set<Itemset> nr_parts;
            for (int pid : item.second) {
                nr_parts.insert(partitions[pid].first);
            }
            bool gen = true;
            auto sp = std::make_pair(psupp, nr_parts.size());
            auto free_map_pair = free_map_.find(sp);
            if (free_map_pair != free_map_.end()) {
                auto free_cands = free_map_pair->second;
                for (const auto &sub_cand : free_cands) {
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
            items.emplace_back(item.first, item.second, psupp);
        }
    }
    MinePatternsDFS(Itemset(), items, lhs, rhs, pair_rhses, partitions, psupps);
}

void CFDDiscovery::MinePatternsDFS(const Itemset &prefix, std::vector<MinerNode<SimpleTidList> > &items,
                                   const Itemset &lhs, int rhs,
                                   std::vector<std::vector<std::pair<int, int> > > & rhses_pair,
                                   std::vector<std::pair<Itemset, std::vector<int> > > &partitions,
                                   std::vector<int> &psupps) {
    for (int ix = (int)items.size() - 1; ix >= 0; ix--) {
        const auto &inode = items[ix];
        Itemset iset = Join(prefix, inode.item);
        auto node_attrs = relation_->GetAttrVectorItems(iset);
        bool lhs_gen = true;
        int out = (iset.size() == lhs.size()) ? GetMaxElem(rhses_pair[inode.tids[0]]) : rhs;
        auto sub = Join(iset, SetDiff(lhs, node_attrs));

        if (out > 0 || !IsConstRulePartition(inode.tids, rhses_pair)) {
            if (rules_.find(out) != rules_.end()) {
                for (const auto & sub_rule : rules_[out]) {
                    if (out < 0 && !Has(sub_rule, [](int si) -> bool { return si < 0; })) continue;
                    if (Precedes(sub_rule, sub)) {
                        lhs_gen = false;
                    }
                }
            }
            if (lhs_gen) {
                int e = GetPartitionError(inode.tids, partitions);
                double conf = 1.0 - ((double) e / (double) inode.node_supp);
                if (conf >= min_conf_) {
                    cfd_list_.emplace_back(sub, out);
                }
                if (conf >= 1) {
                    rules_[out].push_back(sub);
                    if (out > 0) rules_[-1 - relation_->GetAttrIndex(out)].push_back(sub);
                }
            }
        }
        std::vector<MinerNode<SimpleTidList> > suffix;
        for (size_t j = ix + 1; j < items.size(); j++) {
            const auto &jnode = items[j];
            if (std::binary_search(node_attrs.begin(), node_attrs.end(), -1- relation_->GetAttrIndex(jnode.item))) continue;
            Itemset newset = Join(iset, jnode.item);
            SimpleTidList ij_tids = Intersection(inode.tids, jnode.tids);
            int ij_supp = GetPartitionSupport(ij_tids, psupps);
            if (ij_supp >= (int)min_supp_){
                bool gen = true;
                auto nas = relation_->GetAttrVectorItems(newset);
                Itemset ns = Join(newset, SetDiff(lhs, nas));
                std::set<Itemset> nr_parts;
                for (int pid : ij_tids) {
                    nr_parts.insert(partitions[pid].first);
                }
                auto sp = std::make_pair(ij_supp, nr_parts.size());
                auto free_map_pair = free_map_.find(sp);
                if (free_map_pair != free_map_.end()) {
                    auto free_cands = free_map_pair->second;
                    for (const auto &sub_cand : free_cands) {
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
                suffix.emplace_back(jnode.item, std::move(ij_tids), ij_supp);
            }
        }
        if (suffix.size()) {
            std::sort(suffix.begin(), suffix.end(), [](const MinerNode<SimpleTidList>& a, const MinerNode<SimpleTidList>& b) {
                return support(a.tids) < support(b.tids);
            });
            MinePatternsDFS(iset, suffix, lhs, rhs, rhses_pair, partitions, psupps);
        }
    }
}

// This method is used to get all partition singletons using MinerNode<PartitionTidList>
std::vector<MinerNode<PartitionTidList> > CFDDiscovery::GetPartitionSingletons() {
    std::vector<std::vector<SimpleTidList> > partitions(relation_->GetAttrsNumber());
    std::unordered_map<int, std::pair<int,int> > attr_indices;
    for (size_t a = 0; a < relation_->GetAttrsNumber(); a++) {
        const auto& dom = relation_->GetDomain(a);
        partitions[a] = std::vector<SimpleTidList>(dom.size());
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
    std::vector<MinerNode<PartitionTidList> > singletons;
    for (size_t a = 0; a < relation_->GetAttrsNumber(); a++) {
        int attr_item = -1 - a;
        singletons.emplace_back(attr_item);
        const auto& dom = relation_->GetDomain(a);
        singletons.back().tids.tids.reserve(relation_->Size()+dom.size()-1);
        singletons.back().tids.sets_number = dom.size();
        for (unsigned i = 0; i < dom.size(); i++) {
            auto& ts = singletons.back().tids.tids;
            ts.insert(ts.end(), partitions[a][i].begin(), partitions[a][i].end());
            if (i < dom.size() - 1) {
                ts.push_back(PartitionTidList::SEP);
            }
        }
        singletons.back().HashTids();
    }
    return singletons;
}
} //namespace algos
