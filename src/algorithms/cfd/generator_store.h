#pragma once

// see ./LICENSE

#include <map>

#include "algorithms/cfd/util/prefix_tree.h"
#include "model/cfd_types.h"
namespace algos {

template <typename T>
class GeneratorStore {
public:
    [[maybe_unused]] bool AddMinGen(const Itemset& new_set, int supp, const T& hash) {
        const auto& gen_it = generator_map_.find(hash);
        if (gen_it != generator_map_.end()) {
            if (gen_it->second.HasStrictSubset(new_set, supp)) {
                return false;
            }
            gen_it->second.Insert(new_set, supp);
        } else {
            generator_map_[hash].Insert(new_set, supp);
        }
        return true;
    }

    [[maybe_unused]] bool AddMinGen(const Itemset& newset, int supp, const T& hash,
                                    std::vector<Itemset>& subs) {
        const auto& gen_it = generator_map_.find(hash);
        if (gen_it != generator_map_.end()) {
            if (gen_it->second.HasSubset(newset, supp)) {
                const auto& ptSubs = gen_it->second.GetSets();
                subs.insert(subs.begin(), ptSubs.begin(), ptSubs.end());
                return false;
            }
            gen_it->second.Insert(newset, supp);
        } else {
            generator_map_[hash].Insert(newset, supp);
        }
        return true;
    }

    [[maybe_unused]] bool IsMinGen(const Itemset& newset, int supp, const T& hash) {
        const auto& gen_it = generator_map_.find(hash);
        if (gen_it != generator_map_.end()) {
            if (gen_it->second.HasStrictSubset(newset, supp)) {
                return false;
            }
        }
        return true;
    }

    [[maybe_unused]] std::vector<Itemset> GetMinGens(const Itemset& newset, int supp,
                                                     const T& hash) {
        const auto& gen_it = generator_map_.find(hash);
        if (gen_it != generator_map_.end()) {
            return gen_it->second.GetSubsets(newset, supp);
        }
    }

    [[maybe_unused]] std::map<T, PrefixTree<Itemset, int>> GetGenMap() {
        return generator_map_;
    }

private:
    std::map<T, PrefixTree<Itemset, int>> generator_map_;
};
} //namespace algos
