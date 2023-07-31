#include "order.h"

#include <algorithm>

#include "config/names_and_descriptions.h"
#include "config/tabular_data/input_table/option.h"
#include "model/types/types.h"

namespace algos::order {

Order::Order() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName()});
}

void Order::RegisterOptions() {
    using namespace config::names;
    using namespace config::descriptions;
    using config::Option;

    RegisterOption(config::kTableOpt(&input_table_));
}

void Order::LoadDataInternal() {
    typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, false);
}

void Order::ResetState() {}

void Order::CreateSortedPartitions() {
    std::vector<model::TypedColumnData> const& data = typed_relation_->GetColumnData();
    for (unsigned int i = 0; i < data.size(); ++i) {
        if (!model::Type::IsOrdered(data[i].GetTypeId())) {
            continue;
        }
        std::vector<std::pair<unsigned long, std::byte const*>> indexed_byte_data;
        indexed_byte_data.reserve(data[i].GetNumRows());
        std::vector<std::byte const*> const& byte_data = data[i].GetData();
        for (size_t k = 0; k < byte_data.size(); ++k) {
            indexed_byte_data.emplace_back(k, byte_data[k]);
        }
        model::Type const& type = data[i].GetType();
        auto less = [&type](std::pair<unsigned long, std::byte const*> l,
                            std::pair<unsigned long, std::byte const*> r) {
            return type.Compare(l.second, r.second) == model::CompareResult::kLess;
        };
        std::sort(indexed_byte_data.begin(), indexed_byte_data.end(), less);
        std::vector<std::unordered_set<unsigned long>> equivalence_classes;
        equivalence_classes.push_back({indexed_byte_data.front().first});
        auto equal = [&type](std::pair<unsigned long, std::byte const*> l,
                             std::pair<unsigned long, std::byte const*> r) {
            return type.Compare(l.second, r.second) == model::CompareResult::kEqual;
        };
        for (size_t k = 1; k < indexed_byte_data.size(); ++k) {
            if (equal(indexed_byte_data[k - 1], indexed_byte_data[k])) {
                equivalence_classes.back().insert(indexed_byte_data[k].first);
            } else {
                equivalence_classes.push_back({indexed_byte_data[k].first});
            }
        }
        sorted_partitions_[{i}] = SortedPartition(std::move(equivalence_classes));
    }
}

bool SubsetSetDifference(std::unordered_set<unsigned long>& a,
                                                        std::unordered_set<unsigned long>& b) {
    auto const not_found = b.end();
    for (auto const& element: a)
        if (b.find(element) == not_found) {
            return false;
        } else {
            b.erase(element);
        }
    return true;
}

ValidityType Order::CheckForSwap(SortedPartition const& l, SortedPartition const& r) {
    ValidityType res = ValidityType::valid;
    size_t l_i = 0, r_i = 0;
    bool next_l = true, next_r = true;
    std::unordered_set<unsigned long> l_eq_class;
    std::unordered_set<unsigned long> r_eq_class;
    while (l_i < l.sorted_partition.size() && r_i < r.sorted_partition.size()) {
        if (next_l) {
            l_eq_class = l.sorted_partition[l_i];
        }    
        if (next_r) {
            r_eq_class = r.sorted_partition[r_i];
        }
        if (l_eq_class.size() < r_eq_class.size()) {
            if(!SubsetSetDifference(l_eq_class, r_eq_class)) {
                return ValidityType::swap;
            } else {
                res = ValidityType::merge;
                ++l_i;
                next_l = true;
                next_r = false;
            }
        } else {
            if (!SubsetSetDifference(r_eq_class, l_eq_class)) {
                return ValidityType::swap;
            } else {
                ++r_i;
                next_r = true;
                if (l_eq_class.empty()) {
                    ++l_i;
                    next_l = true;
                } else {
                    next_l = false;
                }
            }
        }
    }
    return res;
}

std::vector<Order::AttributeList> GetPrefixes(Order::Node const& node) {
    std::vector<Order::AttributeList> res;
    res.reserve(node.size() - 1);
    for (size_t i = 1; i < node.size(); ++i) {
        res.emplace_back(node.begin(), node.begin() + i);
    }
    return res;
}

void Order::Prune(LatticeLevel& lattice_level) {
    for (auto node_it = lattice_level.begin(); node_it != lattice_level.end();) {
        bool all_candidates_empty = false;
        std::vector<AttributeList> prefixes = GetPrefixes(*node_it);
        for (AttributeList const& lhs : prefixes) {
            if (!candidate_sets_[lhs].empty()) {
                all_candidates_empty = false;
                break;
            } else {
                all_candidates_empty = true;
            }
        }
        if (all_candidates_empty) {
            node_it = lattice_level.erase(node_it);
        } else {
            ++node_it;
        }
    }
    for (auto candidate_it = candidate_sets_.begin(); candidate_it != candidate_sets_.end();) {
        if (candidate_it->second.empty()) {
            candidate_it = candidate_sets_.erase(candidate_it);
        } else {
            ++candidate_it;
        }
    }
}

std::vector<unsigned int> MaxPrefix(std::vector<unsigned int> const& attribute_list) {
    return std::vector<unsigned int>(attribute_list.begin(), attribute_list.end() - 1);
}

using PrefixBlocks = std::unordered_map<Order::AttributeList, std::vector<Order::Node>,
                                        boost::hash<std::vector<unsigned int>>>;

PrefixBlocks GetPrefixBlocks(Order::LatticeLevel const& l) {
    PrefixBlocks res;
    for (Order::Node const& node : l) {
        std::vector<unsigned int> node_prefix = MaxPrefix(node);
        if (res.find(node_prefix) == res.end()) {
            res[node_prefix] = {};
        }
        res[node_prefix].push_back(node);
    }
    return res;
}

Order::Node JoinNodes(Order::Node const& l, Order::Node const& r) {
    Order::Node res(l);
    res.push_back(r.back());
    return res;
}

Order::LatticeLevel Order::GenerateNextLevel(LatticeLevel const& l) {
    LatticeLevel next;
    PrefixBlocks prefix_blocks = GetPrefixBlocks(l);
    for (auto const& [prefix, prefix_block] : prefix_blocks) {
        for (Node const& node : prefix_block) {
            for (Node const& join_node : prefix_block) {
                if (node == join_node) {
                    continue;
                }
                Node joined = JoinNodes(node, join_node);
                sorted_partitions_[joined] =
                        sorted_partitions_[prefix] * sorted_partitions_[{join_node.back()}];
                next.insert(joined);
            }
        }
    }
    for (Node const& node : l) {
        candidate_sets_[node] = {};
    }
    return next;
}

unsigned long long Order::ExecuteInternal() {}

}  // namespace algos::order
