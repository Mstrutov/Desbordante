#include "list_lattice.h"

#include "order_utility.h"

namespace algos::order {

ListLattice::ListLattice(CandidateSets& candidate_sets,
                         std::vector<AttributeList> single_attributes) {
    lattice_level_ = {};
    for (AttributeList const& single_attribute : single_attributes) {
        lattice_level_.insert(single_attribute);
        candidate_sets[single_attribute] = {};
        for (AttributeList const& rhs_single_attribute : single_attributes) {
            if (single_attribute != rhs_single_attribute) {
                candidate_sets[single_attribute].insert(rhs_single_attribute);
            }
        }
    }
}

CandidatePairs ListLattice::ObtainCandidates(Node const& node) const {
    CandidatePairs res;
    res.reserve(node.size() - 1);
    for (size_t i = 1; i < node.size(); ++i) {
        AttributeList lhs(node.begin(), node.begin() + i);
        AttributeList rhs(node.begin() + i, node.end());
        res.emplace_back(lhs, rhs);
    }
    return res;
}

void ListLattice::Prune(CandidateSets& candidate_sets) {
    if (level_num_ < 2) {
        return;
    }
    for (auto node_it = lattice_level_.begin(); node_it != lattice_level_.end();) {
        bool all_candidates_empty = false;
        std::vector<AttributeList> prefixes = GetPrefixes(*node_it);
        for (AttributeList const& lhs : prefixes) {
            if (!candidate_sets[lhs].empty()) {
                all_candidates_empty = false;
                break;
            } else {
                all_candidates_empty = true;
            }
        }
        if (all_candidates_empty) {
            node_it = lattice_level_.erase(node_it);
        } else {
            ++node_it;
        }
    }
    /* TODO: Try to make iteration from metanome */
    for (auto candidate_it = candidate_sets.begin(); candidate_it != candidate_sets.end();) {
        if (candidate_it->second.empty()) {
            candidate_it = candidate_sets.erase(candidate_it);
        } else {
            ++candidate_it;
        }
    }
}

ListLattice::PrefixBlocks ListLattice::GetPrefixBlocks() const {
    PrefixBlocks res;
    for (Node const& node : lattice_level_) {
        std::vector<unsigned int> node_prefix = MaxPrefix(node);
        if (res.find(node_prefix) == res.end()) {
            res[node_prefix] = {};
        }
        res[node_prefix].push_back(node);
    }
    return res;
}

Node ListLattice::JoinNodes(Node const& l, Node const& r) const {
    Node res(l);
    res.push_back(r.back());
    return res;
}

void ListLattice::GenerateNextLevel(CandidateSets& candidate_sets) {
    LatticeLevel next;
    PrefixBlocks prefix_blocks = GetPrefixBlocks();
    for (auto const& [prefix, prefix_block] : prefix_blocks) {
        for (Node const& node : prefix_block) {
            for (Node const& join_node : prefix_block) {
                if (node == join_node) {
                    continue;
                }
                Node joined = JoinNodes(node, join_node);
                next.insert(joined);
            }
        }
    }
    if (level_num_ > 1 && !candidate_sets.empty()) {
        for (Node const& node : lattice_level_) {
            candidate_sets[node] = {};
        }
    }
    lattice_level_ = std::move(next);
    ++level_num_;
}

}  // namespace algos::order
