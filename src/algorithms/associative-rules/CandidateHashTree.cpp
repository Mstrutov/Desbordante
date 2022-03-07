#include "CandidateHashTree.h"

#include <algorithm>
#include "cassert"

void CandidateHashTree::AppendRow(LeafRow row, HashTreeNode& subtree_root) {
    if (!subtree_root.children.empty()) {
        auto const hash = HashFunction(row, subtree_root.level_number);
        AppendRow(std::move(row), subtree_root.children[hash]);
    } else {
        subtree_root.candidates.push_back(std::move(row));

        /* If the number of candidates in a leaf node is more than min_thresold, a leaf node becomes
         * an internal node and the tree expands. But if there is no more levels to expand (maximum
         * level number equals to the cardinality of a candidates), min_threshold is ignored
         * and a new candidates are just appended without trying to further grow the tree.*/
        if (subtree_root.candidates.size() > min_threshold_ &&
            subtree_root.level_number <= row.candidate_node->items.size()
            ) {
            AddLevel(subtree_root);
        }
    }
}

void CandidateHashTree::AddLevel(HashTreeNode& leaf_node) {
    unsigned const next_level_number = leaf_node.level_number + 1;

    //by this we make leaf node interior
    leaf_node.children.reserve(branching_degree_);
    for (unsigned siblingNum = 0; siblingNum < branching_degree_; ++siblingNum) {
        //generate new leaf nodes
        leaf_node.children.emplace_back(next_level_number);
    }
    leaf_node.children.shrink_to_fit();

    //distribute rows of an old leaf between new leaves
    for (auto& row : leaf_node.candidates) {
        AppendRow(std::move(row), leaf_node);
    }

    leaf_node.candidates.clear();
}

void CandidateHashTree::AddCandidate(NodeIterator candidate, Node* parent) {
    AppendRow(LeafRow(candidate, parent), root_);
    ++total_row_count_;
}

unsigned CandidateHashTree::HashFunction(LeafRow const& node_row, unsigned level_num) const {
    auto const& node_items = node_row.candidate_node->items;
    assert(level_num <= node_items.size());
    auto const curr_level_element_id = node_items.at(level_num - 1);
    return ItemHash(curr_level_element_id);
}

void CandidateHashTree::FindAndVisitLeaves(HashTreeNode& subtree_root,
                                           std::vector<unsigned>::const_iterator start,
                                           std::vector<unsigned> const& transaction_items,
                                           int tid) {
    unsigned const next_branch_number = ItemHash(*start);
    auto& next_node = subtree_root.children[next_branch_number];
    if (next_node.children.empty()) {
        //if nextNode is a leaf itself, then we visit it and terminate the recursion
        VisitLeaf(next_node, transaction_items, tid);
    } else {
        for (auto new_start = std::next(start); new_start != transaction_items.end(); ++new_start) {
            FindAndVisitLeaves(next_node, new_start, transaction_items, tid);
        }
    }
}

void CandidateHashTree::VisitLeaf(HashTreeNode& leaf,
                                  std::vector<unsigned> const& transaction_items, int tID) {
    if (leaf.last_visited_transaction_id == tID) { return; }

    leaf.last_visited_transaction_id = tID;

    for (auto& row : leaf.candidates) {
        auto const& candidate_items = row.candidate_node->items;
        if (std::includes(transaction_items.begin(), transaction_items.end(),
                          candidate_items.begin(), candidate_items.end())) {
            row.transaction_count++;
        }
    }
}

void CandidateHashTree::PerformCounting() {
    for (auto const& transaction : transactional_data_->GetTransactions()) {
        auto const& items = transaction.second.GetItemsIDs();
        auto const tid = transaction.first;
        if (root_.children.empty()) {
            //if the root is a leaf itself
            VisitLeaf(root_, items, tid);
        } else {
            for (auto start = items.begin(); start != items.end(); ++start) {
                FindAndVisitLeaves(root_, start, items, tid);
            }
        }
    }
}

void CandidateHashTree::Prune(double minsup, HashTreeNode& subtree_root) {
    if (subtree_root.children.empty()) {
        for (auto& row : subtree_root.candidates) {
            double const support = static_cast<double>(row.transaction_count) /
                                   transactional_data_->GetNumTransactions();
            if (support < minsup) {
                candidates_[row.parent].erase(row.candidate_node);
            } else {
                row.candidate_node->support = support;
            }
        }
    } else {
        for (auto& sibling : subtree_root.children) {
            Prune(minsup, sibling);
        }
    }
}

void CandidateHashTree::PruneNodes(double minsup) {
    Prune(minsup, root_);
}

void CandidateHashTree::AddCandidates() {
    for (auto& [node, candidate_children] : candidates_) {
        for (auto child = candidate_children.begin(); child != candidate_children.end(); ++child) {
            AddCandidate(child, node);
        }
    }
}
