#include "Apriori.h"

#include <algorithm>
#include <cassert>

#include "easylogging++.h"

namespace algos {

void Apriori::GenerateCandidates(std::vector<Node>& children) {
    auto const last_child_iter = std::prev(children.end());
    for (auto child_iter = children.begin(); child_iter != last_child_iter; ++child_iter) {
        for (auto child_right_sibling_iter = std::next(child_iter);
                  child_right_sibling_iter != children.end(); ++child_right_sibling_iter) {
            std::vector<unsigned> items = child_iter->items;
            items.push_back(child_right_sibling_iter->items.back());

            if (!CanBePruned(items)) {
                candidates_[&(*child_iter)].emplace_back(std::move(items));
            }
        }
    }
}

void Apriori::CreateFirstLevelCandidates() {
    for (unsigned item_id = 0; item_id < transactional_data_->GetUniverseSize(); ++item_id) {
        candidates_[&root_].emplace_back(item_id);
    }
    ++level_num_;
}

bool Apriori::GenerateNextCandidateLevel() {
    std::stack<Node*> path;
    path.push(&root_);

    assert(level_num_ >= 2);
    while (!path.empty()) {
        auto node = path.top();
        path.pop();
        if (node->items.size() == level_num_ - 2 && !node->children.empty()) {
            GenerateCandidates(node->children);
        } else {
            UpdatePath(path, node->children);
        }
    }

    ++level_num_;
    return candidate_hash_tree_->Size() > 0;
}

void Apriori::UpdatePath(std::stack<Node*>& path, std::vector<Node>& vertices) {
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter) {
        Node* node_ptr = &(*iter);
        path.push(node_ptr);
    }
}

void Apriori::UpdatePath(std::stack<Node const*>& path, std::vector<Node> const& vertices) {
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter) {
        Node const* node_ptr = &(*iter);
        path.push(node_ptr);
    }
}

void Apriori::UpdatePath(std::queue<Node const*>& path, std::vector<Node> const& vertices) {
    for (auto const& vertex : vertices) {
        Node const* node_ptr = &vertex;
        path.push(node_ptr);
    }
}

bool Apriori::CanBePruned(std::vector<unsigned> const& itemset) {
    // we are able not to skip the last element, because without it the itemset is surely frequent
    assert(itemset.size() >= 2);
    for (unsigned index_to_skip = 0; index_to_skip < itemset.size() - 1; ++index_to_skip) {
        // TODO(alexandrsmirn) можно просто идти по листам вместо стека, пока не попадется пустой
        // std::list<Node> const& nodesToVisit = root.children;
        std::stack<Node*> nodes_to_visit;
        UpdatePath(nodes_to_visit, root_.children);

        unsigned item_index = 0;
        bool found_subset = false;

        while (!nodes_to_visit.empty()) {
            if (item_index == index_to_skip) {
                ++item_index;
            }

            unsigned next_item_id = itemset[item_index];  // item we want to find
            auto node = nodes_to_visit.top();
            nodes_to_visit.pop();

            if (node->items.back() == next_item_id) {
                // we found an item, so we go a level deeper
                // TODO(alexandrsmirn): вот тут кажется можно очистить стек, и заполнить новым, а не дополнять старое
                ++item_index;
                if (item_index == itemset.size()) {
                    // if we reached the final level
                    found_subset = true;
                    break;
                }
                nodes_to_visit = std::stack<Node*>();
                UpdatePath(nodes_to_visit, node->children);
            }
        }

        if (!found_subset) {
            // if at least one of subsets is not frequent, current itemset is not frequent too
            return true;
        }
    }

    return false;
}

unsigned long long Apriori::FindFrequent() {
    auto start_time = std::chrono::system_clock::now();

    CreateFirstLevelCandidates();
    while (!candidates_.empty()) {
        unsigned candidates_count = 0;
        for (auto const& [parent, candidate_children] : candidates_) {
            candidates_count += candidate_children.size();
        }
        auto const branching_degree = level_num_;
        auto const min_treshold = candidates_count / branching_degree + 1;

        candidate_hash_tree_ = std::make_unique<CandidateHashTree>(transactional_data_.get(),
                                                                   candidates_,
                                                                   branching_degree, min_treshold);
        candidate_hash_tree_->PerformCounting();
        candidate_hash_tree_->PruneNodes(minsup_);
        AppendToTree();
        candidates_.clear();
        GenerateNextCandidateLevel();
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    long long millis = elapsed_milliseconds.count();

    return millis;
}

unsigned long long Apriori::GenerateAllRules() {
    auto start_time = std::chrono::system_clock::now();

    std::queue<Node const*> path;
    UpdatePath(path, root_.children);
    unsigned long long frequent_count = 0;

    while (!path.empty()) {
        auto curr_node = path.front();
        path.pop();

        ++frequent_count;
        if (curr_node->items.size() >= 2) {
            GenerateRulesFrom(curr_node->items, curr_node->support);
        }
        UpdatePath(path, curr_node->children);
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    long long millis = elapsed_milliseconds.count();

    LOG(INFO) << "> Count of frequent itemsets: " << frequent_count;
    return millis;
}

std::list<std::set<std::string>> Apriori::GetFrequentList() const {
    std::list<std::set<std::string>> frequent_itemsets;

    std::queue<Node const*> path;
    UpdatePath(path, root_.children);

    while (!path.empty()) {
        auto const curr_node = path.front();
        path.pop();

        std::set<std::string> item_names;
        for (unsigned int item : curr_node->items) {
            item_names.insert(transactional_data_->GetItemUniverse()[item]);
        }

        frequent_itemsets.push_back(std::move(item_names));
        UpdatePath(path, curr_node->children);
    }

    return frequent_itemsets;
}

double Apriori::GetSupport(std::vector<unsigned int> const& frequent_itemset) const {
    auto const* path = &(root_.children);
    unsigned item_index = 0;
    auto node_comparator = [&item_index](Node const& node, std::vector<unsigned> const& items) {
        return node.items[item_index] < items[item_index];
    };

    while (item_index != frequent_itemset.size()) {
        auto const& node_vector = *path;
        auto next_node = std::lower_bound(node_vector.begin(), node_vector.end(), frequent_itemset,
                                          node_comparator);
        if (next_node == node_vector.end()) {
            break;
        } else if (item_index == frequent_itemset.size() - 1) {
            return next_node->support;
        } else {
            path = &(next_node->children);
        }
        ++item_index;
    }
    return -1;
}

void Apriori::AppendToTree() {
    for (auto& [node, children] : candidates_) {
        for (auto& child : children) {
            node->children.push_back(std::move(child));
        }
    }
}

} //namespace algos
