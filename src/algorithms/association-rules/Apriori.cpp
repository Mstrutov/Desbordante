#include "Apriori.h"

//#include "easylogging++.h"
#include <iostream>
#include <algorithm>

void Apriori::GenerateCandidates(std::vector<Node>& children) {
    auto const last_child_iter = children.end() - 1;
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

    while (!path.empty()) {
        auto node = path.top();
        path.pop();
        if (node->items.size() == level_num_ - 2 && !node->children.empty()) { //levelNumber is at least 2
            GenerateCandidates(node->children);
        } else {
            UpdatePath(path, node->children);
        }
    }

    ++level_num_;
    return candidate_hash_tree_->size() > 0;
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
    //последний элемент можем не убирать, так как мы добавили его к существующему
    for (unsigned index_to_skip = 0; index_to_skip < itemset.size() - 1; ++index_to_skip) { //itemset.size() is at least 2
        //std::list<Node> const& nodesToVisit = root.children; //TODO можно просто идти по листам вместо стека, пока не попадется пустой?????????
        std::stack<Node*> nodes_to_visit;
        UpdatePath(nodes_to_visit, root_.children);

        unsigned item_index = 0;
        bool found_subset = false;

        while (!nodes_to_visit.empty()) {
            if (item_index == index_to_skip) {
                ++item_index;
            }

            unsigned next_item_id = itemset[item_index];   //что хотим найти
            auto node = nodes_to_visit.top();
            nodes_to_visit.pop();

            if (node->items.back() == next_item_id) {
                //we found an item, so we go a level deeper
                //TODO вот тут кажется можно очистить стек, и заполнить новым, а не дополнять старое
                ++item_index;
                if (item_index == itemset.size()) {      //прошли необходимое количество вершин
                    found_subset = true;
                    break;
                }
                nodes_to_visit = std::stack<Node*>();
                UpdatePath(nodes_to_visit, node->children);
            }
        }

        if (!found_subset) {
            return true;                                   //если хотя бы одно подмножество не нашли, то бан
        }
    }

    return false;
}

unsigned long long Apriori::FindFrequent() {
    CreateFirstLevelCandidates();
    while (!candidates_.empty()) {
        unsigned candidates_count = 0;
        for (auto const& [parent, candidate_children] : candidates_){
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
    return 0;
}

unsigned long long Apriori::GenerateAllRules() {
    std::queue<Node const*> path;
    UpdatePath(path, root_.children);
    unsigned frequent_count = 0;

    while (!path.empty()) {
        auto curr_node = path.front();
        path.pop();

        ++frequent_count;
        std::cout << curr_node->support << " ";
        for (unsigned int item : curr_node->items) {
            std::cout << '<' << transactional_data_->GetItemUniverse()[item] << '>' << ' ';
        }
        std::cout << '\n';

        if (curr_node->items.size() >= 2) {
            GenerateRulesFrom(curr_node->items, curr_node->support);
        }
        UpdatePath(path, curr_node->children);
    }
    std::cout << frequent_count << '\n';
    return 0;
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
