#include "ARAlgorithm.h"

#include <algorithm>
#include <cassert>
#include <iostream>

namespace algos {

unsigned long long ARAlgorithm::Execute() {
    transactional_data_ = model::TransactionalData::CreateFrom(input_generator_, *input_format_);
    if (transactional_data_->GetNumTransactions() == 0) {
        throw std::runtime_error("Got an empty .csv file: AR mining is meaningless.");
    }

    auto time = FindFrequent();
    time += GenerateAllRules();

    for (auto const& rule : ar_collection_) {
        std::cout << model::ARStrings(rule, transactional_data_.get()).ToString() << '\n';
    }

    return time;
}

void ARAlgorithm::UpdatePath(std::stack<RuleNode*>& path, std::list<RuleNode>& vertices) {
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter) {
        RuleNode* node_ptr = &(*iter);
        path.push(node_ptr);
    }
}

void ARAlgorithm::GenerateRulesFrom(std::vector<unsigned> const& frequent_itemset, double support) {
    root_.children.clear();
    for (auto item_id : frequent_itemset) {
        std::vector<unsigned> rhs{item_id};
        std::vector<unsigned> lhs;
        std::set_difference(frequent_itemset.begin(), frequent_itemset.end(),
                            rhs.begin(), rhs.end(),
                            std::back_inserter(lhs));
        auto const lhs_support = GetSupport(lhs);
        auto const confidence = support / lhs_support;
        if (confidence >= minconf_) {
            auto const& new_ar = ar_collection_.emplace_back(std::move(lhs), std::move(rhs),
                                                             confidence);
            root_.children.emplace_back(new_ar);
        }
    }
    if (root_.children.empty()) {
        return;
    }

    unsigned level_number = 2;
    while (GenerateRuleLevel(frequent_itemset, support, level_number)) {
        ++level_number;
    }
}

bool ARAlgorithm::GenerateRuleLevel(std::vector<unsigned> const& frequent_itemset, double support,
                                    unsigned level_number) {
    bool generated_any = false;
    std::stack<RuleNode*> path;
    path.push(&root_);

    assert(level_number >= 2);
    while (!path.empty()) {
        // TODO(alexandrsmirn) попробовать как-то синхронизировать путь с генерацией
        auto node = path.top();
        path.pop();
        if (node->rule.right.size() == level_number - 2) {  // levelNumber is at least 2
            generated_any = MergeRules(frequent_itemset, support, node);
        } else {
            UpdatePath(path, node->children);
        }
    }

    return generated_any;
}

bool ARAlgorithm::MergeRules(std::vector<unsigned> const& frequent_itemset, double support,
                             RuleNode* node) {
    auto& children = node->children;
    bool is_rule_produced = false;

    auto const last_child_iter = std::prev(children.end());
    for (auto child_iter = children.begin(); child_iter != last_child_iter; ++child_iter) {
        for (auto child_right_sibling_iter = std::next(child_iter);
                  child_right_sibling_iter != children.end(); ++child_right_sibling_iter) {
            std::vector<unsigned> rhs = child_iter->rule.right;
            rhs.push_back(child_right_sibling_iter->rule.right.back());
            if (rhs.size() == frequent_itemset.size()) {
                // in this case the LHS is empty
                continue;
            }

            std::vector<unsigned> lhs;
            std::set_difference(frequent_itemset.begin(), frequent_itemset.end(),
                                rhs.begin(), rhs.end(),
                                std::back_inserter(lhs));

            auto const lhs_support = GetSupport(lhs);
            auto const confidence = support / lhs_support;
            if (confidence >= minconf_) {
                auto const& new_ar = ar_collection_.emplace_back(std::move(lhs), std::move(rhs),
                                                                 confidence);
                child_iter->children.emplace_back(new_ar);
                is_rule_produced = true;
            }
        }
    }
    return is_rule_produced;
}

std::list<model::ARStrings> ARAlgorithm::GetArStringsList() const {
    std::list<model::ARStrings> ar_strings;
    for (auto const& ar : ar_collection_) {
        ar_strings.emplace_back(ar, transactional_data_.get());
    }
    return ar_strings;
}

} // namespace algos
