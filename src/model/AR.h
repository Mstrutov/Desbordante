#pragma once

#include "Itemset.h"
#include "TransactionalData.h"

namespace model {

struct ArIDs {
    std::vector<unsigned> left;   // antecedent
    std::vector<unsigned> right;  // consequent
    double confidence = -1;

    ArIDs() = default;
    ArIDs(std::vector<unsigned> left, std::vector<unsigned> right, double confidence)
        : left(std::move(left)), right(std::move(right)), confidence(confidence) {}

    ArIDs(ArIDs const& other) = default;
    ArIDs& operator=(ArIDs const& other) = default;
    ArIDs(ArIDs&& other) = default;
    ArIDs& operator=(ArIDs&& other) = default;
};

struct ARStrings {
    std::list<std::string> left;   // antecedent
    std::list<std::string> right;  // consequent
    double confidence = -1;

    ARStrings() = default;
    ARStrings(std::list<std::string> left, std::list<std::string> right, double confidence)
        : left(std::move(left)), right(std::move(right)), confidence(confidence) {}

    ARStrings(ArIDs const& id_format_rule, TransactionalData const* transactional_data)
        : confidence(id_format_rule.confidence) {
        std::vector<std::string> const& item_names_map = transactional_data->GetItemUniverse();

        for (auto itemID : id_format_rule.left) {
            this->left.push_back(item_names_map[itemID]);
        }
        for (auto itemID : id_format_rule.right) {
            this->right.push_back(item_names_map[itemID]);
        }
    }

    ARStrings(ARStrings const& other) = default;
    ARStrings& operator=(ARStrings const& other) = default;
    ARStrings(ARStrings&& other) = default;
    ARStrings& operator=(ARStrings&& other) = default;

    std::string ToString() const {
        std::string result;
        result.append(std::to_string(confidence));
        result.append("\t{");
        for (auto const& item_name : left) {
            result.append(item_name);
            result.append(", ");
        }
        result.erase(result.size() - 2, 2);
        result.append("} -> {");
        for (auto const& item_name : right) {
            result.append(item_name);
            result.append(", ");
        }
        result.erase(result.size() - 2, 2);
        result.push_back('}');

        return result;
    }
};

} // namespace model
