#pragma once

#include "Itemset.h"

class AR {
public:
    std::vector<unsigned> left;   //antecedent
    std::vector<unsigned> right;  //consequent
    double confidence = -1;
    //TODO может быть не индексы, а сразу строки?

    AR() = default;
    AR(std::vector<unsigned> && left, std::vector<unsigned> && right, double confidence)
        :left(std::move(left)), right(std::move(right)), confidence(confidence) {}
};

class ARStrings {
public:
    std::list<std::string> left;   //antecedent
    std::list<std::string> right;  //consequent
    double confidence = -1;

    ARStrings() = default;
    ARStrings(std::list<std::string> && left, std::list<std::string> && right, double confidence)
            :left(std::move(left)), right(std::move(right)), confidence(confidence) {}

    ARStrings(AR const& idsRule, TransactionalData const* transactionalData)
            : confidence(idsRule.confidence) {
        auto const& itemNamesMap = transactionalData->getItemUniverse();

        for (auto itemID : idsRule.left) {
            this->left.push_back(itemNamesMap[itemID]);
        }
        for (auto itemID : idsRule.right) {
            this->right.push_back(itemNamesMap[itemID]);
        }
    }

    std::string toString() const {
        std::string result;
        result.append(std::to_string(confidence));
        result.append("\t{");
        for (auto const& itemName : left) {
            result.append(itemName);
            result.append(", ");
        }
        result.erase(result.size() - 2, 2);
        result.append("} -> {");
        for (auto const& itemName : right) {
            result.append(itemName);
            result.append(", ");
        }
        result.erase(result.size() - 2, 2);
        result.push_back('}');

        return result;
    }
};