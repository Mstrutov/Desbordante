#include <algorithm>
#include <stack>

#include "ARAlgorithm.h"

unsigned long long ARAlgorithm::execute() {
    transactionalData = TransactionalData::createFrom(inputGenerator, inputFormat, hasTransactionID);
    if (transactionalData->getNumTransactions() == 0) {
        throw std::runtime_error("Got an empty .csv file: AR mining is meaningless.");
    }

    auto time = findFrequent();
    time += generateAllRules();

    return time;
}

void ARAlgorithm::registerARStrings(AR const& rule) {
    arCollection.emplace_back(rule, transactionalData.get());
}

void ARAlgorithm::updatePath(std::stack<RuleNode*> & path, std::list<RuleNode> & vertices) {
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter) {
        RuleNode* nodePtr = &(*iter);
        path.push(nodePtr);
    }
}

void ARAlgorithm::generateRulesFrom(std::vector<unsigned int> const& frequentItemset, double support) {
    root.children.clear();
    for (auto itemID : frequentItemset) {
        std::vector<unsigned> rhs {itemID};
        std::vector<unsigned> lhs;
        std::set_difference(frequentItemset.begin(), frequentItemset.end(),
                            rhs.begin(), rhs.end(),
                            std::back_inserter(lhs));
        auto const lhsSupport = getSupport(lhs);
        auto const confidence = support / lhsSupport;
        if (confidence >= minconf) {
            root.children.emplace_back(std::move(lhs), std::move(rhs), confidence);
            arCollection.emplace_back(root.children.back().rule, transactionalData.get());
        }
    }
    if (root.children.empty()) {
        return;
    }

    unsigned levelNumber = 2;
    while (generateRuleLevel(frequentItemset, support, levelNumber)) {
        ++levelNumber;
    }
}

bool ARAlgorithm::generateRuleLevel(std::vector<unsigned> const& frequentItemset, double support, unsigned levelNumber) {
    bool generatedAny = false;
    std::stack<RuleNode*> path;
    path.push(&root);

    while (!path.empty()) {
        auto node = path.top(); //TODO попробовать как-то синхронизировать пусть с генерацией
        path.pop();
        if (node->rule.right.size() == levelNumber - 2) { //levelNumber is at least 2
            generatedAny = mergeRules(frequentItemset, support, node);
        } else {
            updatePath(path, node->children);
        }
    }

    return generatedAny;
}

bool ARAlgorithm::mergeRules(std::vector<unsigned> const& frequentItemset, double support, RuleNode* node) {
    auto& children = node->children;
    bool produceRule = false;

    auto const lastChildIter = std::next(children.end(), -1);
    for (auto childIter = children.begin(); childIter != lastChildIter; ++childIter) {
        for (auto childRightSiblingIter = std::next(childIter); childRightSiblingIter != children.end(); ++childRightSiblingIter) {
            std::vector<unsigned> rhs = childIter->rule.right;
            rhs.push_back(childRightSiblingIter->rule.right.back());

            std::vector<unsigned> lhs;
            std::set_difference(frequentItemset.begin(), frequentItemset.end(),
                                rhs.begin(), rhs.end(),
                                std::back_inserter(lhs));

            auto const lhsSupport = getSupport(lhs);
            auto const confidence = support / lhsSupport;
            if (confidence >= minconf) {
                childIter->children.emplace_back(std::move(lhs), std::move(rhs), confidence);                 //нужна ли перегрузка от rvalue? или оставить по значению
                arCollection.emplace_back(childIter->children.back().rule, transactionalData.get());
                produceRule = true;
            }
        }
    }
    return produceRule;
}
