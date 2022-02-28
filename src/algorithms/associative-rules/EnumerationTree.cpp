#include "EnumerationTree.h"

//#include "easylogging++.h"
#include <iostream>
#include <algorithm>

void EnumerationTree::generateCandidates(std::list<Node>& children) {
    auto const lastChildIter = std::next(children.end(), -1);
    for (auto childIter = children.begin(); childIter != lastChildIter; ++childIter) {
        for (auto childRightSiblingIter = std::next(childIter); childRightSiblingIter != children.end(); ++childRightSiblingIter) {
            std::vector<unsigned> items = childIter->items;
            items.push_back(childRightSiblingIter->items.back());

            if (!canBePruned(items)) {
                candidates[&(*childIter)].emplace_back(std::move(items));
            }
        }
    }
}

void EnumerationTree::createFirstLevelCandidates() {
    for (unsigned itemID = 0; itemID < transactionalData->getUniverseSize(); ++itemID) {
        candidates[&root].emplace_back(itemID);
    }
    ++levelNumber;
}

bool EnumerationTree::generateNextCandidateLevel() {
    std::stack<Node*> path;
    path.push(&root);

    while (!path.empty()) {
        auto node = path.top();
        path.pop();
        if (node->items.size() == levelNumber - 2) { //levelNumber is at least 2
            generateCandidates(node->children);
        } else {
            updatePath(path, node->children);
        }
    }

    ++levelNumber;
    return candidateHashTree->size() > 0;
}

void EnumerationTree::updatePath(std::stack<Node*> & path, std::list<Node> & vertices) {
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter) {
        Node* nodePtr = &(*iter);
        path.push(nodePtr);
    }
}

void EnumerationTree::updatePath(std::stack<Node const*> & path, std::list<Node> const& vertices) {
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter) {
        Node const* nodePtr = &(*iter);
        path.push(nodePtr);
    }
}

void EnumerationTree::updatePath(std::queue<Node const*> & path, std::list<Node> const& vertices) {
    for (auto const& vertex : vertices) {
        Node const* nodePtr = &vertex;
        path.push(nodePtr);
    }
}

bool EnumerationTree::canBePruned(std::vector<unsigned> const& itemset) {
    //последний элемент можем не убирать, так как мы добавили его к существующему
    for (unsigned indexToSkip = 0; indexToSkip < itemset.size() - 1; ++indexToSkip) { //itemset.size() is at least 2
        //std::list<Node> const& nodesToVisit = root.children; //TODO можно просто идти по листам вместо стека, пока не попадется пустой?????????
        std::stack<Node*> nodesToVisit;
        updatePath(nodesToVisit, root.children);

        unsigned itemIndex = 0;
        bool foundSubset = false;

        while (!nodesToVisit.empty()) {
            if (itemIndex == indexToSkip) {
                ++itemIndex;
            }

            unsigned nextItemID = itemset[itemIndex];   //что хотим найти
            auto node = nodesToVisit.top();
            nodesToVisit.pop();

            if (node->items.back() == nextItemID) {
                //we found an item, so we go a level deeper
                //TODO вот тут кажется можно очистить стек, и заполнить новым, а не дополнять старое
                ++itemIndex;
                if (itemIndex == itemset.size()) {      //прошли необходимое количество вершин
                    foundSubset = true;
                    break;
                }
                nodesToVisit = std::stack<Node*>();
                updatePath(nodesToVisit, node->children);
            }
        }

        if (!foundSubset) {
            return true;                                   //если хотя бы одно подмножество не нашли, то бан
        }
    }

    return false;
}

unsigned long long EnumerationTree::findFrequent() {
    //TODO branching degree и minThreshold сделать зависимыми от номера уровня кандидатов
    createFirstLevelCandidates();
    while (!candidates.empty()) {
        candidateHashTree = std::make_unique<CandidateHashTree>(transactionalData.get(), candidates, 5, 5);
        candidateHashTree->performCounting();
        candidateHashTree->pruneNodes(minsup);
        appendToTree();
        candidates.clear();
        generateNextCandidateLevel();
    }
    return 0;
}

unsigned long long EnumerationTree::generateAllRules() {
    std::queue<Node const*> path;
    updatePath(path, root.children);

    while (!path.empty()) {
        auto currNode = path.front();
        path.pop();

        if (currNode->items.size() >= 2) {
            generateRulesFrom(currNode->items, currNode->support);
        }
        updatePath(path, currNode->children);
    }

    return 0;
}

std::list<std::set<std::string>> EnumerationTree::getAllFrequent() const {
    std::list<std::set<std::string>> frequentItemsets;

    std::queue<Node const*> path;
    updatePath(path, root.children);

    while (!path.empty()) {
        auto const currNode = path.front();
        path.pop();

        std::set<std::string> itemNames;
        for (unsigned int item : currNode->items) {
            itemNames.insert(transactionalData->getItemUniverse()[item]);
        }

        frequentItemsets.push_back(std::move(itemNames));
        updatePath(path, currNode->children);
    }

    return frequentItemsets;
}

double EnumerationTree::getSupport(std::vector<unsigned int> const& frequentItemset) const {
    std::list<Node> const* path = &(root.children);
    unsigned itemIndex = 0;
    while (itemIndex != frequentItemset.size()) {
        for (auto const& node : *path) {
            if (node.items[itemIndex] == frequentItemset[itemIndex]) {
                if (itemIndex == frequentItemset.size() - 1) {
                    return node.support;
                } else {
                    path = &(node.children);
                    break;
                }
            }
        }
        ++itemIndex;
    }
    return -1;
}

void EnumerationTree::appendToTree() {
    for (auto& [node, children] : candidates) {
        for (auto& child : children) {
            node->children.push_back(std::move(child));
        }
    }
}
