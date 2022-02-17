#include "EnumerationTree.h"

//#include "easylogging++.h"
#include <iostream>

void EnumerationTree::generateCandidates(Node* node) {
    auto& children = node->children;

    auto const lastChildIter = std::next(children.end(), -1);
    for (auto childIter = children.begin(); childIter != lastChildIter; ++childIter) {
        for (auto childRightSiblingIter = std::next(childIter); childRightSiblingIter != children.end(); ++childRightSiblingIter) {
            std::vector<unsigned> items = childIter->items;
            items.push_back(childRightSiblingIter->items.back());

            if (!canBePruned(items)) {
                //Node candidate(std::move(items));
                childIter->children.emplace_back(std::move(items));                 //нужна ли перегрузка от rvalue? или оставить по значению
                candidateHashTree->addCandidate(std::next(childIter->children.end(), -1), &(childIter->children));        //добавляем итератор на только что добавленный в дерево кандидатов
            }
        }
    }
}

bool EnumerationTree::generateNextCandidateLevel() {
    if (levelNumber == 1) {
        for (unsigned itemID = 0; itemID < transactionalData->getUniverseSize(); ++itemID) {
            root.children.emplace_back(itemID);
            candidateHashTree->addCandidate(std::next(root.children.end(), -1), &root.children);
        }
    } else {
        std::stack<Node*> path; //TODO может итераторы?
        path.push(&root);

        while (!path.empty()) {
            auto node = path.top();
            path.pop();
            if (node->items.size() == levelNumber - 2) { //levelNumber is at least 2
                generateCandidates(node);
            } else {
                updatePath(path, node->children);
            }
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
    candidateHashTree = std::make_unique<CandidateHashTree>(transactionalData.get(), 10, 100);
    while (generateNextCandidateLevel()) {
        candidateHashTree->performCounting();
        candidateHashTree->pruneNodes(minsup);
        candidateHashTree = std::make_unique<CandidateHashTree>(transactionalData.get(), 10, 100);
    }
    return 0;
}

unsigned long long EnumerationTree::generateAllRules() {
    //TODO обходим дерево frequentов и для каждого вызываем generateAR, полученное складываем в arCollection
    //где применяем minconf?
    std::queue<Node const*> path;
    unsigned rowCount = 0;
    updatePath(path, root.children);

    while (!path.empty()) {
        auto currNode = path.front();
        path.pop();

        //generateRules(currNode->items);
        for (unsigned int item : currNode->items) {
            //LOG(DEBUG) << item;
            std::cout << '<' << transactionalData->getItemUniverse()[item] << '>' << ' ';
        }
        //LOG(DEBUG) << '\n';
        std::cout << '\n';
        ++rowCount;

        updatePath(path, currNode->children);
    }
    std::cout << rowCount << '\n';
    return 0;
}
