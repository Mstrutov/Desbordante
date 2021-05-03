//
// Created by alexandrsmirn
//

#include "DFD.h"

#include <algorithm>
#include <random>

#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "PositionListIndex.h"
#include "PLICache.h"

unsigned long long DFD::execute() {
    //shared_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);//second parameter?
    shared_ptr<RelationalSchema> schema = relation->getSchema();

    //std::list<shared_ptr<Column>> possibleRHSs(schema->getColumns().begin(), schema->getColumns().end());

    /*std::list<shared_ptr<Column>> possibleRHSs;
    for (auto const& column : schema->getColumns()) {
        possibleRHSs.push_back(column);
    }

    //first loop of DFD
    for (auto columnIter = possibleRHSs.begin(); columnIter != possibleRHSs.end(); ) {
        shared_ptr<ColumnData> columnData = relation->getColumnData((*columnIter)->getIndex());
        shared_ptr<PositionListIndex> columnPLI = columnData->getPositionListIndex();
        auto nextIter = std::next(columnIter);

        //if current column is unique
        if (columnPLI->getNumNonSingletonCluster() == 0) {
            Vertical const& lhs = Vertical(**columnIter);//наверно стоит убрать, потому что при первом запуске же будет мув
            possibleRHSs.erase(columnIter);

            for (auto const& rhs : possibleRHSs) {
                this->registerFD(lhs, *rhs); //TODO ptrs???
            }
        }
        columnIter = nextIter;
    }*/
    //TODO может переделать на сет, потому что проверям лежит ли элемент там
    std::vector<Vertical> uniqueVerticals;

    for (auto const& column : schema->getColumns()) {
        shared_ptr<ColumnData> const columnData = relation->getColumnData(column->getIndex());
        shared_ptr<PositionListIndex> const columnPLI = columnData->getPositionListIndex();

        if (columnPLI->getNumNonSingletonCluster() == 0) {
            Vertical const& lhs = Vertical(*column);
            uniqueVerticals.push_back(lhs);
            /*for (auto const& rhs : schema->getColumns()) {
                if (rhs->getIndex() != column->getIndex()) {
                    registerFD(lhs, *rhs);
                }
            }*/
        }
    }

    //second loop of DFD
    for (auto const& rhs : schema->getColumns()) {
        //тут строим новую решетку, соответственно нужно обнулить/завести некоторые структуры данных
        minimalDeps.clear();
        maximalNonDeps.clear();
        dependenciesMap = DependenciesMap(relation->getSchema());
        nonDependenciesMap = NonDependenciesMap(relation->getSchema());
        observations.clear();

        //в метаноме немного по другому
        for (auto const& lhs : uniqueVerticals) {
            if (!lhs.contains(*rhs)) {
                observations[lhs] = NodeCategory::minimalDependency;
                //разобраться с шаред поинтерами
                dependenciesMap.addNewDependency(std::make_shared<Vertical>(lhs));//TODO костыль
                minimalDeps.insert(std::make_shared<Vertical>(lhs));
            }
        }

        findLHSs(rhs, schema);

        //TODO не заьыть зареать все зависимости из списка
        for (auto const& minimalDependencyLHS : minimalDeps) {
            registerFD(std::move(*minimalDependencyLHS), *rhs);
        }
    }
    std::cout << "====JSON-FD========\r\n" << FDAlgorithm::getJsonFDs() << std::endl;
    std::cout << "HASH: " << FDAlgorithm::fletcher16() << std::endl;
}

shared_ptr<Vertical> DFD::takeRandom(std::list<shared_ptr<Vertical>> & nodeList) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeList.begin(), nodeList.end()) - 1);
    auto iterator = nodeList.begin();
    std::advance(iterator, dis(this->gen));
    shared_ptr<Vertical> node = move(*iterator);
    nodeList.erase(iterator);
    return node;
}

shared_ptr<Vertical> DFD::takeRandom(std::vector<shared_ptr<Vertical>> const& nodeVector) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeVector.begin(), nodeVector.end()) - 1);
    auto iterator = nodeVector.begin();
    std::advance(iterator, dis(this->gen));
    return *iterator;
}

void DFD::findLHSs(shared_ptr<Column const> const& rhs, shared_ptr<RelationalSchema> schema) {

    //std::vector<shared_ptr<Vertical>> minimalDeps; //TODO мб их определять либо в функции execute, либо полями класса
    //std::vector<shared_ptr<Vertical>> maximalNonDeps;

    std::list<shared_ptr<Vertical>> seeds; //TODO лист или вектор?

    //initialize seeds nodes
    for (auto const& column : schema->getColumns()) {
        //add to seeds all columns except rhs
        if (column->getIndex() != rhs->getIndex())
            seeds.push_back(std::make_shared<Vertical>(*column)); //TODO проверить, поменял emplace на push
    }
    do {
        while (!seeds.empty()) {
            shared_ptr<Vertical> node = takeRandom(seeds);
            do {
                auto const nodeObservationIter = observations.find(*node); //const?

                if (nodeObservationIter != observations.end()) {
                    NodeCategory &nodeCategory = nodeObservationIter->second;

                    if (nodeCategory == NodeCategory::candidateMinimalDependency) {
                        nodeCategory = observations.updateDependencyCategory(node);
                        if (nodeCategory == NodeCategory::minimalDependency) {
                            minimalDeps.insert(node);
                        }
                    } else if (nodeCategory == NodeCategory::candidateMaximalNonDependency) {
                        nodeCategory = observations.updateNonDependencyCategory(node);
                        if (nodeCategory == NodeCategory::maximalNonDependency) {
                            maximalNonDeps.insert(node);
                        }
                    }

                } else if (!observations.inferCategory(node)) { //TODO переделать inferCategory
                    auto nodePartition = partitionStorage->getOrCreateFor(*node);
                    auto rhsPartition = relation->getColumnData(rhs->getIndex())->getPositionListIndex();
                    auto nodeIntersectedWithRHSPartition = nodePartition->intersect(rhsPartition); //может ещё раз вызвать getOrCreateFor вместо пересечения?

                    if (nodePartition->getNumNonSingletonCluster() ==
                        nodeIntersectedWithRHSPartition->getNumNonSingletonCluster()
                    ) {
                        observations.insert(std::make_pair(*node, NodeCategory::candidateMinimalDependency));
                    } else {
                        observations.insert(std::make_pair(*node, NodeCategory::candidateMaximalNonDependency));
                    }
                }

                node = pickNextNode(node);
            } while (node != nullptr);
        }
        seeds = generateNextSeeds(rhs);
    } while (!seeds.empty());
}

DFD::DFD(const std::filesystem::path &path, char separator, bool hasHeader)
        : FDAlgorithm(path, separator, hasHeader), gen(rd()), observations() {
    relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    partitionStorage = std::make_shared<PartitionStorage>(relation, CachingMethod::ALLCACHING, CacheEvictionMethod::MEDAINUSAGE);
    dependenciesMap = DependenciesMap(relation->getSchema());
    nonDependenciesMap = NonDependenciesMap(relation->getSchema());
}

shared_ptr<Vertical> DFD::pickNextNode(shared_ptr<Vertical> const& node) {
    dynamic_bitset<> columnIndices = node->getColumnIndices();
    auto nodeIter = observations.find(*node);

    //можно зарефакторить, если сделать категорию undefined?
    if (nodeIter != observations.end()) {
        if (nodeIter->second == NodeCategory::candidateMinimalDependency) {
            vector<shared_ptr<Vertical>> uncheckedSubsets = dependenciesMap.getUncheckedSubsets(node, observations); //TODO переписать observations в конструктор?
            if (uncheckedSubsets.empty()) {
                minimalDeps.insert(node);
                //observations[*node] = NodeCategory::minimalDependency;
                dependenciesMap.addNewDependency(node);
            } else {
                shared_ptr<Vertical> nextNode = takeRandom(uncheckedSubsets);
                trace.push(node);
                return nextNode;
            }
        } else if (nodeIter->second == NodeCategory::candidateMaximalNonDependency) {
            vector<shared_ptr<Vertical>> uncheckedSupersets = nonDependenciesMap.getUncheckedSupersets(node, observations);
            if (uncheckedSupersets.empty()) {
                maximalNonDeps.insert(node);
                //observations[*node] = NodeCategory::maximalNonDependency;
                nonDependenciesMap.addNewNonDependency(node);
            } else {
                shared_ptr<Vertical> nextNode = takeRandom(uncheckedSupersets);
                trace.push(node);
                return nextNode;
            }
        }
    }

    shared_ptr<Vertical> nextNode = nullptr;
    if (!trace.empty()) {
        nextNode = trace.top();
        trace.pop();
    }
    return nextNode;
}


std::list<shared_ptr<Vertical>> DFD::generateNextSeeds(shared_ptr<Column const> const& currentRHS) {
    /*struct qwe {
        bool operator()(shared_ptr<Vertical> const& v1, shared_ptr<Vertical> const& v2) const {
            return *v1 == *v2;
        }
    };*/

    using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;

    vertical_set seeds;
    vertical_set newSeeds;
    dynamic_bitset<> singleColumnBitset(relation->getNumColumns(), 0);

    //TODO переписать под метод getColumns
    for (auto const& node : maximalNonDeps) {
        //shared_ptr<Vertical> complementNode = node->invert();
        dynamic_bitset<> nodeIndices = node->getColumnIndices();
        nodeIndices[currentRHS->getIndex()] = true;
        shared_ptr<Vertical> complementNode = std::make_shared<Vertical>(node->getSchema(), std::move(nodeIndices.flip()));
        singleColumnBitset.reset();
        //dynamic_bitset<> complementNodeIndices = complementNode->getColumnIndices();

        if (seeds.empty()) {
            for (size_t columnIndex = complementNode->getColumnIndices().find_first();
                 columnIndex < complementNode->getColumnIndices().size();
                 columnIndex = complementNode->getColumnIndices().find_next(columnIndex)
            ) {
                singleColumnBitset[columnIndex] = true;
                seeds.insert(std::make_shared<Vertical>(relation->getSchema(), singleColumnBitset));
                singleColumnBitset[columnIndex] = false;
            }
        } else {
            for (auto const& dependency : seeds) {
                for (size_t columnIndex = complementNode->getColumnIndices().find_first();
                     columnIndex < complementNode->getColumnIndices().size();
                     columnIndex = complementNode->getColumnIndices().find_next(columnIndex)
                ) {
                    //TODO дикие костыли
                    singleColumnBitset[columnIndex] = true;
                    singleColumnBitset |= dependency->getColumnIndices();
                    newSeeds.insert(std::make_shared<Vertical>(relation->getSchema(), singleColumnBitset));
                    singleColumnBitset.reset();
                }
            }

            minimize(newSeeds);
            seeds.clear();
            for (auto const& newSeed : newSeeds) {
                seeds.insert(newSeed);
            }
            newSeeds.clear();
        }
    }

    //TODO может быть затратно?
    vertical_set discoveredMinimalDepsSet(minimalDeps.begin(), minimalDeps.end());
    /*for (auto const& seed : seeds) {
        if (discoveredMinimalDepsSet.find(seed) != discoveredMinimalDepsSet.end()) {
            seeds.erase(seed);
        }
    }*/

    for (auto seedIter = seeds.begin(); seedIter != seeds.end(); ) {
        if (discoveredMinimalDepsSet.find(*seedIter) != discoveredMinimalDepsSet.end()) {
            seedIter = seeds.erase(seedIter);
        } else {
            seedIter++;
        }
    }

    return std::list(seeds.begin(), seeds.end());
}

//TODO пока что дикий костыль за квадрат
void DFD::minimize(std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator> & nodeList) {
    /*for (auto const& node : nodeList) {
        for (auto const& nodeToCheck : nodeList) {
            if (!(*nodeToCheck == *node) && nodeToCheck->contains(*node)) {
                //TODO переписать на итераторы
                nodeList.erase(nodeToCheck);
            }
        }
    }*/

    for (auto nodeIter = nodeList.begin(); nodeIter != nodeList.end(); nodeIter++) {
        for (auto nodeToCheckIter = nodeList.begin(); nodeToCheckIter != nodeList.end(); ) {
            if ( !(**nodeIter == **nodeToCheckIter) && (*nodeToCheckIter)->contains(**nodeIter)) {
                nodeToCheckIter = nodeList.erase(nodeToCheckIter);
            } else {
                nodeToCheckIter++;
            }
        }
    }
}