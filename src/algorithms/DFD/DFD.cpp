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

//using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;

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
            Vertical const lhs = Vertical(*column);
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
        trace = std::stack<shared_ptr<Vertical>>();//TODO clear trace?

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

shared_ptr<Vertical> DFD::takeRandom(vertical_set & nodeSet) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeSet.begin(), nodeSet.end()) - 1);
    auto iterator = nodeSet.begin();
    std::advance(iterator, dis(this->gen));
    shared_ptr<Vertical> node = move(*iterator);
    nodeSet.erase(iterator);
    return node;
}

shared_ptr<Vertical> DFD::takeRandom(std::vector<shared_ptr<Vertical>> const& nodeVector) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeVector.begin(), nodeVector.end()) - 1);
    auto iterator = nodeVector.begin();
    std::advance(iterator, dis(this->gen));
    return *iterator;
}

void DFD::findLHSs(shared_ptr<Column const> const& rhs, shared_ptr<RelationalSchema> schema) {

    std::list<shared_ptr<Vertical>> seeds; //TODO лист или вектор?

    //initialize seeds nodes
    for (auto const& column : schema->getColumns()) {
        //add to seeds all columns except rhs
        //seeds.push_back(schema->emptyVertical);
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
                        nodeCategory = observations.updateNonDependencyCategory(node, rhs->getIndex());
                        if (nodeCategory == NodeCategory::maximalNonDependency) {
                            maximalNonDeps.insert(node);
                        }
                    }
                    //node = pickNextNode(node);
                } else if (!inferCategory(node, rhs->getIndex())) { //TODO переделать inferCategory
                    auto nodePartition = partitionStorage->getOrCreateFor(*node);
                    //auto rhsPartition = relation->getColumnData(rhs->getIndex())->getPositionListIndex();
                    //auto nodeIntersectedWithRHSPartition = nodePartition->intersect(rhsPartition); //может ещё раз вызвать getOrCreateFor вместо пересечения?
                    auto nodeIntersectedWithRHSPartition = partitionStorage->getOrCreateFor(*(node->Union(*rhs)));

                    if (nodePartition->getNumNonSingletonCluster() ==
                        nodeIntersectedWithRHSPartition->getNumNonSingletonCluster()
                    ) {
                        //observations.insert(std::make_pair(*node, NodeCategory::candidateMinimalDependency));
                        observations[*node] = observations.updateDependencyCategory(node);
                        if (observations[*node] == NodeCategory::minimalDependency) {
                            minimalDeps.insert(node);
                        }
                        dependenciesMap.addNewDependency(node);
                    } else {
                        //observations.insert(std::make_pair(*node, NodeCategory::candidateMaximalNonDependency));
                        observations[*node] = observations.updateNonDependencyCategory(node, rhs->getIndex());
                        if (observations[*node] == NodeCategory::maximalNonDependency) {
                            maximalNonDeps.insert(node);
                        }
                        nonDependenciesMap.addNewNonDependency(node);
                    }
                    //node = pickNextNode(node);
                }

                node = pickNextNode(node, rhs->getIndex());
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

shared_ptr<Vertical> DFD::pickNextNode(shared_ptr<Vertical> const &node, size_t rhsIndex) {
    dynamic_bitset<> columnIndices = node->getColumnIndices();
    auto nodeIter = observations.find(*node);

    //можно зарефакторить, если сделать категорию undefined?
    if (nodeIter != observations.end()) {
        if (nodeIter->second == NodeCategory::candidateMinimalDependency) {
            vertical_set uncheckedSubsets = observations.getUncheckedSubsets(node, rhsIndex); //TODO переписать observations в конструктор?
            vertical_set prunedNonDepSubsets = nonDependenciesMap.getPrunedSupersets(uncheckedSubsets);
            for (auto const& prunedSubset : prunedNonDepSubsets) {
                observations[*prunedSubset] = NodeCategory::nonDependency;
                //dependenciesMap.addNewDependency(node);
            }
            substractSets(uncheckedSubsets, prunedNonDepSubsets);

            if (uncheckedSubsets.empty()) { //TODO в метаноме немнго по другому
                minimalDeps.insert(node);
                observations[*node] = NodeCategory::minimalDependency;
                //dependenciesMap.addNewDependency(node);
            } else {
                shared_ptr<Vertical> nextNode = takeRandom(uncheckedSubsets);
                trace.push(node);
                return nextNode;
            }
        } else if (nodeIter->second == NodeCategory::candidateMaximalNonDependency) {
            vertical_set uncheckedSupersets = observations.getUncheckedSupersets(node, rhsIndex);
            vertical_set prunedNonDepSupersets = nonDependenciesMap.getPrunedSupersets(uncheckedSupersets);
            vertical_set prunedDepSupersets = dependenciesMap.getPrunedSubsets(uncheckedSupersets);

            for (auto const& prunedSuperset : prunedNonDepSupersets) {
                observations[*prunedSuperset] = NodeCategory::nonDependency;
                //nonDependenciesMap.addNewNonDependency(node);
            }
            for (auto const& prunedSuperset : prunedDepSupersets) {
                observations[*prunedSuperset] = NodeCategory::dependency;
                //dependenciesMap.addNewDependency(node);
            }

            substractSets(uncheckedSupersets, prunedDepSupersets);
            substractSets(uncheckedSupersets, prunedNonDepSupersets);

            if (uncheckedSupersets.empty()) {
                maximalNonDeps.insert(node);
                observations[*node] = NodeCategory::maximalNonDependency;
                //nonDependenciesMap.addNewNonDependency(node);
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

    vertical_set seeds;
    vertical_set newSeeds;
    dynamic_bitset<> singleColumnBitset(relation->getNumColumns(), 0);

    //TODO переписать под метод getColumns
    for (auto const& nonDep : maximalNonDeps) {
        //shared_ptr<Vertical> complementNode = nonDep->invert();
        dynamic_bitset<> nodeIndices = nonDep->getColumnIndices();
        nodeIndices[currentRHS->getIndex()] = true;
        shared_ptr<Vertical> complementNode = std::make_shared<Vertical>(nonDep->getSchema(), std::move(nodeIndices.flip()));
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

void DFD::substractSets(vertical_set & set, vertical_set const& setToSubstract) {
    for (const auto & nodeToDelete : setToSubstract) {
        auto foundElementIter = set.find(nodeToDelete);
        if (foundElementIter != set.end()) {
            set.erase(foundElementIter);
        }
    }
}

bool DFD::inferCategory(shared_ptr<Vertical> const &node, size_t rhsIndex) {
    if (nonDependenciesMap.canBePruned(*node)) {
        observations[*node] = observations.updateNonDependencyCategory(node, rhsIndex);
        nonDependenciesMap.addNewNonDependency(node);
        if (observations[*node] == NodeCategory::minimalDependency) { //хз надо или нет
            minimalDeps.insert(node);
        }
        return true;
    } else if (dependenciesMap.canBePruned(*node)) {
        observations[*node] = observations.updateDependencyCategory(node);
        dependenciesMap.addNewDependency(node);
        if (observations[*node] == NodeCategory::maximalNonDependency) { //хз надо или нет
            maximalNonDeps.insert(node);
        }
        return true;
    }

    return false;
}