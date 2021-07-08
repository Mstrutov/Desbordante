//
// Created by alexandrsmirn
//

#include "DFD.h"

#include <algorithm>
#include <random>
#include <list>

#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "PositionListIndex.h"
#include "PLICache.h"

//using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;

unsigned long long DFD::execute() {
    //shared_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);//second parameter?
    RelationalSchema const* const schema = relation->getSchema();
    if (relation->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }

    auto startTime = std::chrono::system_clock::now();

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
        ColumnData& columnData = relation->getColumnData(column->getIndex());
        PositionListIndex const* const columnPLI = columnData.getPositionListIndex();

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
    for (auto & rhs : schema->getColumns()) {
        //тут строим новую решетку, соответственно нужно обнулить/завести некоторые структуры данных
        minimalDeps.clear();
        maximalNonDeps.clear();
        dependenciesMap = DependenciesMap(relation->getSchema());
        nonDependenciesMap = NonDependenciesMap(relation->getSchema());
        observations.clear();
        trace = std::stack<Vertical>();//TODO clear trace?

        ColumnData const& rhsData = relation->getColumnData(rhs->getIndex());
        PositionListIndex const* const rhsPLI = rhsData.getPositionListIndex();

        if (rhsPLI->getNepAsLong() == relation->getNumTuplePairs()) {
            this->registerFD(*(schema->emptyVertical), *rhs);
            continue;
        }

        //в метаноме немного по другому
        for (auto const& lhs : uniqueVerticals) {
            if (!lhs.contains(*rhs)) {
                observations[lhs] = NodeCategory::minimalDependency;
                //разобраться с шаред поинтерами
                dependenciesMap.addNewDependency(Vertical(lhs));//TODO костыль
                minimalDeps.insert(lhs);
            }
        }

        findLHSs(rhs.get(), schema);

        //TODO не заьыть зареать все зависимости из списка
        for (auto const& minimalDependencyLHS : minimalDeps) {
            registerFD(minimalDependencyLHS, *rhs);
        }
    }
    std::chrono::milliseconds elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    long long aprioriMillis = elapsed_milliseconds.count();

    std::cout << "====JSON-FD========\r\n" << getJsonFDs() << std::endl;
    std::cout << "HASH: " << FDAlgorithm::fletcher16() << std::endl;

    return aprioriMillis;
}

Vertical DFD::takeRandom(std::list<Vertical> & nodeList) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeList.begin(), nodeList.end()) - 1);
    auto iterator = nodeList.begin();
    std::advance(iterator, dis(this->gen));
    Vertical node = std::move(*iterator);
    nodeList.erase(iterator);
    return node;
}

Vertical DFD::takeRandom(std::unordered_set<Vertical> & nodeSet) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeSet.begin(), nodeSet.end()) - 1);
    auto iterator = nodeSet.begin();
    std::advance(iterator, dis(this->gen));
    Vertical node = *iterator;
    nodeSet.erase(iterator);
    return node;
}

Vertical DFD::takeRandom(std::vector<Vertical> const& nodeVector) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeVector.begin(), nodeVector.end()) - 1);
    auto iterator = nodeVector.begin();
    std::advance(iterator, dis(this->gen));
    return *iterator;
}

void DFD::findLHSs(Column const* const  rhs, RelationalSchema const* const schema) {

    std::list<Vertical> seeds; //TODO лист или вектор?

    //initialize seeds nodes
    for (auto const& column : schema->getColumns()) {
        //add to seeds all columns except rhs
        //seeds.push_back(schema->emptyVertical);
        if (column->getIndex() != rhs->getIndex())
            seeds.push_back(Vertical(*column)); //TODO проверить, поменял emplace на push
    }

    //auto

    do {
        while (!seeds.empty()) {
            Vertical node = takeRandom(seeds);
            do {
                auto const nodeObservationIter = observations.find(node); //const?

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
                    /*auto nodePartition = partitionStorage->getOrCreateFor(node);
                    auto rhsPartition = relation->getColumnData(rhs->getIndex()).getPositionListIndex();
                    auto nodeIntersectedWithRHSPartition = nodePartition->intersect(rhsPartition); //может ещё раз вызвать getOrCreateFor вместо пересечения?
                    //auto nodeIntersectedWithRHSPartition = partitionStorage->getOrCreateFor(*(node->Union(*rhs)));
                     */

                    auto nodePLI = partitionStorage->getOrCreateFor(node);
                    auto nodePliPointer = std::holds_alternative<PositionListIndex*>(nodePLI)
                                      ? std::get<PositionListIndex*>(nodePLI)
                                      : std::get<std::unique_ptr<PositionListIndex>>(nodePLI).get();
                    auto intersectedPLI = partitionStorage->getOrCreateFor(node.Union(*rhs));
                    auto intersectrdPLIPointer = std::holds_alternative<PositionListIndex*>(intersectedPLI)
                                      ? std::get<PositionListIndex*>(intersectedPLI)
                                      : std::get<std::unique_ptr<PositionListIndex>>(intersectedPLI).get();

                    if (nodePliPointer->getNepAsLong() ==
                        intersectrdPLIPointer->getNepAsLong()
                    ) {
                        //observations.insert(std::make_pair(*node, NodeCategory::candidateMinimalDependency));
                        observations[node] = observations.updateDependencyCategory(node);
                        if (observations[node] == NodeCategory::minimalDependency) {
                            minimalDeps.insert(node);
                        }
                        dependenciesMap.addNewDependency(node);
                    } else {
                        //observations.insert(std::make_pair(*node, NodeCategory::candidateMaximalNonDependency));
                        observations[node] = observations.updateNonDependencyCategory(node, rhs->getIndex());
                        if (observations[node] == NodeCategory::maximalNonDependency) {
                            maximalNonDeps.insert(node);
                        }
                        nonDependenciesMap.addNewNonDependency(node);
                    }
                    //node = pickNextNode(node);
                }

                node = pickNextNode(node, rhs->getIndex());
            } while (node != *node.getSchema()->emptyVertical);
        }
        seeds = generateNextSeeds(rhs);
    } while (!seeds.empty());
}

DFD::DFD(const std::filesystem::path &path, char separator, bool hasHeader)
        : FDAlgorithm(path, separator, hasHeader), gen(rd()), observations() {
    relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    //TODO прааильно ли передаем uniqur_ptr?
    partitionStorage = std::make_unique<PartitionStorage>(relation.get(), CachingMethod::ALLCACHING, CacheEvictionMethod::MEDAINUSAGE);
    dependenciesMap = DependenciesMap(relation->getSchema());
    nonDependenciesMap = NonDependenciesMap(relation->getSchema());
}

Vertical DFD::pickNextNode(Vertical const &node, size_t rhsIndex) {
    boost::dynamic_bitset<> columnIndices = node.getColumnIndices();
    auto nodeIter = observations.find(node);

    //можно зарефакторить, если сделать категорию undefined?
    if (nodeIter != observations.end()) {
        if (nodeIter->second == NodeCategory::candidateMinimalDependency) {
            std::unordered_set<Vertical> uncheckedSubsets = observations.getUncheckedSubsets(node, rhsIndex); //TODO переписать observations в конструктор?
            std::unordered_set<Vertical> prunedNonDepSubsets = nonDependenciesMap.getPrunedSupersets(uncheckedSubsets);
            for (auto const& prunedSubset : prunedNonDepSubsets) {
                observations[prunedSubset] = NodeCategory::nonDependency;
                //dependenciesMap.addNewDependency(node);
            }
            substractSets(uncheckedSubsets, prunedNonDepSubsets);

            if (uncheckedSubsets.empty()) { //TODO в метаноме немнго по другому
                minimalDeps.insert(node);
                observations[node] = NodeCategory::minimalDependency;
                //dependenciesMap.addNewDependency(node);
            } else {
                Vertical nextNode = takeRandom(uncheckedSubsets);
                trace.push(node);
                return nextNode;
            }
        } else if (nodeIter->second == NodeCategory::candidateMaximalNonDependency) {
            std::unordered_set<Vertical> uncheckedSupersets = observations.getUncheckedSupersets(node, rhsIndex);
            std::unordered_set<Vertical> prunedNonDepSupersets = nonDependenciesMap.getPrunedSupersets(uncheckedSupersets);
            std::unordered_set<Vertical> prunedDepSupersets = dependenciesMap.getPrunedSubsets(uncheckedSupersets);

            for (auto const& prunedSuperset : prunedNonDepSupersets) {
                observations[prunedSuperset] = NodeCategory::nonDependency;
                //nonDependenciesMap.addNewNonDependency(node);
            }
            for (auto const& prunedSuperset : prunedDepSupersets) {
                observations[prunedSuperset] = NodeCategory::dependency;
                //dependenciesMap.addNewDependency(node);
            }

            substractSets(uncheckedSupersets, prunedDepSupersets);
            substractSets(uncheckedSupersets, prunedNonDepSupersets);

            if (uncheckedSupersets.empty()) {
                maximalNonDeps.insert(node);
                observations[node] = NodeCategory::maximalNonDependency;
                //nonDependenciesMap.addNewNonDependency(node);
            } else {
                    Vertical nextNode = takeRandom(uncheckedSupersets);
                trace.push(node);
                return nextNode;
            }
        }
    }

    //TODO тут проверить
    Vertical nextNode = *(node.getSchema()->emptyVertical);
    if (!trace.empty()) {
        nextNode = trace.top();
        trace.pop();
    }
    return nextNode;
}


std::list<Vertical> DFD::generateNextSeeds(Column const* const currentRHS) {

    std::unordered_set<Vertical> seeds;
    std::unordered_set<Vertical> newSeeds;
    boost::dynamic_bitset<> singleColumnBitset(relation->getNumColumns(), 0);

    //TODO переписать под метод getColumns
    for (auto const& nonDep : maximalNonDeps) {
        //shared_ptr<Vertical> complementNode = nonDep->invert();
        boost::dynamic_bitset<> nodeIndices = nonDep.getColumnIndices();
        nodeIndices[currentRHS->getIndex()] = true;
        Vertical complementNode = Vertical(nonDep.getSchema(), std::move(nodeIndices.flip()));
        singleColumnBitset.reset();
        //dynamic_bitset<> complementNodeIndices = complementNode->getColumnIndices();

        if (seeds.empty()) {
            for (size_t columnIndex = complementNode.getColumnIndices().find_first();
                 columnIndex < complementNode.getColumnIndices().size();
                 columnIndex = complementNode.getColumnIndices().find_next(columnIndex)
            ) {
                singleColumnBitset[columnIndex] = true;
                seeds.insert(Vertical(relation->getSchema(), singleColumnBitset));
                singleColumnBitset[columnIndex] = false;
            }
        } else {
            for (auto const& dependency : seeds) {
                for (size_t columnIndex = complementNode.getColumnIndices().find_first();
                     columnIndex < complementNode.getColumnIndices().size();
                     columnIndex = complementNode.getColumnIndices().find_next(columnIndex)
                ) {
                    //TODO дикие костыли
                    singleColumnBitset[columnIndex] = true;
                    singleColumnBitset |= dependency.getColumnIndices();
                    newSeeds.insert(Vertical(relation->getSchema(), singleColumnBitset));
                    singleColumnBitset.reset();
                }
            }

            std::vector<Vertical> minimizedNewSeeds = minimize(newSeeds);
            seeds.clear();
            for (auto & newSeed : minimizedNewSeeds) {
                seeds.insert(std::move(newSeed));
            }
            newSeeds.clear();
        }
    }

    //TODO может быть затратно?
    std::unordered_set<Vertical> discoveredMinimalDepsSet(minimalDeps.begin(), minimalDeps.end());

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
std::vector<Vertical> DFD::minimize(std::unordered_set<Vertical> const& nodeList) {
    long long maxCardinality = 0;
    std::unordered_map<long long, std::list<Vertical>> seedsBySize;
    for (auto const& seed : nodeList) {
        long long cardinalityOfSeed = seed.getArity(); //TODO check
        maxCardinality = std::max(maxCardinality, cardinalityOfSeed);
        if (seedsBySize.find(cardinalityOfSeed) == seedsBySize.end()) {
            seedsBySize[cardinalityOfSeed] = std::list<Vertical>();
        }
        seedsBySize[cardinalityOfSeed].push_back(seed);
    }

    for (long long lowerBound = 1; lowerBound < maxCardinality; lowerBound++) {
        if (seedsBySize.find(lowerBound) != seedsBySize.end()) {
            std::list<Vertical> const& lowerBoundSeeds = seedsBySize.find(lowerBound)->second;
            for (long long upperBound = maxCardinality; upperBound > lowerBound; upperBound--) {
                if (seedsBySize.find(upperBound) != seedsBySize.end()) {
                    std::list<Vertical> & upperBoundSeeds = seedsBySize.find(upperBound)->second;
                    for (auto lowerIt = lowerBoundSeeds.begin(); lowerIt != lowerBoundSeeds.end(); lowerIt++) {
                        //Vertical const& lowerSeed = *lowerIt;
                        for (auto upperIt = upperBoundSeeds.begin(); upperIt != upperBoundSeeds.end();) {
                            if (upperIt->contains(*lowerIt)) {
                                upperIt = upperBoundSeeds.erase(upperIt);
                            } else {
                                upperIt++;
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<Vertical> newSeeds();
    for (auto seedList : seedsBySize) {
        for (Vertical seed : seedList.second) {
            newSeeds.push_back(std::move(seed));
        }
    }
    return newSeeds;
}


void DFD::substractSets(std::unordered_set<Vertical> & set, std::unordered_set<Vertical> const& setToSubstract) {
    for (const auto & nodeToDelete : setToSubstract) {
        auto foundElementIter = set.find(nodeToDelete);
        if (foundElementIter != set.end()) {
            set.erase(foundElementIter);
        }
    }
}

bool DFD::inferCategory(Vertical const& node, size_t rhsIndex) {
    if (nonDependenciesMap.canBePruned(node)) {
        observations[node] = observations.updateNonDependencyCategory(node, rhsIndex);
        nonDependenciesMap.addNewNonDependency(node);
        if (observations[node] == NodeCategory::minimalDependency) { //хз надо или нет
            minimalDeps.insert(node);
        }
        return true;
    } else if (dependenciesMap.canBePruned(node)) {
        observations[node] = observations.updateDependencyCategory(node);
        dependenciesMap.addNewDependency(node);
        if (observations[node] == NodeCategory::maximalNonDependency) { //хз надо или нет
            maximalNonDeps.insert(node);
        }
        return true;
    }

    return false;
}