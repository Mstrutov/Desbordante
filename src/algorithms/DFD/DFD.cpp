#include "DFD.h"

#include <algorithm>
#include <list>
#include <stack>

#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "PositionListIndex.h"

unsigned long long DFD::execute() {
    RelationalSchema const* const schema = relation->getSchema();

    //обработка случая, когда пустая таблица
    if (relation->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }

    auto startTime = std::chrono::system_clock::now();

    std::vector<Vertical> uniqueVerticals;

    //ищем уникальные столбцы
    for (auto const& column : schema->getColumns()) {
        ColumnData& columnData = relation->getColumnData(column->getIndex());
        PositionListIndex const* const columnPLI = columnData.getPositionListIndex();

        if (columnPLI->getNumNonSingletonCluster() == 0) {
            Vertical const lhs = Vertical(*column);
            uniqueVerticals.push_back(lhs);
            //в метаноме регистрируем зависимость сразу, а тут нет, чтобы сначала рассмотреть пустые
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
        trace = std::stack<Vertical>();

        ColumnData const& rhsData = relation->getColumnData(rhs->getIndex());
        PositionListIndex const* const rhsPLI = rhsData.getPositionListIndex();

        //если все строки имеют одинаковое значение, то добавляем зависимость с пустым lhs
        if (rhsPLI->getNepAsLong() == relation->getNumTuplePairs()) {
            this->registerFD(*(schema->emptyVertical), *rhs);
            continue;
            //минимальная зависимость []->RHS, меньше точно не найти, поэтому берем следующий RHS
        }

        //в метаноме немного по-другому, но суть такая же
        //тут обрабатываем найденные уникальные столбцы
        for (auto const& lhs : uniqueVerticals) {
            if (!lhs.contains(*rhs)) {
                observations[lhs] = NodeCategory::minimalDependency;
                dependenciesMap.addNewDependency(lhs);
                minimalDeps.insert(lhs); //вот теперь добавляем
            }
        }

        //находим минимальные зависимости для текущего RHS
        findLHSs(rhs.get());

        //регистрируем полученные зависимости
        for (auto const& minimalDependencyLHS : minimalDeps) {
            registerFD(minimalDependencyLHS, *rhs);
        }
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    long long aprioriMillis = elapsed_milliseconds.count();

    //можно вывести найденные зависимости в формате Json:
    //std::cout << "====JSON-FD========\r\n" << getJsonFDs() << std::endl;
    std::cout << "HASH: " << FDAlgorithm::fletcher16() << std::endl;

    return aprioriMillis;
}

Vertical const& DFD::takeRandom(std::unordered_set<Vertical> & nodeSet) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeSet.begin(), nodeSet.end()) - 1);
    auto iterator = nodeSet.begin();
    std::advance(iterator, dis(this->gen));
    Vertical const& node = *iterator;
    return node;
}

void DFD::findLHSs(Column const* const  rhs) {
    RelationalSchema const* const schema = relation->getSchema();
    std::stack<Vertical> seeds;

    for (int partitionIndex : columnOrder.getOrderHighDistinctCount(Vertical(*rhs).invert())) {
        if (partitionIndex != rhs->getIndex()) {
            seeds.push(Vertical(*schema->getColumn(partitionIndex)));
        }
    }

    do {
        while (!seeds.empty()) {
            Vertical node;
            if (!seeds.empty()) {
                node = std::move(seeds.top());
                seeds.pop();
            } else {
                node = *schema->emptyVertical;
            }

            do {
                auto const nodeObservationIter = observations.find(node);

                if (nodeObservationIter != observations.end()) {
                    NodeCategory& nodeCategory = nodeObservationIter->second;

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
                } else if (!inferCategory(node, rhs->getIndex())) {
                    //не смогли определить категорию --- значит считаем партиции
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
                        observations.updateDependencyCategory(node);
                        if (observations[node] == NodeCategory::minimalDependency) {
                            minimalDeps.insert(node);
                        }
                        dependenciesMap.addNewDependency(node);
                    } else {
                        observations.updateNonDependencyCategory(node, rhs->getIndex());
                        if (observations[node] == NodeCategory::maximalNonDependency) {
                            maximalNonDeps.insert(node);
                        }
                        nonDependenciesMap.addNewNonDependency(node);
                    }
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
    partitionStorage = std::make_unique<PartitionStorage>(relation.get(), CachingMethod::ALLCACHING, CacheEvictionMethod::MEDAINUSAGE);
    dependenciesMap = DependenciesMap(relation->getSchema());
    nonDependenciesMap = NonDependenciesMap(relation->getSchema());
    columnOrder = ColumnOrder(relation.get());
}

Vertical DFD::pickNextNode(Vertical const &node, size_t rhsIndex) {
    auto nodeIter = observations.find(node);

    if (nodeIter != observations.end()) {
        if (nodeIter->second == NodeCategory::candidateMinimalDependency) {
            auto uncheckedSubsets = observations.getUncheckedSubsets(node, columnOrder);
            auto prunedNonDepSubsets = nonDependenciesMap.getPrunedSupersets(uncheckedSubsets);
            for (auto const& prunedSubset : prunedNonDepSubsets) {
                observations[prunedSubset] = NodeCategory::nonDependency;
            }
            substractSets(uncheckedSubsets, prunedNonDepSubsets);

            if (uncheckedSubsets.empty() && prunedNonDepSubsets.empty()) {
                minimalDeps.insert(node);
                observations[node] = NodeCategory::minimalDependency;
            } else if (!uncheckedSubsets.empty()) {
                auto const& nextNode = takeRandom(uncheckedSubsets);
                //чтобы при каждом запуске выбирать одинаковый путь, нужно заменить строчку выше на строку ниже
                //auto const& nextNode = *uncheckedSubsets.begin();
                trace.push(node);
                return nextNode;
            }
        } else if (nodeIter->second == NodeCategory::candidateMaximalNonDependency) {
            auto uncheckedSupersets = observations.getUncheckedSupersets(node, rhsIndex, columnOrder);
            auto prunedNonDepSupersets = nonDependenciesMap.getPrunedSupersets(uncheckedSupersets);
            auto prunedDepSupersets = dependenciesMap.getPrunedSubsets(uncheckedSupersets);

            for (auto const& prunedSuperset : prunedNonDepSupersets) {
                observations[prunedSuperset] = NodeCategory::nonDependency;
            }
            for (auto const& prunedSuperset : prunedDepSupersets) {
                observations[prunedSuperset] = NodeCategory::dependency;
            }

            substractSets(uncheckedSupersets, prunedDepSupersets);
            substractSets(uncheckedSupersets, prunedNonDepSupersets);

            if (uncheckedSupersets.empty() && prunedNonDepSupersets.empty()) {
                maximalNonDeps.insert(node);
                observations[node] = NodeCategory::maximalNonDependency;
            } else if (!uncheckedSupersets.empty()) {
                auto const& nextNode = takeRandom(uncheckedSupersets);
                //чтобы при каждом запуске выбирать одинаковый путь, нужно заменить строчку выше на строку ниже
                //auto const& nextNode = *uncheckedSupersets.begin();
                trace.push(node);
                return nextNode;
            }
        }
    }

    Vertical nextNode = *(node.getSchema()->emptyVertical);
    if (!trace.empty()) {
        nextNode = trace.top();
        trace.pop();
    }
    return nextNode;
}

std::stack<Vertical> DFD::generateNextSeeds(Column const* const currentRHS) {
    std::unordered_set<Vertical> seeds;
    std::unordered_set<Vertical> newSeeds;

    for (auto const& nonDep : maximalNonDeps) {
        auto complementIndices = nonDep.getColumnIndicesRef();
        complementIndices[currentRHS->getIndex()] = true;
        complementIndices.flip();

        if (seeds.empty()) {
            boost::dynamic_bitset<> singleColumnBitset(relation->getNumColumns(), 0);
            singleColumnBitset.reset();

            for (size_t columnIndex = complementIndices.find_first();
                 columnIndex < complementIndices.size();
                 columnIndex = complementIndices.find_next(columnIndex)
            ) {
                singleColumnBitset[columnIndex] = true;
                seeds.emplace(relation->getSchema(), singleColumnBitset);
                singleColumnBitset[columnIndex] = false;
            }
        } else {
            for (auto const& dependency : seeds) {
                auto newCombination = dependency.getColumnIndicesRef();

                for (size_t columnIndex = complementIndices.find_first();
                     columnIndex < complementIndices.size();
                     columnIndex = complementIndices.find_next(columnIndex)
                ) {
                    newCombination[columnIndex] = true;
                    newSeeds.emplace(relation->getSchema(), newCombination);
                    newCombination[columnIndex] = dependency.getColumnIndicesRef()[columnIndex];
                }
            }

            std::list<Vertical> minimizedNewSeeds = minimize(newSeeds);
            seeds.clear();
            for (auto & newSeed : minimizedNewSeeds) {
                seeds.insert(std::move(newSeed));
            }
            newSeeds.clear();
        }
    }

    for (auto seedIter = seeds.begin(); seedIter != seeds.end(); ) {
        if (minimalDeps.find(*seedIter) != minimalDeps.end()) {
            seedIter = seeds.erase(seedIter);
        } else {
            seedIter++;
        }
    }

    std::stack<Vertical> remainingSeeds;

    for (auto const& newSeed : seeds) {
        remainingSeeds.push(newSeed);
    }

    return remainingSeeds;
}

std::list<Vertical> DFD::minimize(std::unordered_set<Vertical> const& nodeList) {
    long long maxCardinality = 0;
    std::unordered_map<long long, std::list<Vertical const*>> seedsBySize(nodeList.size() / relation->getNumColumns());
    for (auto const& seed : nodeList) {
        long long cardinalityOfSeed = seed.getArity();
        maxCardinality = std::max(maxCardinality, cardinalityOfSeed);
        if (seedsBySize.find(cardinalityOfSeed) == seedsBySize.end()) {
            seedsBySize[cardinalityOfSeed] = std::list<Vertical const*>();
        }
        seedsBySize[cardinalityOfSeed].push_back(&seed);
    }

    for (long long lowerBound = 1; lowerBound < maxCardinality; lowerBound++) {
        if (seedsBySize.find(lowerBound) != seedsBySize.end()) {
            auto const& lowerBoundSeeds = seedsBySize.find(lowerBound)->second;
            for (long long upperBound = maxCardinality; upperBound > lowerBound; upperBound--) {
                if (seedsBySize.find(upperBound) != seedsBySize.end()) {
                    auto& upperBoundSeeds = seedsBySize.find(upperBound)->second;
                    for (auto const& lowerSeed : lowerBoundSeeds) {
                        for (auto upperIt = upperBoundSeeds.begin(); upperIt != upperBoundSeeds.end();) {
                            if ((*upperIt)->contains(*lowerSeed)) {
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

    std::list<Vertical> newSeeds;
    for (auto & seedList : seedsBySize) {
        for (auto& seed : seedList.second) {
            newSeeds.push_back(*seed);
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
        observations.updateNonDependencyCategory(node, rhsIndex);
        nonDependenciesMap.addNewNonDependency(node);
        if (observations[node] == NodeCategory::minimalDependency) {
            minimalDeps.insert(node);
        }
        return true;
    } else if (dependenciesMap.canBePruned(node)) {
        observations.updateDependencyCategory(node);
        dependenciesMap.addNewDependency(node);
        if (observations[node] == NodeCategory::maximalNonDependency) {
            maximalNonDeps.insert(node);
        }
        return true;
    }

    return false;
}
