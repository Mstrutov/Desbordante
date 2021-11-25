#include "LatticeTraversal.h"

#include <random>

#include "PositionListIndex.h"

LatticeTraversal::LatticeTraversal(const Column *const rhs,
                                   const ColumnLayoutRelationData *const relation,
                                   const std::vector<Vertical> &uniqueVerticals,
                                   PartitionStorage *const partitionStorage)
    : rhs(rhs), relation(relation), uniqueColumns(uniqueVerticals), gen(rd()),
      dependenciesMap(relation->getSchema()), nonDependenciesMap(relation->getSchema()),
      columnOrder(relation), partitionStorage(partitionStorage)
    {}

std::unordered_set<Vertical> LatticeTraversal::findLHSs() {
    RelationalSchema const* const schema = relation->getSchema();

    //processing of found unique columns
    for (auto const &lhs: uniqueColumns) {
        if (!lhs.contains(*rhs)) {
            observations[lhs] = NodeCategory::minimalDependency;
            dependenciesMap.addNewDependency(lhs);
            minimalDeps.insert(lhs);
        }
    }

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
                    //if we were not able to infer category, we calculate the partitions
                    auto nodePLI = partitionStorage->getOrCreateFor(node);
                    auto nodePliPointer = std::holds_alternative<util::PositionListIndex*>(nodePLI)
                                          ? std::get<util::PositionListIndex*>(nodePLI)
                                          : std::get<std::unique_ptr<util::PositionListIndex>>(nodePLI).get();
                    auto intersectedPLI = partitionStorage->getOrCreateFor(node.Union(*rhs));
                    auto intersectedPLIPointer = std::holds_alternative<util::PositionListIndex*>(intersectedPLI)
                                                 ? std::get<util::PositionListIndex*>(intersectedPLI)
                                                 : std::get<std::unique_ptr<util::PositionListIndex>>(intersectedPLI).get();

                    if (nodePliPointer->getNepAsLong() ==
                        intersectedPLIPointer->getNepAsLong()
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

    return minimalDeps;
}

bool LatticeTraversal::inferCategory(Vertical const& node, unsigned int rhsIndex) {
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

Vertical const& LatticeTraversal::takeRandom(std::unordered_set<Vertical> & nodeSet) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeSet.begin(), nodeSet.end()) - 1);
    auto iterator = nodeSet.begin();
    std::advance(iterator, dis(this->gen));
    Vertical const& node = *iterator;
    return node;
}

Vertical LatticeTraversal::pickNextNode(Vertical const &node, unsigned int rhsIndex) {
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

std::stack<Vertical> LatticeTraversal::generateNextSeeds(Column const* const currentRHS) {
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

std::list<Vertical> LatticeTraversal::minimize(std::unordered_set<Vertical> const& nodeList) const {
    unsigned int maxCardinality = 0;
    std::unordered_map<unsigned int, std::list<Vertical const*>> seedsBySize(nodeList.size() / relation->getNumColumns());

    for (auto const& seed : nodeList) {
        unsigned int const cardinalityOfSeed = seed.getArity();
        maxCardinality = std::max(maxCardinality, cardinalityOfSeed);
        if (seedsBySize.find(cardinalityOfSeed) == seedsBySize.end()) {
            seedsBySize[cardinalityOfSeed] = std::list<Vertical const*>();
        }
        seedsBySize[cardinalityOfSeed].push_back(&seed);
    }

    for (unsigned int lowerBound = 1; lowerBound < maxCardinality; ++lowerBound) {
        if (seedsBySize.find(lowerBound) != seedsBySize.end()) {
            auto const& lowerBoundSeeds = seedsBySize.find(lowerBound)->second;
            for (unsigned int upperBound = maxCardinality; upperBound > lowerBound; --upperBound) {
                if (seedsBySize.find(upperBound) != seedsBySize.end()) {
                    auto& upperBoundSeeds = seedsBySize.find(upperBound)->second;
                    for (auto const& lowerSeed : lowerBoundSeeds) {
                        for (auto upperIt = upperBoundSeeds.begin(); upperIt != upperBoundSeeds.end();) {
                            if ((*upperIt)->contains(*lowerSeed)) {
                                upperIt = upperBoundSeeds.erase(upperIt);
                            } else {
                                ++upperIt;
                            }
                        }
                    }
                }
            }
        }
    }

    std::list<Vertical> newSeeds;
    for (auto & seedList : seedsBySize) {
        for (auto & seed : seedList.second) {
            newSeeds.push_back(*seed);
        }
    }
    return newSeeds;
}

void LatticeTraversal::substractSets(std::unordered_set<Vertical> & set, std::unordered_set<Vertical> const& setToSubstract) {
    for (const auto & nodeToDelete : setToSubstract) {
        auto foundElementIter = set.find(nodeToDelete);
        if (foundElementIter != set.end()) {
            set.erase(foundElementIter);
        }
    }
}
