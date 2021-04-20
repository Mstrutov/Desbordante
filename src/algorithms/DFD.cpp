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

//using std::shared_ptr;

unsigned long long DFD::execute() {
    shared_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);//second parameter?
    shared_ptr<RelationalSchema> schema = relation->getSchema();
    partitionCache = std::make_shared<PLICache>(relation, CachingMethod::ALLCACHING, CacheEvictionMethod::MEDAINUSAGE, 0, 0, 0, 0, 0, 0, 0);

    std::list<shared_ptr<Column>> possibleRHSs(schema->getColumns().begin(), schema->getColumns().end());

    //first loop of DFD
    for (auto columnIter = possibleRHSs.begin(); columnIter != possibleRHSs.end(); columnIter++) {
        shared_ptr<ColumnData> columnData = relation->getColumnData((*columnIter)->getIndex());
        shared_ptr<PositionListIndex> columnPLI = columnData->getPositionListIndex();

        //if current column is unique
        if (columnPLI->getNumNonSingletonCluster() == 0) {
            possibleRHSs.erase(columnIter);
            auto lhs = Vertical(**columnIter);//наверно стоит убрать, потому что при первом запуске же будет мув

            for (auto rhs : possibleRHSs) {
                this->registerFD(lhs, *rhs); //TODO ptrs???
            }
        }
    }

    //second loop of DFD
    for (auto rhs : possibleRHSs) {
        //тут строим новую решетку, соответственно нужно обнулить/завести некоторые структуры данных
        findLHSs(rhs, schema);
    }
}

shared_ptr<Vertical> DFD::takeRandom(std::list<shared_ptr<Vertical>> const& nodeList) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeList.begin(), nodeList.end()) - 1);
    auto iterator = nodeList.begin();
    std::advance(iterator, dis(this->gen));
    return *iterator;
}

shared_ptr<Vertical> DFD::takeRandom(std::vector<shared_ptr<Vertical>> const& nodeVector) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeVector.begin(), nodeVector.end()) - 1);
    auto iterator = nodeVector.begin();
    std::advance(iterator, dis(this->gen));
    return *iterator;
}

void DFD::findLHSs(shared_ptr<Column> rhs, shared_ptr<RelationalSchema> schema) {

    //std::vector<shared_ptr<Vertical>> minimalDeps; //TODO мб их определять либо в функции execute, либо полями класса
    //std::vector<shared_ptr<Vertical>> maximalNonDeps;

    std::list<shared_ptr<Vertical>> seeds; //TODO лист или вектор?

    //initialize seeds nodes
    for (auto column : schema->getColumns()) {
        //add to seeds all columns except rhs
        if (column->getIndex() != rhs->getIndex())
            seeds.emplace_back(std::make_shared<Vertical>(*column)); //TODO проверить
    }



    while (!seeds.empty()) {
        shared_ptr<Vertical> node = takeRandom(seeds);
        do {
            auto const nodeObservation = observations.find(*node); //const?

            if (nodeObservation != observations.end()) {
                NodeCategory& nodeCategory = nodeObservation->second;

                if (nodeCategory == NodeCategory::candidateMinimalDependency) {
                    nodeCategory = observations.updateDependencyCategory(node);
                    if (nodeCategory == NodeCategory::minimalDependency) {
                        minimalDeps.push_back(node);
                    }
                } else if (nodeCategory == NodeCategory::candidateMaximalNonDependency) {
                    nodeCategory = observations.updateNonDependencyCategory(node);
                    if (nodeCategory == NodeCategory::maximalNonDependency) {
                        maximalNonDeps.push_back(node);
                    }
                }

            } else if (!observations.inferCategory(node)) {
                //auto nodePartition = partitionCache->getOrCreateFor()

            }
            //node = pickNextNode(node);
        } while (node != nullptr);
    }
}

DFD::DFD(const std::filesystem::path &path, char separator, bool hasHeader)
        : FDAlgorithm(path, separator, hasHeader), gen(rd()), observations() {
    relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    dependencies = DependenciesSet(relation->getSchema());
}

shared_ptr<Vertical> DFD::pickNextNode(shared_ptr<Vertical> node) {
    dynamic_bitset<> columnIndices = node->getColumnIndices();
    auto nodeIter = observations.find(*node);

    //можно зарефакторить, если сделать категорию undefined?
    if (nodeIter != observations.end()) {
        if (nodeIter->second == NodeCategory::candidateMinimalDependency) {
            vector<shared_ptr<Vertical>> uncheckedSubsets = dependencies.getUncheckedSubsets(node, observations); //TODO переписать observations в конструктор?
            if (uncheckedSubsets.empty()) {
                minimalDeps.push_back(node);
                dependencies.addNewDependency(node);
            } else {
                //TODO сделать чтобы можно у функции takeRandom на вход подавать итераторы вместо перегрузки
                shared_ptr<Vertical> nextNode = takeRandom(uncheckedSubsets);
                trace.push(node);
                return nextNode;
            }
        } else if (nodeIter->second == NodeCategory::candidateMaximalNonDependency) {
            vector<shared_ptr<Vertical>> uncheckedSupersets = nonDependenciesMap.getUncheckedSupersets(node, observations);
            if (uncheckedSupersets.empty()) {
                maximalNonDeps.push_back(node);
                nonDependenciesMap.addNewNonDependency(node);
            } else {
                shared_ptr<Vertical> nextNode = takeRandom(uncheckedSupersets);
                trace.push(node);
                return nextNode;
            }
        }
    }

    shared_ptr<Vertical> nextNode = trace.empty() ? nullptr : trace.top();
    trace.pop();
    return nextNode;
}