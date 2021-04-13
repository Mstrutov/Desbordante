//
// Created by alexandrsmirn
//

#include "DFD.h"

#include <algorithm>
#include <random>

#include "LatticeNode.h"
#include "../model/ColumnLayoutRelationData.h"
#include "../model/RelationalSchema.h"
#include "../util/PositionListIndex.h"

//using std::shared_ptr;

unsigned long long DFD::execute() {
    shared_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);//second parameter?
    shared_ptr<RelationalSchema> schema = relation->getSchema();

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

shared_ptr<LatticeNode> DFD::takeRandom(std::list<shared_ptr<LatticeNode>> const& nodeList) {
    std::uniform_int_distribution<> dis(0, std::distance(nodeList.begin(), nodeList.end()) - 1);
    auto iterator = nodeList.begin();
    std::advance(iterator, dis(this->gen));
    return *iterator;
}

void DFD::findLHSs(shared_ptr<Column> rhs, shared_ptr<RelationalSchema> schema) {
    /*std::list<shared_ptr<Vertical>> seeds;

    //initialize seeds nodes
    for (auto column : schema->getColumns()) {
        //add to seeds all columns except rhs
        if (column->getIndex() != rhs->getIndex())
            seeds.emplace_back(std::make_shared<Vertical>(*column)); //TO DO проверить
    }*/

    std::vector<shared_ptr<LatticeNode>> minimalDeps; //TODO мб их определять либо в функции execute, либо полями класса
    std::vector<shared_ptr<LatticeNode>> maximalNonDeps;

    std::list<shared_ptr<LatticeNode>> seeds; //TODO лист или вектор?

    for (auto column : schema->getColumns()) {
        //add to seeds all columns except rhs
        if (column->getIndex() != rhs->getIndex()) {
            seeds.emplace_back(std::make_shared<LatticeNode>(*column)); //TODO проверить
        }
    }

    while (!seeds.empty()) {
        shared_ptr<LatticeNode> node = takeRandom(seeds);
        do {
            if (node->isVisited() && node->isCandidate()) {
                if (node->isMinimalDependency()) {
                    minimalDeps.push_back(node);
                } else if (node->isMaximalNonDependency()) {
                    maximalNonDeps.push_back(node);
                }
                node.updateDependencyType();
            } else {
                if (!inferCategory(node)) {
                    computePartitions(node);
                }
            }

            node = pickNextNode();
        } while (node != nullptr);
    }
}

shared_ptr<LatticeNode> pickNextNode(shared_ptr<LatticeNode> node) {

}