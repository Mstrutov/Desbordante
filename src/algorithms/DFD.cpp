//
// Created by alexandrsmirn
//

#include "DFD.h"

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
            auto lhs = Vertical(**columnIter);

            for (auto rhs : possibleRHSs) {
                this->registerFD(lhs, *rhs); //TODO ptrs???
            }
        }
    }

}
