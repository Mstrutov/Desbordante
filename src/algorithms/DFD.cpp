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

    //first loop of DFD
    for (auto column : schema->getColumns()) {
        shared_ptr<ColumnData> columnData = relation->getColumnData(column->getIndex());
        shared_ptr<PositionListIndex> columnPLI = columnData->getPositionListIndex();

        //if current column is unique
        if (columnPLI->getNumNonSingletonCluster() == 0) {
            auto lhs = Vertical(*column);

            for (auto rhs : schema->getColumns()) {
                if (rhs->getIndex() != column->getIndex()) {
                    this->registerFD(lhs, *rhs); //TODO ptrs???
                }
            }
        }
    }

}
