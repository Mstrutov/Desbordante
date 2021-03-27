//
// Created by alexandrsmirn
//

#include "DFD.h"

#include "../model/ColumnLayoutRelationData.h"
#include "../model/RelationalSchema.h"

unsigned long long DFD::execute() {
    shared_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);//second parameter?
    shared_ptr<RelationalSchema> schema = relation->getSchema();

    for (auto column : schema->getColumns()) {
        shared_ptr<ColumnData> columnData = relation->getColumnData(column->getIndex());
    }

}
