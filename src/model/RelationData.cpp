//
// Created by kek on 26.07.2019.
//

#include "RelationData.h"
#include "../util/PositionListIndex.h"

const int RelationData::nullValueId = -1;
const int RelationData::singletonValueId = PositionListIndex::singletonValueId;

RelationData::RelationData(shared_ptr<RelationalSchema>& schema): schema(schema) {}

int RelationData::getNumColumns() {
    return schema->getNumColumns();
}

shared_ptr<RelationalSchema> RelationData::getSchema() {
    return schema;
}

int RelationData::getNumTuplePairs() {
    return this->getNumRows() * (this->getNumRows() - 1) / 2;
}
