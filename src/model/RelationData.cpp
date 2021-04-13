//
// Created by kek on 26.07.2019.
//

#include "RelationData.h"
#include "../util/PositionListIndex.h"

const int RelationData::nullValueId = -1;
const int RelationData::singletonValueId = PositionListIndex::singletonValueId;

RelationData::RelationData(shared_ptr<RelationalSchema> const& schema): schema(schema) {}

unsigned int RelationData::getNumColumns() const {
    return schema->getNumColumns();
}

shared_ptr<RelationalSchema> RelationData::getSchema() const {
    return schema;
}

int RelationData::getNumTuplePairs() {
    return this->getNumRows() * (this->getNumRows() - 1) / 2;
}
