//
// Created by kek on 26.07.2019.
//

#include "RelationData.h"

const int RelationData::nullValueId = -1;

RelationData::RelationData(shared_ptr<RelationalSchema> const& schema): schema(schema) {}

unsigned int RelationData::getNumColumns() const {
    return schema->getNumColumns();
}

shared_ptr<RelationalSchema> RelationData::getSchema() const {
    return schema;
}