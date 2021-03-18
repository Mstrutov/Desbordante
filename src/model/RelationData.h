//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <iostream>
#include <vector>

#include "ColumnData.h"
#include "RelationalSchema.h"

class RelationData {
public:
    static constexpr int nullValueId = -1;

    virtual unsigned int getNumRows() const = 0;
    virtual std::vector<ColumnData>& getColumnData() = 0;
    virtual std::vector<ColumnData> const& getColumnData() const = 0;
    virtual ColumnData& getColumnData(int columnIndex) = 0;
    virtual ColumnData const& getColumnData(int columnIndex) const = 0;
    virtual std::vector<int> getTuple(int tupleIndex) const = 0;

    //virtual void shuffleColumns()=0;

    unsigned int getNumColumns() const { return schema->getNumColumns(); }
    double getMaximumNip() const { return getNumRows() * (getNumRows() - 1) / 2.0; }
    unsigned long long getNumTuplePairs() const { return (unsigned long long) getNumRows() * (getNumRows() - 1) / 2; }
    RelationalSchema const* getSchema() const { return schema.get(); }

protected:
    std::unique_ptr<RelationalSchema> schema;

    explicit RelationData(std::unique_ptr<RelationalSchema> schema) : schema(std::move(schema)) {}
};
