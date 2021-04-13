//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <iostream>
#include <vector>

#include "ColumnData.h"
#include "RelationalSchema.h"

using std::vector;

class RelationData {
public:
    //static int singletonValueId;

    //c++17: inline initialization || constexpr
    static const int nullValueId;

    virtual unsigned int getNumRows() const = 0;
    unsigned int getNumColumns() const;
    virtual vector<shared_ptr<ColumnData>> getColumnData() const = 0;
    virtual shared_ptr<ColumnData> getColumnData(int columnIndex) const = 0;
    virtual vector<int> getTuple(int tupleIndex) const = 0;
    virtual void shuffleColumns()=0;
    double getMaximumNip() const { return getNumRows() * (getNumRows() - 1) / 2.0; }
    unsigned long long getNumTuplePairs() const { return (unsigned long long) getNumRows() * (getNumRows() - 1) / 2; }

    shared_ptr<RelationalSchema> getSchema() const;
protected:
    explicit RelationData(shared_ptr<RelationalSchema> const& schema);
    shared_ptr<RelationalSchema> schema;
};
