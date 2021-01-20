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
    static const int singletonValueId;

    //c++17: inline initialization || constexpr
    static const int nullValueId;

    virtual unsigned int getNumRows()=0;
    unsigned int getNumColumns();
    virtual vector<shared_ptr<ColumnData>> getColumnData()=0;
    virtual shared_ptr<ColumnData> getColumnData(int columnIndex)=0;
    virtual vector<int> getTuple(int tupleIndex)=0;
    virtual void shuffleColumns()=0;
    double getMaximumNip() { return getNumRows() * (getNumRows() - 1) / 2.0; }
    unsigned long long getNumTuplePairs() { return (unsigned long long) getNumRows() * (getNumRows() - 1) / 2; }

    shared_ptr<RelationalSchema> getSchema();
protected:
    explicit RelationData(shared_ptr<RelationalSchema>& schema);
    shared_ptr<RelationalSchema> schema;
};
