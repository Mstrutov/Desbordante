//
// Created by kek on 26.07.2019.
//

#pragma once

#include <iostream>
#include <vector>

#include "model/ColumnData.h"
#include "model/RelationalSchema.h"

using std::vector;

class RelationData {
public:
    //static int singletonValueId;

    //c++17: inline initialization || constexpr
    static const int nullValueId;

    virtual unsigned int getNumRows()=0;
    int getNumColumns();
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
