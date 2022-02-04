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
    static constexpr int kNullValueId = -1;

    virtual unsigned int GetNumRows() const = 0;
    virtual std::vector<ColumnData>& GetColumnData() = 0;
    virtual std::vector<ColumnData> const& GetColumnData() const = 0;
    virtual ColumnData& GetColumnData(int column_index) = 0;
    virtual ColumnData const& GetColumnData(int column_index) const = 0;
    virtual std::vector<int> GetTuple(int tuple_index) const = 0;

    virtual ~RelationData();

    //virtual void shuffleColumns()=0;

    unsigned int GetNumColumns() const { return schema_->GetNumColumns(); }
    double GetMaximumNip() const { return GetNumRows() * (GetNumRows() - 1) / 2.0; }
    unsigned long long GetNumTuplePairs() const { return (unsigned long long) GetNumRows() * (GetNumRows() - 1) / 2; }
    RelationalSchema const* GetSchema() const { return schema_.get(); }

protected:
    std::unique_ptr<RelationalSchema> schema_;

    explicit RelationData(std::unique_ptr<RelationalSchema> schema) : schema_(std::move(schema)) {}
};
