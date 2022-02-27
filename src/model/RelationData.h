//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <vector>

#include "ColumnData.h"
#include "RelationalSchema.h"

template <typename T>
class AbstractRelationData {
public:
    using ColumnType = T;

protected:
    std::unique_ptr<RelationalSchema> schema_;
    std::vector<ColumnType> column_data_;

public:
    static constexpr char const* kNullValue = "NULL";

    explicit AbstractRelationData(std::unique_ptr<RelationalSchema> schema,
                                  std::vector<ColumnType> column_data) noexcept
        : schema_(std::move(schema)), column_data_(std::move(column_data)) {}

    virtual unsigned int GetNumRows() const = 0;

    virtual std::vector<ColumnType>& GetColumnData() { return column_data_; }
    virtual std::vector<ColumnType> const& GetColumnData() const { return column_data_; }
    virtual ColumnType& GetColumnData(unsigned int column_index) {
        return column_data_[column_index];
    }
    virtual ColumnType const& GetColumnData(unsigned int column_index) const {
        return column_data_[column_index];
    }

    virtual ~AbstractRelationData() = default;

    double GetMaximumEntropy() const { return std::log(GetNumRows()); }
    unsigned int GetNumColumns() const { return schema_->GetNumColumns(); }
    double GetMaximumNip() const { return GetNumRows() * (GetNumRows() - 1) / 2.0; }
    unsigned long long GetNumTuplePairs() const {
        return (unsigned long long)GetNumRows() * (GetNumRows() - 1) / 2;
    }
    RelationalSchema const* GetSchema() const { return schema_.get(); }
};

using RelationData = AbstractRelationData<ColumnData>;

