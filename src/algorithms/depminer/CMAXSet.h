#pragma once

#include <unordered_set>

#include "ColumnLayoutRelationData.h"
#include "custom/CustomHashes.h"

class CMAXSet {
private:
    Column column_;
    std::unordered_set<Vertical> column_combinations_;
public:
    explicit CMAXSet(Column const& column) : column_(column) {};
    void MakeNewCombinations(std::unordered_set<Vertical> comb) {
        this->column_combinations_ = std::move(comb);
    }
    void AddCombination(Vertical const& combination) { column_combinations_.insert(combination); }
    std::unordered_set<Vertical> const& GetCombinations() const { return column_combinations_; }
    Column const& GetColumn() const { return column_; }
};
