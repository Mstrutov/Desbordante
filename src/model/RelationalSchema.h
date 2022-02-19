//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>

#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>

class Column;

class Vertical;

class RelationalSchema {
private:
    std::vector<std::unique_ptr<Column>> columns_;
    std::string name_;
    bool is_null_eq_null_;

public:
    std::unique_ptr<Vertical> empty_vertical_;

    RelationalSchema(std::string name, bool is_null_eq_null);
    void Init();
    std::string GetName() const { return name_; }
    std::vector<std::unique_ptr<Column>> const& GetColumns() const { return columns_; };
    Column const* GetColumn(const std::string& col_name) const;
    Column const* GetColumn(int index) const;
    size_t GetNumColumns() const;
    Vertical GetVertical(boost::dynamic_bitset<> indices) const;
    bool IsNullEqualNull() const;

    void AppendColumn(const std::string& col_name);
    void AppendColumn(Column column);

    std::unordered_set<Vertical> CalculateHittingSet(
        std::vector<Vertical> verticals,
        boost::optional<std::function<bool(Vertical const&)>> pruning_function) const;

    ~RelationalSchema();

    friend inline bool operator==(RelationalSchema const& l, RelationalSchema const& r);
};

inline bool operator==(RelationalSchema const& l, RelationalSchema const& r) {
    return (l.name_ == r.name_ && l.is_null_eq_null_ == r.is_null_eq_null_ &&
            l.columns_.size() == r.columns_.size());
}
inline bool operator!=(RelationalSchema const& l, RelationalSchema const& r) {
    return !(l == r);
}

