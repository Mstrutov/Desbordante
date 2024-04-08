//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <functional>
#include <list>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>

#include "bitset_utils.h"

class Column;

class Vertical;

class RelationalSchema {
private:
    std::vector<std::unique_ptr<Column>> columns_;
    std::string name_;

public:
    std::unique_ptr<Vertical> empty_vertical_;

    RelationalSchema(std::string name);

    RelationalSchema(RelationalSchema const& other) = delete;
    RelationalSchema& operator=(RelationalSchema const& rhs) = delete;
    RelationalSchema(RelationalSchema&& other) noexcept = default;
    RelationalSchema& operator=(RelationalSchema&& rhs) noexcept = default;

    void Init();

    std::string GetName() const {
        return name_;
    }

    std::vector<std::unique_ptr<Column>> const& GetColumns() const {
        return columns_;
    };

    bool IsColumnInSchema(std::string const& col_name) const;

    Column const* GetColumn(std::string const& col_name) const;
    Column const* GetColumn(size_t index) const;
    size_t GetNumColumns() const;
    Vertical GetVertical(boost::dynamic_bitset<> indices) const;

    void AppendColumn(std::string const& col_name);
    void AppendColumn(Column column);

    template <typename Container>
    boost::dynamic_bitset<> IndicesToBitset(Container const& indices) const;
    template <typename ForwardIt>
    boost::dynamic_bitset<> IndicesToBitset(ForwardIt begin, ForwardIt end) const;

    std::unordered_set<Vertical> CalculateHittingSet(
            std::vector<Vertical> verticals,
            boost::optional<std::function<bool(Vertical const&)>> pruning_function) const;

    ~RelationalSchema();

    friend inline bool operator==(RelationalSchema const& l, RelationalSchema const& r);
};

inline bool operator==(RelationalSchema const& l, RelationalSchema const& r) {
    return (l.name_ == r.name_ && l.columns_.size() == r.columns_.size());
}

inline bool operator!=(RelationalSchema const& l, RelationalSchema const& r) {
    return !(l == r);
}

template <typename Container>
boost::dynamic_bitset<> RelationalSchema::IndicesToBitset(Container const& indices) const {
    return util::IndicesToBitset(indices, GetNumColumns());
}

template <typename ForwardIt>
boost::dynamic_bitset<> RelationalSchema::IndicesToBitset(ForwardIt begin, ForwardIt end) const {
    return util::IndicesToBitset(begin, end, GetNumColumns());
}
