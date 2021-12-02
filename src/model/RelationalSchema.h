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
    std::vector<std::unique_ptr<Column>> columns;
    std::string name;
    bool isNullEqNull;

public:
    std::unique_ptr<Vertical> emptyVertical;

    RelationalSchema(std::string name, bool isNullEqNull);
    void init();
    std::string getName() const { return name; }
    std::vector<std::unique_ptr<Column>> const& getColumns() const { return columns; };
    Column const* getColumn(const std::string &colName) const;
    Column const* getColumn(int index) const;
    size_t getNumColumns() const;
    Vertical getVertical(boost::dynamic_bitset<> indices) const;
    bool isNullEqualNull() const;

    void appendColumn(const std::string& colName);
    void appendColumn(Column column);

    std::unordered_set<Vertical> calculateHittingSet(std::vector<Vertical> verticals,
            boost::optional<std::function<bool (Vertical const&)>> pruningFunction) const;

    ~RelationalSchema();
};
