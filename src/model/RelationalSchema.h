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

using boost::dynamic_bitset, std::string, std::vector;

class RelationalSchema : public std::enable_shared_from_this<RelationalSchema> {
private:
    RelationalSchema(string name, bool isNullEqNull);
    vector<std::shared_ptr<Column>> columns;
    string name;
    bool isNullEqNull;

public:
    std::shared_ptr<Vertical> emptyVertical;

    void init();
    static std::shared_ptr<RelationalSchema> create(string name, bool isNullEqNull);
    string getName();
    vector<std::shared_ptr<Column>> getColumns();
    std::shared_ptr<Column> getColumn(const string &colName);
    std::shared_ptr<Column> getColumn(int index);
    void appendColumn(const string& colName);
    void appendColumn(std::shared_ptr<Column> column);
    int getNumColumns();
    std::shared_ptr<Vertical> getVertical(dynamic_bitset<> indices);
    bool isNullEqualNull();

    std::unordered_set<std::shared_ptr<Vertical>> calculateHittingSet(std::list<std::shared_ptr<Vertical>>&& verticals,
            boost::optional<std::function<bool (Vertical const&)>> pruningFunction);
};
