//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "model/Column.h"

using boost::dynamic_bitset, std::string, std::weak_ptr, std::shared_ptr, std::vector;

class Vertical {
private:
    //Vertical(shared_ptr<RelationalSchema>& relSchema, int indices);

    dynamic_bitset<> columnIndices;
    weak_ptr<RelationalSchema> schema;

public:
    //TODO: think about passing dynamic_bitset.
    Vertical(shared_ptr<RelationalSchema> relSchema, dynamic_bitset<> const & indices);
    Vertical() = default;
    explicit Vertical(Column & col);

    Vertical(Vertical& other) = default;
    Vertical(Vertical const& other) = default;
    Vertical& operator=(const Vertical& rhs) = default;
    Vertical(Vertical&& other) noexcept;        //just = default
    Vertical& operator=(Vertical&& rhs) noexcept;
    virtual ~Vertical() = default;
    bool operator==(Vertical const& other) const { return columnIndices == other.columnIndices; }

    dynamic_bitset<> getColumnIndices() const;
    shared_ptr<RelationalSchema> getSchema();
    bool contains(Vertical& that);
    bool intersects(Vertical& that);
    Vertical Union(Vertical const& that);
    Vertical project(Vertical& that);
    Vertical without (Vertical const & that) const;
    Vertical invert();
    Vertical invert(Vertical& scope);
    static Vertical emptyVertical(shared_ptr<RelationalSchema> relSchema);
    int getArity() const;
    vector<shared_ptr<Column>> getColumns() const;
    vector<shared_ptr<Vertical>> getParents();
    //possible to use list or set?

    string toString() const;
    explicit operator std::string() const { return toString(); }
};
