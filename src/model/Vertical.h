//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Column.h"

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

    //Vertical(Vertical& other) = default;
    Vertical(Vertical const& other) = default;
    Vertical& operator=(const Vertical& rhs) = default;
    Vertical(Vertical&& other) = default;        //just = default
    Vertical& operator=(Vertical&& rhs) = default;
    virtual ~Vertical() = default;
    bool operator==(Vertical const& other) const { return columnIndices == other.columnIndices; }

    dynamic_bitset<> getColumnIndices() const { return columnIndices; }
    shared_ptr<RelationalSchema> getSchema() { return schema.lock(); }
    bool contains(Vertical& that);
    bool intersects(Vertical& that);
    std::shared_ptr<Vertical> Union(Vertical const& that);
    std::shared_ptr<Vertical> project(Vertical& that);
    std::shared_ptr<Vertical> without (Vertical const & that) const;
    std::shared_ptr<Vertical> invert();
    std::shared_ptr<Vertical> invert(Vertical& scope);
    static Vertical emptyVertical(shared_ptr<RelationalSchema> relSchema);
    unsigned int getArity() const { return columnIndices.count(); }
    vector<shared_ptr<Column>> getColumns() const;
    vector<shared_ptr<Vertical>> getParents();
    //possible to use list or set?

    string toString() const;
    explicit operator std::string() const { return toString(); }
};
