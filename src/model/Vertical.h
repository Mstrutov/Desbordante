//
// Created by kek on 17.07.19.
//

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

using boost::dynamic_bitset, std::string, std::weak_ptr, std::shared_ptr, std::vector;

class RelationalSchema;

class Vertical {
protected:
    //Vertical(shared_ptr<RelationalSchema>& relSchema, int indices);

    dynamic_bitset<> columnIndices;
    weak_ptr<RelationalSchema> schema;

public:
    Vertical(shared_ptr<RelationalSchema>& relSchema, int indices);
    Vertical() = default;
    Vertical(Vertical& other) = default;
    Vertical& operator=(const Vertical& rhs) = default;
    Vertical(Vertical&& other) noexcept;        //just = default
    Vertical& operator=(Vertical&& rhs) noexcept ;

    dynamic_bitset<> getColumnIndices();    //return const &
    shared_ptr<RelationalSchema> getSchema();
    bool contains(Vertical& that);
    bool intersects(Vertical& that);
    Vertical Union(Vertical& that);
    Vertical project(Vertical& that);
    Vertical without (Vertical& that);
    Vertical invert();
    Vertical invert(Vertical& scope);
    static Vertical emptyVertical(shared_ptr<RelationalSchema> relSchema);
    int getArity();
    vector<shared_ptr<Vertical>> getColumns();          //Very strange - orig demands vector<column> - how???
    //possible to use list or set?

    virtual string toString();
};
