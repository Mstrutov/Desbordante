#pragma once

#include <set>

class Itemset {
private:
    std::set<unsigned> indices;
public:
    std::set<unsigned> const& getItemsIDs() const noexcept { return indices; }
    void addItemID(unsigned itemID) { indices.insert(itemID); }

    //TODO конструктор от одного айдишника?
};
