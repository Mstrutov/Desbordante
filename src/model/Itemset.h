#pragma once

#include <vector>

class Itemset {
private:
    std::vector<unsigned> indices;
public:
    std::vector<unsigned> const& getItemsIDs() const noexcept { return indices; }
    void addItemID(unsigned itemID) { indices.push_back(itemID); }
    void sort();
    //TODO конструктор от одного айдишника?
};
