#pragma once

#include <vector>

class Itemset {
private:
    std::vector<unsigned> indices_;
public:
    std::vector<unsigned> const& GetItemsIDs() const noexcept { return indices_; }
    void AddItemId(unsigned item_id) { indices_.push_back(item_id); }
    void Sort();
    //TODO конструктор от одного айдишника?
};
