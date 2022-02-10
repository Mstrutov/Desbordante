#pragma once

#include <vector>
#include <list>

struct Node {
    std::vector<unsigned> items;
    std::list<Node> children;

    Node() = default;
    Node(Node const& node) = delete;

    explicit Node(unsigned itemID)
            : items({itemID}) {}

    explicit Node(std::vector<unsigned> && itemsToAdd)
            : Node() {
        std::swap(itemsToAdd, items);
    }
};