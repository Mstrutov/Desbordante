#pragma once

#include <vector>
#include <list>

struct Node {
    std::vector<unsigned> items;
    double support = 0;
    std::vector<Node> children;

    Node() = default;
    Node(Node&& other) = default;
    Node(Node const& node) = delete;
    Node& operator=(Node const&) = delete;

    explicit Node(unsigned itemID)
            : items({itemID}), support(0) {}

    explicit Node(std::vector<unsigned> && itemsToAdd)
            : Node() {
        std::swap(itemsToAdd, items);
    }
};