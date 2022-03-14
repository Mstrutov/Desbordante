#pragma once

#include <list>
#include <vector>

namespace algos {

struct Node {
    std::vector<unsigned> items;
    double support = 0;
    std::vector<Node> children;

    Node() = default;
    Node(Node&& other) = default;
    Node& operator=(Node&& other) = default;

    Node(Node const& node) = delete;
    Node& operator=(Node const&) = delete;

    explicit Node(unsigned item_id)
        : items({item_id}), support(0) {}

    explicit Node(std::vector<unsigned>&& items_to_add)
        : Node() {
        std::swap(items_to_add, items);
    }
};

} // namespace algos
