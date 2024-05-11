#include "search_tree.h"

namespace algos {

SearchTreeEulerFD::Node::Node(size_t bit, SearchTreeEulerFD::Bitset set, SearchTreeEulerFD::Bitset sets_union,
                       SearchTreeEulerFD::Bitset sets_inter, std::shared_ptr<Node> const& parent,
                       std::shared_ptr<Node> left, std::shared_ptr<Node> right)
    : bit_(bit),
      set_(std::move(set)),
      union_(std::move(sets_union)),
      inter_(std::move(sets_inter)),
      left_(std::move(left)),
      right_(std::move(right)),
      parent_(parent) {}

SearchTreeEulerFD::Node::Node(size_t bit, SearchTreeEulerFD::Bitset sets_union, SearchTreeEulerFD::Bitset sets_inter,
                       std::shared_ptr<Node> const& parent, std::shared_ptr<Node> left,
                       std::shared_ptr<Node> right)
    : bit_(bit),
      union_(std::move(sets_union)),
      inter_(std::move(sets_inter)),
      left_(std::move(left)),
      right_(std::move(right)),
      parent_(parent) {}

SearchTreeEulerFD::Node::Node(size_t bit, SearchTreeEulerFD::Bitset set, std::shared_ptr<Node> const& parent)
    : bit_(bit), set_(std::move(set)), parent_(parent) {}

bool SearchTreeEulerFD::Node::IsLeaf() const {
    assert((left_ != nullptr) ^ (right_ == nullptr));
    return !left_ && !right_;
}

SearchTreeEulerFD::Bitset const& SearchTreeEulerFD::Node::GetUnion() const {
    return IsLeaf() ? set_ : union_;
}

SearchTreeEulerFD::Bitset const& SearchTreeEulerFD::Node::GetInter() const {
    return IsLeaf() ? set_ : inter_;
}

SearchTreeEulerFD::SearchTreeEulerFD(size_t number_of_attributes) : number_of_attributes_(number_of_attributes) {}

SearchTreeEulerFD::SearchTreeEulerFD(Bitset const& set) : number_of_attributes_(set.size()) {
    CreateSingleElementSets(set);
}

void SearchTreeEulerFD::CreateSingleElementSets(Bitset const& set) {
    for (size_t bit = set.find_first(); bit != Bitset::npos; bit = set.find_next(bit)) {
        Bitset bitset_to_add(number_of_attributes_);
        bitset_to_add.set(bit);
        Add(bitset_to_add);
    }
}

bool SearchTreeEulerFD::Add(Bitset const& set) {
    assert(!set.empty());
    assert(set.size() == number_of_attributes_);

    if (!root_) {
        root_ = std::make_shared<Node>(set.find_first(), set, nullptr);
        ++cardinality_;
        return true;
    }

    auto current_node = root_;
    for (size_t prev_node_bit = 0, curr_node_bit = current_node->bit_; !current_node->IsLeaf();
         current_node = set[curr_node_bit] ? current_node->right_ : current_node->left_,
                prev_node_bit = std::exchange(curr_node_bit, current_node->bit_) + 1) {
        assert(curr_node_bit != Bitset::npos);
        assert(prev_node_bit <= curr_node_bit);
        for (size_t set_bit = prev_node_bit; set_bit < curr_node_bit; ++set_bit) {
            if (set[set_bit] != current_node->union_[set_bit]) {
                InsertLeafIntoMiddle(current_node, set, set_bit);
                ++cardinality_;
                return true;
            }
        }
    }

    auto bits = FindNodeAndSetBits(current_node->set_, set);
    size_t node_bit = bits.first;
    size_t set_bit = bits.second;
    if (node_bit == set_bit) {
        return false;
    }

    InsertLeafIntoEnd(current_node, set, node_bit, set_bit);
    ++cardinality_;
    return true;
}

bool SearchTreeEulerFD::Remove(Bitset const& set) {
    assert(!set.empty());
    if (!root_) {
        return false;
    }

    auto node_to_remove = FindNode(set);
    if (!node_to_remove) {
        return false;
    }

    CutLeaf(node_to_remove);
    --cardinality_;
    return true;
}

void SearchTreeEulerFD::CollectSubsets(Bitset const& set, std::shared_ptr<Node> const& current_node,
                                BitsetConsumer const& collect, bool& go_further) const {
    if (current_node->IsLeaf()) {
        auto node_set = current_node->set_;
        if (node_set.is_subset_of(set)) {
            collect(node_set);
        }
        return;
    }

    if (go_further && current_node->inter_.is_subset_of(set)) {
        CollectSubsets(set, current_node->left_, collect, go_further);
        if (set[current_node->bit_]) {
            CollectSubsets(set, current_node->right_, collect, go_further);
        }
    }
}

bool SearchTreeEulerFD::ContainsAnySubsetOf(Bitset const& set) const {
    bool go_further = true;
    auto contains_func = [&go_further](Bitset const&) { go_further = false; };

    if (root_) {
        CollectSubsets(set, root_, contains_func, go_further);
    }

    return !go_further;
}

void SearchTreeEulerFD::ForEach(std::shared_ptr<Node> const& current_node,
                         BitsetConsumer const& collect) const {
    if (current_node->IsLeaf()) {
        collect(current_node->set_);
        return;
    }

    ForEach(current_node->left_, collect);
    ForEach(current_node->right_, collect);
}

void SearchTreeEulerFD::ForEach(BitsetConsumer const& collect) const {
    if (root_) {
        ForEach(root_, collect);
    }
}

void SearchTreeEulerFD::ForEachSubset(SearchTreeEulerFD::Bitset const& set, BitsetConsumer const& collect) const {
    if (root_) {
        bool go_further = true;
        CollectSubsets(set, root_, collect, go_further);
    }
}

void SearchTreeEulerFD::UpdateInterAndUnion(std::shared_ptr<Node> const& node) {
    auto node_copy = node;
    while (node_copy) {
        auto const& left = node_copy->left_;
        auto const& right = node_copy->right_;
        node_copy->union_ = left->GetUnion() | right->GetUnion();
        node_copy->inter_ = left->GetInter() & right->GetInter();
        node_copy = node_copy->parent_.lock();
    }
}

std::shared_ptr<SearchTreeEulerFD::Node> SearchTreeEulerFD::FindNode(Bitset const& set) {
    auto current_node = root_;
    for (size_t set_bit = set.find_first(), node_bit = current_node->bit_; !current_node->IsLeaf();
         node_bit = current_node->bit_) {
        if (!(set.is_subset_of(current_node->union_) && current_node->inter_.is_subset_of(set))) {
            return nullptr;
        }

        if (set[node_bit]) {
            current_node = current_node->right_;
            set_bit = set.find_next(set_bit);
        } else {
            current_node = current_node->left_;
        }
    }

    if (current_node->set_ != set) {
        return nullptr;
    }

    return current_node;
}

void SearchTreeEulerFD::CutLeaf(std::shared_ptr<Node> const& node_to_remove) {
    if (node_to_remove == root_) {
        root_ = nullptr;
        return;
    }

    auto parent_node = node_to_remove->parent_.lock();
    auto another_child_node =
            (parent_node->right_ == node_to_remove) ? parent_node->left_ : parent_node->right_;

    parent_node->left_ = another_child_node->left_;
    parent_node->right_ = another_child_node->right_;
    parent_node->set_ = another_child_node->set_;
    parent_node->union_ = another_child_node->union_;
    parent_node->inter_ = another_child_node->inter_;
    if (!another_child_node->IsLeaf()) {
        parent_node->bit_ = another_child_node->bit_;
        parent_node->left_->parent_ = parent_node;
        parent_node->right_->parent_ = parent_node;
    }

    UpdateInterAndUnion(parent_node->parent_.lock());
}

std::pair<size_t, size_t> SearchTreeEulerFD::FindNodeAndSetBits(Bitset const& node_set,
                                                         Bitset const& set) {
    size_t node_bit = node_set.find_first();
    size_t set_bit = set.find_first();
    for (; node_bit == set_bit;
         node_bit = node_set.find_next(node_bit), set_bit = set.find_next(set_bit)) {
        if (node_bit == Bitset::npos) {
            break;
        }
    }

    return std::make_pair(node_bit, set_bit);
}

void SearchTreeEulerFD::InsertLeafIntoEnd(std::shared_ptr<Node> const& current_node, Bitset const& set,
                                   size_t node_bit, size_t set_bit) {
    Bitset const& node_set = current_node->set_;
    auto new_left = std::make_shared<Node>(node_bit, node_set, current_node);
    auto new_right = std::make_shared<Node>(set_bit, set, current_node);

    if (node_bit < set_bit) {
        new_left->bit_ = node_set.find_next(node_bit);
        std::swap(new_left, new_right);
        current_node->bit_ = node_bit;
    } else {
        new_right->bit_ = set.find_next(set_bit);
        current_node->bit_ = set_bit;
    }

    current_node->left_ = new_left;
    current_node->right_ = new_right;

    UpdateInterAndUnion(current_node);
}

void SearchTreeEulerFD::InsertLeafIntoMiddle(std::shared_ptr<Node> const& current_node,
                                      SearchTreeEulerFD::Bitset const& set, size_t set_bit) {
    auto new_left = std::make_shared<Node>(current_node->bit_, current_node->set_,
                                           current_node->union_, current_node->inter_, current_node,
                                           current_node->left_, current_node->right_);
    auto new_right = std::make_shared<Node>(set.find_next(set_bit), set, current_node);
    auto replacing_node = new_left;
    if (!set[set_bit]) {
        std::swap(new_left, new_right);
        replacing_node = new_right;
    }

    current_node->left_->parent_ = replacing_node;
    current_node->right_->parent_ = replacing_node;

    current_node->bit_ = set_bit;
    current_node->left_ = new_left;
    current_node->right_ = new_right;

    UpdateInterAndUnion(current_node);
}

bool SearchTreeEulerFD::SupersetsTraverse(
  Bitset const& set, std::shared_ptr<Node> const& current_node) const {
    if (current_node->IsLeaf()) {
        return set.is_subset_of(current_node->set_);
    }
    if (!set.is_subset_of(current_node->union_)) {
        return false;
    }
    return SupersetsTraverse(set, current_node->left_) ||
           SupersetsTraverse(set, current_node->right_);
}

bool SearchTreeEulerFD::ContainsAnySupersetOf(Bitset const& set) const {
    return root_ ? SupersetsTraverse(set, root_) : false;
}

} // namespace algos
