#pragma once

// see ../algorithms/cfd/LICENSE

#include <map>
#include <unordered_set>
#include <queue>
#include <iostream>
#include <algorithm>

template<typename T>
size_t HashCollection(const T& xs) {
	size_t res = 0;
	for (const typename T::value_type& x : xs) {
		res ^= std::hash<typename T::value_type>()(x) + 0x9e3779b9 + (res << 6) + (res >> 2);
	}
	return res;
}

template <typename Key, typename Value>
class PrefixTree {
public:
	PrefixTree();
	void Reserve(int);
	void Insert(const Key&, const Value&);
	Value* Find(const Key&) const;
	void Erase(const Key&) const;
	bool HasSubset(const Key&, const Value&) const;
	bool HasStrictSubset(const Key&, const Value&) const;
	std::vector<Key> GetSubsets(const Key&, const Value&) const;
	std::vector<const Value*> GetSubsets(const Key&) const;
	std::vector<Key> GetSets() const;
    [[maybe_unused]] std::map<Key, int> GetSupports(const Key&) const;
private:
	struct PrefixNode {
		Value value;
		std::map<typename Key::value_type, PrefixNode> sub_trees;
		PrefixNode* parent;
		typename Key::value_type key;
		int depth;
	};
private:
	int size_;
	PrefixNode root_;
	std::unordered_set<size_t> hashes_;
	std::unordered_map<size_t, std::vector<PrefixNode*> > jumps_;
};

template <typename Key, typename Value>
PrefixTree<Key,Value>::PrefixTree() {
    root_.depth = 0;
    root_.value = Value();
}

template <typename Key, typename Value>
void PrefixTree<Key,Value>::Reserve(int r) {
    size_ = r;
	//root_.sub_trees.reserve(r);
}

template <typename Key, typename Value>
void PrefixTree<Key,Value>::Insert(const Key& k, const Value& v) {
	typename Key::const_iterator it = k.begin();
	PrefixNode* insertion_point = &root_;
	while (it != k.end()) {
		PrefixNode* next = &insertion_point->sub_trees[*it];
		if (!next->parent) {
			next->key = *it;
			next->depth = insertion_point->depth + 1;
			next->value = Value();
			next->parent = insertion_point;
		}
        insertion_point = next;
		it++;
	}
    jumps_[HashCollection(k)].push_back(insertion_point);
    insertion_point->value = v;
}

template <typename Key, typename Value>
Value* PrefixTree<Key,Value>::Find(const Key& k) const {
	const auto& elem = jumps_.find(HashCollection(k));
	if (elem == jumps_.end()) return 0;

	for (PrefixNode* prefix_node : elem->second) {
		if ((int)k.size() != prefix_node->depth) continue;
		PrefixNode* src = prefix_node;
		typename Key::const_reverse_iterator rit = k.rbegin();
		while (rit != k.rend() && prefix_node && prefix_node->key == *rit) {
            prefix_node = prefix_node->parent;
			rit++;
		}
		if (rit == k.rend() && prefix_node == &root_) {
			if (src->value != Value())
				return &src->value;
			else
				return 0;
		}
	}
	return 0;
}

template <typename Key, typename Value>
void PrefixTree<Key,Value>::Erase(const Key& k) const {
	const auto& elem = jumps_.find(HashCollection(k));
	if (elem == jumps_.end()) return;

	for (PrefixNode* prefix_node : elem->second) {
		if (k.size() != prefix_node->depth) continue;
		PrefixNode* src = prefix_node;
		typename Key::const_reverse_iterator rit = k.rbegin();
		while (rit != k.rend() && prefix_node && prefix_node->key == *rit) {
            prefix_node = prefix_node->parent;
			rit++;
		}
		if (rit == k.rend() && prefix_node == &root_) {
			src->value = Value();
			return;
		}
	}
}

template <typename Key, typename Value>
bool PrefixTree<Key,Value>::HasSubset(const Key& k, const Value& v) const {
	std::queue<const PrefixNode*> fringe;
	fringe.push(&root_);
	while (!fringe.empty()) {
		const PrefixNode* elem = fringe.front();
		fringe.pop();
		if (elem->value == v) return true;
		for (const auto& sub : elem->sub_trees) {
			if (std::binary_search(k.begin(), k.end(), sub.first)) {
				fringe.push(&sub.second);
			}
		}
	}
	return false;
}

template <typename Key, typename Value>
bool PrefixTree<Key,Value>::HasStrictSubset(const Key& k, const Value& v) const {
	std::queue<const PrefixNode*> fringe;
	fringe.push(&root_);
	int d = 0;
	while (!fringe.empty() && d < k.size()) {
		const PrefixNode* elem = fringe.front();
		fringe.pop();
		if (elem->value == v) return true;
		for (const auto& sub : elem->sub_trees) {
			if (std::binary_search(k.begin(), k.end(), sub.first)) {
				fringe.push(&sub.second);
			}
		}
		d++;
	}
	return false;
}

template <typename Key, typename Value>
std::vector<Key> PrefixTree<Key,Value>::GetSubsets(const Key& k, const Value& v) const {
    std::vector<Key> subs;
    std::queue<std::pair<const PrefixNode*, Key> > fringe;
    fringe.push(std::make_pair(&root_, Key()));
    while (!fringe.empty()) {
        const PrefixNode* elem = fringe.front().first;
        Key elem_key = fringe.front().second;
        fringe.pop();
        if (elem->value == v) {
            subs.push_back(elem_key);
        }
        for (const auto& sub : elem->sub_trees) {
            if (std::binary_search(k.begin(), k.end(), sub.first)) {
                Key sub_key(elem_key);
                sub_key.push_back(sub.first);
                fringe.push(std::make_pair(&sub.second, sub_key));
            }
        }
    }
    return subs;
}


template <typename Key, typename Value>
std::vector<const Value*> PrefixTree<Key,Value>::GetSubsets(const Key& k) const {
    std::vector<const Value*> subs;
    std::queue<std::pair<const PrefixNode*, Key> > fringe;
    fringe.push(std::make_pair(&root_, Key()));
    while (!fringe.empty()) {
        const PrefixNode* elem = fringe.front().first;
        Key elem_key = fringe.front().second;
        fringe.pop();
        if (elem->value != Value()) {
            subs.push_back(&elem->value);
        }
        for (const auto& sub : elem->sub_trees) {
            if (std::binary_search(k.begin(), k.end(), sub.first)) {
                Key sub_key(elem_key);
                sub_key.push_back(sub.first);
                fringe.push(std::make_pair(&sub.second, sub_key));
            }
        }
    }
    return subs;
}

template <typename Key, typename Value>
std::vector<Key> PrefixTree<Key,Value>::GetSets() const {
	std::vector<Key> subs;
	std::queue<std::pair<const PrefixNode*, Key> > fringe;
	fringe.push(std::make_pair(&root_, Key()));
	while (!fringe.empty()) {
		const PrefixNode* elem = fringe.front().first;
		Key elem_key = fringe.front().second;
		fringe.pop();
		if (elem->value != Value()) {
			subs.push_back(elem_key);
		}
		for (const auto& sub : elem->sub_trees) {
			Key sub_key(elem_key);
            sub_key.push_back(sub.first);
			fringe.push(std::make_pair(&sub.second, sub_key));
		}
	}
	return subs;
}

template <typename Key, typename Value>
[[maybe_unused]] std::map<Key, int> PrefixTree<Key,Value>::GetSupports(const Key& k) const {
	std::map<Key, int> result;
	std::queue<std::pair<Key, const PrefixNode*> > fringe;
	fringe.push(std::make_pair(Key(), &root_));
	while (!fringe.empty()) {
		const PrefixNode* elem = fringe.front().second;
		const Key& prefix_key = fringe.front().first;
		for (const auto& sub : elem->sub_trees) {
			if (std::binary_search(k.begin(), k.end(), sub.first)) {
				Key new_key(prefix_key.begin(), prefix_key.end());
                new_key.push_back(sub.first);
				result[new_key] = sub.second.value;
				fringe.push(std::make_pair(new_key, &sub.second));
			}
		}
		fringe.pop();
	}
	return result;
}
