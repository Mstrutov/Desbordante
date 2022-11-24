#pragma once

// see ../algorithms/cfd/LICENSE

#include <cmath>
#include <ctime>
#include <bitset>
#include <vector>
#include <algorithm>
#include <iostream>
#include <map>
#include <random>
#include <functional>

struct [[maybe_unused]] SubsetIterator {
    SubsetIterator(int, int);
    std::bitset<64> Next();
    long int GetSubsNumber();

    long int seed;
    const long int max_seed;
    const long int subs_number;
};

[[maybe_unused]] int BinomialCoeff(int n, int k);
[[maybe_unused]] std::vector<std::vector<int>> CartesianProduct(const std::vector<int>&);
[[maybe_unused]] std::vector<std::bitset<32>> AllSubsets(const int);
[[maybe_unused]] std::vector<std::bitset<32>> AllSubsetsIncl(const int);
[[maybe_unused]] void SubsetsLengthK(const int, const int, std::vector<std::bitset<32>>&);
std::vector<int> Range(int, int, int=1);
std::vector<int> Iota(int);

template <typename T>
[[maybe_unused]] typename T::value_type product(const T& items) {
   typename T::value_type prod{};
   for (const typename T::value_type& i : items) {
       prod *= i;
   }
   return prod;
}

template <typename T>
[[maybe_unused]] typename T::value_type implode(const T& collection) {
    typename T::value_type result;
    for (const auto& c : collection) {
        result = Join(result, c);
    }
    return result;
}

template<typename T>
bool IsSubsetOf(const T& sub, const T& super) {
    return std::includes(super.begin(), super.end(), sub.begin(), sub.end());
}

template<typename T>
[[maybe_unused]] bool IsStrictSubsetOf(const T& sub, const T& super) {
    return sub.size() < super.size() && std::includes(super.begin(), super.end(), sub.begin(), sub.end());
}

template <typename T>
[[maybe_unused]] bool ContainsStrictSubsetOf(const T& collection, const typename T::value_type& item) {
    for (const typename T::value_type& s : collection) {
        if (s.size() < item.size() && IsSubsetOf(s, item)) return true;
    }
    return false;
}

template <typename T>
[[maybe_unused]] bool ContainsSubsetOf(const T& collection, const typename T::value_type& item) {
    for (const typename T::value_type& s : collection) {
        if (IsSubsetOf(s, item)) return true;
    }
    return false;
}

template <typename T>
[[maybe_unused]] bool ContainsSupersetOf(const T& collection, const typename T::value_type& item) {
    for (const typename T::value_type& s : collection) {
        if (IsSubsetOf(item, s)) return true;
    }
    return false;
}

template<typename T>
T Subset(const T& items, const std::bitset<64> mask) {
    T sub;
    sub.reserve(mask.count());
    for (unsigned mi = 0; mi < items.size(); mi++) {
        if (mask[mi]) {
            sub.push_back(items[mi]);
        }
    }
    return sub;
}

template<typename T>
T Subset(const T& items, const typename T::value_type leave_out) {
    if (items.size() == 1) return T();
    T sub;
    sub.reserve(items.size() - 1);
    for (unsigned mi = 0; mi < items.size(); mi++) {
        if (items[mi] != leave_out) {
            sub.push_back(items[mi]);
        }
    }
    return sub;
}

template<typename T>
[[maybe_unused]] void InsertSorted(const T& in, const typename T::value_type item, T& out) {
    auto it = in.begin();
    while (it != in.end() && *it < item) { it++; }
    if (it != in.begin()) {
        out.insert(out.begin(), in.begin(), it);
    }
    int offset = (int)(it - in.begin());
    out.insert(out.begin() + offset, item);
    if (in.size() && it != in.end()) {
        out.insert(out.begin() + offset + 1, it, in.end());
    }
}

template<typename T>
T* Intersection(const T* lhs, const T* rhs) {
    T* isect = new T(lhs->size());
    auto it = std::set_intersection(lhs->begin(), lhs->end(), rhs->begin(), rhs->end(), isect->begin());
    isect->resize((int)(it - isect->begin()));
    return isect;
}

template<typename T>
T Intersection(const T& lhs, const T& rhs) {
    T isect(std::min(lhs.size(), rhs.size()));
    auto it = std::set_intersection(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), isect.begin());
    isect.resize((int)(it - isect.begin()));
    return isect;
}

template<typename T>
T* SetDiff(const T* lhs, const T* rhs) {
    T* res = new T();
    std::set_difference(lhs->begin(), lhs->end(), rhs->begin(), rhs->end(), std::inserter(*res, res->begin()));
    return res;
}

template<typename T>
T SetDiff(const T& lhs, const T& rhs) {
    T res;
    std::set_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::inserter(res, res.begin()));
    return res;
}

template<typename T>
T* Join(const T* lhs, const T* rhs) {
    T* uni = new T(lhs->size() + rhs->size());
    auto it = std::set_union(lhs->begin(), lhs->end(), rhs->begin(), rhs->end(), uni->begin());
    uni->resize((int)(it - uni->begin()));
    return uni;
}

template<typename T>
T Join(const T& lhs, const T& rhs) {
    T uni(lhs.size() + rhs.size());
    auto it = std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), uni.begin());
    uni.resize((int)(it - uni.begin()));
    return uni;
}

template<typename T>
T Join(const T& lhs, const typename T::value_type& rhs) {
    T res;
    res.reserve(lhs.size() + 1);
    InsertSorted(lhs, rhs, res);
    return res;
}

template<typename T>
T* Join(const T* lhs, const typename T::value_type& rhs_item) {
    T* res = new T();
    res->reserve(lhs->size() + 1);
    InsertSorted(*lhs, rhs_item, *res);
    return res;
}

template<typename T>
T Remove(const T& lhs, const typename T::value_type& rhs_item) {
    T res;
    res.reserve(lhs.size() - 1);
    const auto& lower = std::lower_bound(lhs.begin(), lhs.end(), rhs_item);
    res.insert(res.begin(), lhs.begin(), lower);
    res.insert(res.begin() + (lower - lhs.begin()), lower + 1, lhs.end());
    return res;
}

template <typename T>
[[maybe_unused]] bool Contains(const T& collection, const typename T::value_type& item) {
    return std::find(collection.begin(), collection.end(), item) != collection.end();
}

template <typename T>
[[maybe_unused]] bool ContainsKey(const T& collection, const typename T::key_type& item) {
    return collection.Find(item) != collection.end();
}

template <typename T>
[[maybe_unused]] typename T::value_type RandElem(const T& collection) {
    std::mt19937 gen(time(0));
    std::uniform_int_distribution<int> dis(0, collection.size()-1);
    return collection[dis(gen)];
}

template <typename T>
[[maybe_unused]] T RandSubset(const T& collection, int size) {
    std::vector<int> indices = Iota(collection.size());
    std::mt19937 gen(time(0));
    std::shuffle(indices.begin(), indices.end(), gen);
    T res(size);
    for (int i = 0; i < size; i++) {
        res[i] = collection[indices[i]];
    }
    return res;
}

template <typename T>
[[maybe_unused]] std::vector<T> RandPartitions(const T& collection, int nr) {
    std::vector<int> indices = Iota(collection.size());
    std::mt19937 gen(time(0));
    std::shuffle(indices.begin(), indices.end(), gen);
    std::vector<T> res(nr);
    int size = collection.size() / nr;
    int j = 0;
    int block = 0;
    for (int i = 0; i < collection.size(); i++) {
        res[j].push_back(collection[indices[i]]);
        block++;
        if (block == size) {
            std::sort(res[j].begin(), res[j].end());
            block = 0;
            j++;
        }
    }
    return res;
}

template <typename T>
bool Has(const T& collection, std::function<bool(typename T::value_type)> f) {
    for (const auto& i : collection) {
        if (f(i)) {
            return true;
        }
    }
    return false;
}

template <typename T>
[[maybe_unused]] T Where(const T& collection, std::function<bool(typename T::value_type)> f) {
    T res;
    for (const auto& i : collection) {
        if (f(i)) {
            res.push_back(i);
        }
    }
    return res;
}

template<typename T>
void Shuffle(T& collection) {
    std::mt19937 gen(time(0));
    std::shuffle(collection.begin(), collection.end(), gen);
}

template<typename T>
void Shuffle(T& collection, unsigned seed) {
    std::mt19937 gen(seed);
    std::shuffle(collection.begin(), collection.end(), gen);
}

template <typename T, typename S>
T Projection(const T& collection, const S& indices) {
    T res;
    for (int i : indices) {
        res.push_back(collection[i]);
    }
    return res;
}

template<typename T>
[[maybe_unused]] T Sorted(const T& collection) {
    T copy = collection;
    std::sort(copy.begin(), copy.end());
    return copy;
}

template<typename T>
[[maybe_unused]] double PrecisionAt(const T& found, const T& expected, int at) {
    int prec = 0;
    for (int i = 0; i < at; i++) {
        int exp_pos = 0;
        for (; exp_pos <= at; exp_pos++) {
            if (exp_pos < at && found[i] == expected[exp_pos]) break;
        }
        prec += std::abs(i - exp_pos);
    }
    return prec;
}

template<typename T>
[[maybe_unused]] std::vector<T> Split(const T& collection, const typename T::value_type& spl) {
    std::vector<T> res(1);
    for (int i : collection) {
        if (i == spl) {
            res.push_back(T());
        }
        else {
            res.back().push_back(i);
        }
    }
    res.pop_back();
    return res;
}

template <typename T>
T GetMaxElem(const std::vector<std::pair<T, int> >& collection) {
    int max = -1;
    int maxI = -1;
    for (int i = 0; i < (int)collection.size(); i++) {
        const auto& elem = collection[i];
        if (elem.second > max) {
            max = elem.second;
            maxI = i;
        }
    }
    return collection[maxI].first;
}
