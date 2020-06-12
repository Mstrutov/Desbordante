//
// Created by maxim on 07.05.2020.
//
#include <exception>
#include <queue>
#include <unordered_set>
#include "VerticalMap.h"
#include "VerticalInfo.h"
#include "DependencyCandidate.h"
/*void linkingProblemsSolver() {
    VerticalMap<std::shared_ptr<Vertical>> obj(nullptr);
}*/

template <class Value>
Value VerticalMap<Value>::SetTrie::associate(bitset const& key, size_t nextBit, Value value) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == dynamic_bitset<>::npos) {
        Value oldValue = value_;
        value_ = value;
        return oldValue;
    }
    return getOrCreateSubTrie(nextBit)->associate(key, nextBit + 1, value);
}

template <class Value>
Value VerticalMap<Value>::SetTrie::get(bitset const&key, size_t nextBit) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == dynamic_bitset<>::npos) {
        return value_;
    }

    std::shared_ptr<SetTrie> subtrie = getSubtrie(nextBit);
    if (subtrie == nullptr) return nullptr;
    return subtrie->get(key, nextBit + 1);
}

template <class Value>
Value VerticalMap<Value>::SetTrie::remove(bitset const&key, size_t nextBit) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == bitset::npos) {
        Value removedValue = value_;
        value_ = nullptr;
        return removedValue;
    }

    std::shared_ptr<SetTrie> subtrie = getSubtrie(nextBit);
    if (subtrie == nullptr) return nullptr;
    Value removedValue = subtrie->remove(key, nextBit + 1);
    if (subtrie->isEmpty()) {
        subtries_[nextBit - offset_] = nullptr;
    }

    return removedValue;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::isEmpty() const {
    if (value_ == nullptr) return false;
    if (subtries_.empty()) return true;
    for (std::shared_ptr<SetTrie> subtrie : subtries_) {
        if (subtrie == nullptr) return false;
    }
    return true;
}

template <class Value>
std::shared_ptr<typename VerticalMap<Value>::SetTrie> VerticalMap<Value>::SetTrie::getOrCreateSubTrie(size_t index) {
    if (subtries_.empty()) {
        subtries_ = std::vector<std::shared_ptr<SetTrie>>(dimension_ - offset_);
    }
    std::shared_ptr<SetTrie> subtrie = getSubtrie(index);
    if (subtrie == nullptr) {
        subtrie = std::make_shared<SetTrie>(index + 1, dimension_);
        subtries_[index - offset_] = subtrie;
    }
    return subtrie;
}

template <class Value>
std::shared_ptr<typename VerticalMap<Value>::SetTrie> VerticalMap<Value>::SetTrie::getSubtrie(size_t index) {
    if (subtries_.empty()) {
        return nullptr;
    }
    if (index < offset_ || index >= dimension_) {
        throw std::runtime_error("Error in getSubtrie: index must be in [offset_; dimension_)");
    }
    return subtries_[index - offset_];

}

template <class Value>
void VerticalMap<Value>::SetTrie::traverseEntries(bitset &subsetKey,
                                                  std::function<void(bitset const&, Value)> collector) {
    if (value_ != nullptr) {
        collector(bitset(subsetKey), value_);
    }   
    for (size_t i = offset_; i < dimension_; i++) {
        std::shared_ptr<SetTrie> subtrie = getSubtrie(i);
        if (subtrie != nullptr) {
            subsetKey.set(i);
            subtrie->traverseEntries(subsetKey, collector);
            subsetKey.reset(i);
        }
    }
}

template <class Value>
bool VerticalMap<Value>::SetTrie::collectSubsetKeys(bitset const& key, size_t nextBit, bitset & subsetKey, std::function<bool (bitset &&, Value)> const& collector) {
    if (value_ != nullptr) {
        if (!collector(bitset(subsetKey), value_)) return false;
    }
    
    for (nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
         nextBit != bitset::npos; 
         nextBit = key.find_next(nextBit)) {
        std::shared_ptr<SetTrie> subtrie = getSubtrie(nextBit);
        if (subtrie != nullptr) {
            subsetKey.set(nextBit);
            if (!subtrie->collectSubsetKeys(key, nextBit + 1, subsetKey, collector)) return false;
            subsetKey.reset(nextBit);
        }
    }
    return true;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::collectSupersetKeys(bitset const& key, size_t nextBit, bitset &supersetKey, std::function<bool (bitset &&, Value)> const& collector) {
    if (nextBit != bitset::npos) {
        nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    }
    if (nextBit == bitset::npos) {
        if (value_ != nullptr) {
            if (!collector(bitset(supersetKey), value_)) return false;
        }
        for (size_t i = offset_; i < dimension_; i++) {
            std::shared_ptr<SetTrie> subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectSupersetKeys(key, nextBit, supersetKey, collector)) return false;
                supersetKey.reset(i);
            }
        }
    }
    else {
        for (size_t i = offset_; i < nextBit; i++) {
            std::shared_ptr<SetTrie> subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectSupersetKeys(key, nextBit, supersetKey, collector)) return false;
                supersetKey.reset(i);
            }
        }
        
        std::shared_ptr<SetTrie> subtrie = getSubtrie(nextBit);
        if (subtrie != nullptr) {
            supersetKey.set(nextBit);
            if (!subtrie->collectSupersetKeys(key, nextBit + 1, supersetKey, collector)) return false;
            supersetKey.reset(nextBit);
        }        
    }
    return true;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::collectRestrictedSupersetKeys(bitset const&key, bitset const& blacklist, size_t nextBit, bitset &supersetKey, std::function<void (bitset &&, Value)>const& collector) {
    if (nextBit != bitset::npos) {
        nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    }
    if (nextBit == bitset::npos) {
        if (value_ != nullptr) {
            collector(bitset(supersetKey), value_);
        }
        for (size_t i = offset_; i < dimension_; i++) {
            if (blacklist.test(i)) continue;
            std::shared_ptr<SetTrie> subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectRestrictedSupersetKeys(key, blacklist, nextBit, supersetKey, collector)) return false;
                supersetKey.reset(i);
            }
        }
    }
    else {
        for (size_t i = offset_; i < nextBit; i++) {
            if (blacklist.test(i)) continue;
            std::shared_ptr<SetTrie> subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectRestrictedSupersetKeys(key, blacklist, nextBit, supersetKey, collector)) return false;
                supersetKey.reset(i);
            }
        }

        std::shared_ptr<SetTrie> subtrie = getSubtrie(nextBit);
        if (subtrie != nullptr) {
            supersetKey.set(nextBit);
            if (!subtrie->collectRestrictedSupersetKeys(key, blacklist, nextBit + 1, supersetKey, collector)) return false;
            supersetKey.reset(nextBit);
        }
    }
    return true;
}

template<class Value>
std::vector<std::shared_ptr<Vertical>> VerticalMap<Value>::getSubsetKeys(const Vertical &vertical) {
    std::vector<std::shared_ptr<Vertical>> subsetKeys;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset subsetKey(relation->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&subsetKeys, relation](auto indices, auto value)
                {
                    subsetKeys.push_back(relation->getVertical(indices));
                    return true;        //??? what's the point in the original code - just continue til the end
                }
            );
    return subsetKeys;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::getSubsetEntries(const Vertical &vertical) {
    std::vector<typename VerticalMap<Value>::Entry> entries;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset subsetKey(relation->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entries, relation](auto indices, auto value) {
                entries.emplace_back(relation->getVertical(indices), value);
                return true;
            }
    );
    return entries;
}

//returns an empty pair if no entry is found
template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySubsetEntry(const Vertical &vertical) {
    typename VerticalMap<Value>::Entry entry;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset subsetKey(relation->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entry, relation](auto indices, auto value) {
                entry = {relation->getVertical(indices), value};
                return false;
            }
    );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySubsetEntry(const Vertical &vertical, std::function<bool(Vertical*, Value)> const& condition) {
    typename VerticalMap<Value>::Entry entry;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset subsetKey(relation->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entry, relation, &condition](auto indices, auto value) {
                auto kv = relation->getVertical(indices);
                if (condition(kv.get(), value)) {
                    entry = {kv, value};
                    return false;
                }
                else {
                    return true;
                }
            }
    );
    return entry;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::getSupersetEntries(Vertical const &vertical) {
    std::vector<typename VerticalMap<Value>::Entry> entries;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset supersetKey(relation->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entries, relation](auto indices, auto value) {
                entries.emplace_back(relation->getVertical(indices), value);
                return true;
            }
    );
    return entries;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySupersetEntry(Vertical const &vertical) {
    typename VerticalMap<Value>::Entry entry;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset supersetKey(relation->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entry, relation](auto indices, auto value) {
                entry = {relation->getVertical(indices), value};
                return false;
            }
    );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry
VerticalMap<Value>::getAnySupersetEntry(Vertical const &vertical, std::function<bool(Vertical*, Value)> condition) {
    typename VerticalMap<Value>::Entry entry;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset supersetKey(relation->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entry, relation, &condition](auto indices, auto value) {
                auto kv = relation->getVertical(indices);
                if (condition(kv.get(), value)) {
                    entry = {kv, value};
                    return false;
                }
                else {
                    return true;
                }
            }
    );
    return entry;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::getRestrictedSupersetEntries(Vertical const &vertical, Vertical const &exclusion) {
    if (vertical.getColumnIndices().intersects(exclusion.getColumnIndices()))
        throw std::runtime_error("Error in getRestrictedSupersetEntries: a vertical shouldn't intersect with a restriction");

    std::vector<typename VerticalMap<Value>::Entry> entries;
    std::shared_ptr<RelationalSchema> relation = relation_.lock();
    bitset supersetKey(relation->getNumColumns());
    setTrie_.collectRestrictedSupersetKeys(
            vertical.getColumnIndices(),
            exclusion.getColumnIndices(),
            0,
            supersetKey,
            [&entries, relation](auto indices, auto value) {
                entries.emplace_back(relation->getVertical(indices), value);
                return true;
            }
    );
    return entries;
}

template<class Value>
bool VerticalMap<Value>::removeSupersetEntries(Vertical const& key) {
    std::vector<typename VerticalMap<Value>::Entry> supersetEntries = getSupersetEntries(key);
    for (auto supersetEntry : supersetEntries) {
        remove(*supersetEntry.first);
    }
    return !supersetEntries.empty();
}

template<class Value>
bool VerticalMap<Value>::removeSubsetEntries(Vertical const& key) {
    std::vector<typename VerticalMap<Value>::Entry> subsetEntries = getSubsetEntries(key);
    for (auto subsetEntry : subsetEntries) {
        remove(*subsetEntry.first);
    }
    return !subsetEntries.empty();
}

template<class Value>
std::unordered_set<std::shared_ptr<Vertical>> VerticalMap<Value>::keySet() {
    std::unordered_set<std::shared_ptr<Vertical>> keySet;
    auto relation_ptr = relation_.lock();

    bitset subsetKey(relation_ptr->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keySet, relation_ptr](auto k, auto v) { keySet.insert(relation_ptr->getVertical(k)); }
            );
    return keySet;
}

template<class Value>
std::vector<Value> VerticalMap<Value>::values() {
    std::vector<Value> values;
    auto relation_ptr = relation_.lock();
    bitset subsetKey(relation_ptr->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&values, relation_ptr](auto k, auto v) -> void { values.push_back(v); }
    );
    return values;
}

template<class Value>
std::unordered_set<typename VerticalMap<Value>::Entry> VerticalMap<Value>::entrySet() {
    std::unordered_set<typename VerticalMap<Value>::Entry> entrySet;
    auto relation_ptr = relation_.lock();
    bitset subsetKey(relation_ptr->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&entrySet, relation_ptr](auto k, auto v) -> void { entrySet.emplace(relation_ptr->getVertical(k), v); }
    );
    return entrySet;
}

template<class Value>
int VerticalMap<Value>::removeFromUsageCounter(std::unordered_map<Vertical, int>& usageCounter, Vertical key) {
    return usageCounter.erase(key);
}

template<class Value>
Value VerticalMap<Value>::remove(Vertical const &key)     {
    Value removedValue = setTrie_.remove(key.getColumnIndices(), 0);
    if (removedValue != nullptr) size_--;
    return removedValue;
}

template<class Value>
Value VerticalMap<Value>::remove(const VerticalMap::bitset &key)  {
    Value removedValue = setTrie_.remove(key, 0);
    if (removedValue != nullptr) size_--;
    return removedValue;
}

//comparator is of Compare type - check ascending/descending issues
template<class Value>
void VerticalMap<Value>::shrink(double factor, std::function<bool(Entry, Entry)> const &compare,
                                std::function<bool(Entry)> const &canRemove, ProfilingContext::ObjectToCache cacheObject) {
    //some logging

    std::priority_queue<Entry, std::vector<Entry>, std::function<bool(Entry, Entry)>> keyQueue(compare, std::vector<Entry>(size_));
    auto relation_ptr = relation_.lock();
    bitset subsetKey(relation_ptr->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keyQueue, relation_ptr, &canRemove](auto k, auto v) {
                if (Entry entry(relation_ptr->getVertical(k), v); canRemove(entry)) {
                    keyQueue.push(entry);
                }
            }
            );
    unsigned int numOfRemoved = 0;
    unsigned int targetSize = size_ * factor;
    while (!keyQueue.empty() && size_ > targetSize) {
        auto key = keyQueue.top().first;
        keyQueue.pop();

        //insert additional logging

        numOfRemoved++;
        remove(*key);
    }
    shrinkInvocations_++;
    timeSpentOnShrinking_+= 1; //haven't implemented time measuring yet
}

template<class Value>
void VerticalMap<Value>::shrink(std::unordered_map<Vertical, int> &usageCounter,
                                std::function<bool(Entry)> const &canRemove) {
    //some logging

    std::vector<int> usageCounters(usageCounter.size());
    for (auto& [first, second] : usageCounter) {
        usageCounters.push_back(second);
    }
    std::sort(usageCounters.begin(), usageCounters.end());
    unsigned int medianOfUsage = usageCounters.size() % 2 == 0 ?
                                 (usageCounters[usageCounters.size() / 2 + 1] + usageCounters[usageCounters.size() / 2]) / 2 :
                                 usageCounters[usageCounters.size() / 2];

    std::queue<Entry> keyQueue;
    auto relation_ptr = relation_.lock();
    bitset subsetKey(relation_ptr->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keyQueue, relation_ptr, &canRemove, &usageCounter, medianOfUsage](auto k, auto v) -> void {
                if (Entry entry(relation_ptr->getVertical(k), v); canRemove(entry) && usageCounter.at(*entry.first) <= medianOfUsage) {
                    keyQueue.push(entry);
                }
            }
    );
    unsigned int numOfRemoved = 0;
    while (!keyQueue.empty()) {
        auto key = keyQueue.front().first;
        keyQueue.pop();

        //insert additional logging

        numOfRemoved++;
        remove(*key);
        removeFromUsageCounter(usageCounter, *key);
    }

    //TODO: what do we want to accomplish here? - looks ok btw
    for (auto& [first, second] : usageCounter) {
        second = 0;
    }

    shrinkInvocations_++;
    timeSpentOnShrinking_+= 1; //haven't implemented time measuring yet
}

template<class Value>
Value VerticalMap<Value>::put(Vertical const &key, Value value) {
    Value oldValue = setTrie_.associate(key.getColumnIndices(), 0, value);
    if (oldValue == nullptr) size_++;

    return oldValue;
}

template<class Value>
Value VerticalMap<Value>::get(Vertical const &key) {
    return setTrie_.get(key.getColumnIndices(), 0); ;
}

template<class Value>
Value VerticalMap<Value>::get(bitset const &key) {
    return setTrie_.get(key, 0); ;
}

//explicitly instantiate to solve template implementation linking issues
template class VerticalMap<std::shared_ptr<PositionListIndex>>;
template class VerticalMap<std::shared_ptr<AgreeSetSample>>;
template class VerticalMap<std::shared_ptr<DependencyCandidate>>;
template class VerticalMap<std::shared_ptr<VerticalInfo>>;
template class VerticalMap<std::shared_ptr<Vertical>>;
