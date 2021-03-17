#include <exception>
#include <queue>
#include <unordered_set>
#include "VerticalMap.h"
#include "VerticalInfo.h"
#include "DependencyCandidate.h"

template <class Value>
std::unique_ptr<Value> VerticalMap<Value>::SetTrie::associate(
        bitset const& key, size_t nextBit, std::unique_ptr<Value> value) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == dynamic_bitset<>::npos) {
        auto oldValue = std::move(value_);
        value_ = std::move(value);
        return oldValue;
    }
    return getOrCreateSubTrie(nextBit)->associate(key, nextBit + 1, std::move(value));
}

template <class Value>
Value const* VerticalMap<Value>::SetTrie::get(bitset const&key, size_t nextBit) const {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == dynamic_bitset<>::npos) {
        return value_.get();
    }

    auto subtrie = getSubtrie(nextBit);
    if (subtrie == nullptr) return nullptr;
    return subtrie->get(key, nextBit + 1);
}

template <class Value>
std::unique_ptr<Value> VerticalMap<Value>::SetTrie::remove(bitset const&key, size_t nextBit) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == bitset::npos) {
        auto removedValue = std::move(value_);
        value_ = nullptr;
        return removedValue;
    }

    auto subtrie = getSubtrie(nextBit);
    if (subtrie == nullptr) return nullptr;
    auto removedValue = subtrie->remove(key, nextBit + 1);
    if (subtrie->isEmpty()) {
        subtries_[nextBit - offset_] = nullptr;
    }

    return removedValue;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::isEmpty() const {
    if (value_ == nullptr) return false;
    return std::all_of(subtries_.begin(), subtries_.end(),
                       [](auto& subtrie_ptr) { return subtrie_ptr == nullptr; });
}

template <class Value>
typename VerticalMap<Value>::SetTrie* VerticalMap<Value>::SetTrie::getOrCreateSubTrie(size_t index) {
    if (subtries_.empty()) {
        subtries_ = std::vector<std::unique_ptr<SetTrie>>(dimension_ - offset_);
    }
    SetTrie* subtrie = getSubtrie(index);
    if (subtrie == nullptr) {
        subtries_[index - offset_] = std::make_unique<SetTrie>(index + 1, dimension_);
    }
    return subtrie;
}

template <class Value>
typename VerticalMap<Value>::SetTrie* VerticalMap<Value>::SetTrie::getSubtrie(size_t index) {
    return const_cast<VerticalMap<Value>::SetTrie*>(
            (const_cast<VerticalMap<Value>::SetTrie const*>(this))->getSubtrie(index)
            );
}

template <class Value>
typename VerticalMap<Value>::SetTrie const* VerticalMap<Value>::SetTrie::getSubtrie(size_t index) const {
    if (subtries_.empty()) {
        return nullptr;
    }
    if (index < offset_ || index >= dimension_) {
        throw std::runtime_error("Error in getSubtrie: index must be in [offset_; dimension_)");
    }
    return subtries_[index - offset_].get();
}

template <class Value>
void VerticalMap<Value>::SetTrie::traverseEntries(bitset &subsetKey,
                                                  std::function<void(bitset const&, Value const*)> collector) const {
    if (value_ != nullptr) {
        collector(bitset(subsetKey), value_.get());
    }   
    for (size_t i = offset_; i < dimension_; i++) {
        auto subtrie = getSubtrie(i);
        if (subtrie != nullptr) {
            subsetKey.set(i);
            subtrie->traverseEntries(subsetKey, collector);
            subsetKey.reset(i);
        }
    }
}

template <class Value>
bool VerticalMap<Value>::SetTrie::collectSubsetKeys(
        bitset const& key, size_t nextBit, bitset & subsetKey,
        std::function<bool (bitset &&, Value const*)> const& collector) const {
    if (value_ != nullptr) {
        if (!collector(bitset(subsetKey), value_.get())) return false;
    }
    
    for (nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
         nextBit != bitset::npos; 
         nextBit = key.find_next(nextBit)) {
        auto subtrie = getSubtrie(nextBit);
        if (subtrie != nullptr) {
            subsetKey.set(nextBit);
            if (!subtrie->collectSubsetKeys(key, nextBit + 1, subsetKey, collector)) return false;
            subsetKey.reset(nextBit);
        }
    }
    return true;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::collectSupersetKeys(
        bitset const& key, size_t nextBit, bitset &supersetKey,
        std::function<bool (bitset &&, Value const*)> const& collector) const {
    if (nextBit != bitset::npos) {
        nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    }
    if (nextBit == bitset::npos) {
        if (value_ != nullptr) {
            if (!collector(bitset(supersetKey), value_.get())) return false;
        }
        for (size_t i = offset_; i < dimension_; i++) {
            auto subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectSupersetKeys(key, nextBit, supersetKey, collector)) return false;
                supersetKey.reset(i);
            }
        }
    }
    else {
        for (size_t i = offset_; i < nextBit; i++) {
            auto subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectSupersetKeys(key, nextBit, supersetKey, collector)) return false;
                supersetKey.reset(i);
            }
        }
        
        auto subtrie = getSubtrie(nextBit);
        if (subtrie != nullptr) {
            supersetKey.set(nextBit);
            if (!subtrie->collectSupersetKeys(key, nextBit + 1, supersetKey, collector)) return false;
            supersetKey.reset(nextBit);
        }        
    }
    return true;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::collectRestrictedSupersetKeys(
        bitset const&key, bitset const& blacklist, size_t nextBit, bitset &supersetKey,
        std::function<void (bitset &&, Value const*)> const& collector) const {
    if (nextBit != bitset::npos) {
        nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    }
    if (nextBit == bitset::npos) {
        if (value_ != nullptr) {
            collector(bitset(supersetKey), value_.get());
        }
        for (size_t i = offset_; i < dimension_; i++) {
            if (blacklist.test(i)) continue;
            auto subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectRestrictedSupersetKeys(key, blacklist, nextBit, supersetKey, collector)) {
                    return false;
                }
                supersetKey.reset(i);
            }
        }
    }
    else {
        for (size_t i = offset_; i < nextBit; i++) {
            if (blacklist.test(i)) continue;
            auto subtrie = getSubtrie(i);
            if (subtrie != nullptr) {
                supersetKey.set(i);
                if (!subtrie->collectRestrictedSupersetKeys(key, blacklist, nextBit, supersetKey, collector)) {
                    return false;
                }
                supersetKey.reset(i);
            }
        }

        auto subtrie = getSubtrie(nextBit);
        if (subtrie != nullptr) {
            supersetKey.set(nextBit);
            if (!subtrie->collectRestrictedSupersetKeys(key, blacklist, nextBit + 1, supersetKey, collector)) {
                return false;
            }
            supersetKey.reset(nextBit);
        }
    }
    return true;
}

template<class Value>
std::vector<Vertical> VerticalMap<Value>::getSubsetKeys(Vertical const& vertical) const {
    std::vector<Vertical> subsetKeys;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&subsetKeys, this](auto indices, auto value)
                {
                    subsetKeys.push_back(relation_->getVertical(indices));
                    return true;
                }
            );
    return subsetKeys;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::getSubsetEntries(const Vertical &vertical) const {
    std::vector<typename VerticalMap<Value>::Entry> entries;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entries, this](auto indices, auto value) {
                entries.emplace_back(relation_->getVertical(indices), value);
                return true;
            }
    );
    return entries;
}

//returns an empty pair if no entry is found
template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySubsetEntry(Vertical const& vertical) const {
    typename VerticalMap<Value>::Entry entry;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entry, this](auto indices, auto value) {
                entry = {relation_->getVertical(indices), value};
                return false;
            }
    );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySubsetEntry(
        const Vertical &vertical,std::function<bool(Vertical const*, Value const*)> const& condition) const {
    typename VerticalMap<Value>::Entry entry;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entry, this, &condition](auto indices, auto value) {
                auto kv = relation_->getVertical(indices);
                if (condition(&kv, value)) {
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
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::getSupersetEntries(Vertical const &vertical) const {
    std::vector<typename VerticalMap<Value>::Entry> entries;
    bitset supersetKey(relation_->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entries, this](auto indices, auto value) {
                entries.emplace_back(relation_->getVertical(indices), value);
                return true;
            }
    );
    return entries;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySupersetEntry(Vertical const &vertical) const {
    typename VerticalMap<Value>::Entry entry;
    bitset supersetKey(relation_->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entry, this](auto indices, auto value) {
                entry = {relation_->getVertical(indices), value};
                return false;
            }
    );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry
VerticalMap<Value>::getAnySupersetEntry(
        Vertical const &vertical, std::function<bool(Vertical const*, Value const*)> condition) const {
    typename VerticalMap<Value>::Entry entry;
    bitset supersetKey(relation_->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entry, this, &condition](auto indices, auto value) {
                auto kv = relation_->getVertical(indices);
                if (condition(&kv, value)) {
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
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::getRestrictedSupersetEntries(
        Vertical const &vertical, Vertical const &exclusion) const {
    if (vertical.getColumnIndices().intersects(exclusion.getColumnIndices()))
        throw std::runtime_error("Error in getRestrictedSupersetEntries: a vertical shouldn't intersect with a restriction");

    std::vector<typename VerticalMap<Value>::Entry> entries;
    bitset supersetKey(relation_->getNumColumns());
    setTrie_.collectRestrictedSupersetKeys(
            vertical.getColumnIndices(),
            exclusion.getColumnIndices(),
            0,
            supersetKey,
            [&entries, this](auto indices, auto value) {
                entries.emplace_back(relation_->getVertical(indices), value);
                return true;
            }
    );
    return entries;
}

template<class Value>
bool VerticalMap<Value>::removeSupersetEntries(Vertical const& key) {
    std::vector<typename VerticalMap<Value>::Entry> supersetEntries = getSupersetEntries(key);
    for (auto supersetEntry : supersetEntries) {
        remove(supersetEntry.first);
    }
    return !supersetEntries.empty();
}

template<class Value>
bool VerticalMap<Value>::removeSubsetEntries(Vertical const& key) {
    std::vector<typename VerticalMap<Value>::Entry> subsetEntries = getSubsetEntries(key);
    for (auto subsetEntry : subsetEntries) {
        remove(subsetEntry.first);
    }
    return !subsetEntries.empty();
}

template<class Value>
std::unordered_set<Vertical> VerticalMap<Value>::keySet() {
    std::unordered_set<Vertical> keySet;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keySet, this](auto k, auto v) { keySet.insert(relation_->getVertical(k)); }
            );
    return keySet;
}

template<class Value>
std::vector<Value const*> VerticalMap<Value>::values() {
    std::vector<Value const*> values;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&values](auto k, auto v) -> void { values.push_back(v); }
    );
    return values;
}

template<class Value>
std::unordered_set<typename VerticalMap<Value>::Entry> VerticalMap<Value>::entrySet() {
    std::unordered_set<typename VerticalMap<Value>::Entry> entrySet;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&entrySet, this](auto k, auto v) -> void { entrySet.emplace(relation_->getVertical(k), v); }
    );
    return entrySet;
}

template<class Value>
unsigned int VerticalMap<Value>::removeFromUsageCounter(
        std::unordered_map<Vertical, unsigned int>& usageCounter, const Vertical& key) {
    return usageCounter.erase(key);
}

template<class Value>
std::unique_ptr<Value> VerticalMap<Value>::remove(Vertical const &key) {
    auto removedValue = setTrie_.remove(key.getColumnIndices(), 0);
    if (removedValue != nullptr) size_--;
    return removedValue;
}

template<class Value>
std::unique_ptr<Value> VerticalMap<Value>::remove(const VerticalMap::bitset &key)  {
    auto removedValue = setTrie_.remove(key, 0);
    if (removedValue != nullptr) size_--;
    return removedValue;
}

//comparator is of Compare type - check ascending/descending issues
template<class Value>
void VerticalMap<Value>::shrink(double factor, std::function<bool(Entry, Entry)> const &compare,
                                std::function<bool(Entry)> const &canRemove, ProfilingContext::ObjectToCache cacheObject) {
    //some logging

    std::priority_queue<Entry, std::vector<Entry>, std::function<bool(Entry, Entry)>> keyQueue(
            compare, std::vector<Entry>(size_));
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keyQueue, this, &canRemove](auto k, auto v) {
                if (Entry entry(relation_->getVertical(k), v); canRemove(entry)) {
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
        remove(key);
    }
    shrinkInvocations_++;
    timeSpentOnShrinking_+= 1; //haven't implemented time measuring yet
}

template<class Value>
void VerticalMap<Value>::shrink(std::unordered_map<Vertical, unsigned int> &usageCounter,
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
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keyQueue, this, &canRemove, &usageCounter, medianOfUsage](auto k, auto v) -> void {
                if (Entry entry(relation_->getVertical(k), v); canRemove(entry) && usageCounter.at(entry.first) <= medianOfUsage) {
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
        remove(key);
        removeFromUsageCounter(usageCounter, key);
    }

    //TODO: what do we want to accomplish here? - looks ok btw
    for (auto& [first, second] : usageCounter) {
        second = 0;
    }

    shrinkInvocations_++;
    timeSpentOnShrinking_+= 1; //haven't implemented time measuring yet
}

template<class Value>
std::unique_ptr<Value> VerticalMap<Value>::put(Vertical const &key, std::unique_ptr<Value> value) {
    auto oldValue = setTrie_.associate(key.getColumnIndices(), 0, std::move(value));
    if (oldValue == nullptr) size_++;

    return oldValue;
}

template<class Value>
Value const* VerticalMap<Value>::get(Vertical const &key) const {
    return setTrie_.get(key.getColumnIndices(), 0); ;
}

template<class Value>
Value const* VerticalMap<Value>::get(bitset const &key) const {
    return setTrie_.get(key, 0); ;
}

//explicitly instantiate to solve template implementation linking issues
template class VerticalMap<PositionListIndex>;
template class VerticalMap<AgreeSetSample>;
template class VerticalMap<DependencyCandidate>;
template class VerticalMap<VerticalInfo>;
template class VerticalMap<Vertical>;
