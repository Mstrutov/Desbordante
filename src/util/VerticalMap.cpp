#include <exception>
#include <queue>
#include <unordered_set>
#include "VerticalMap.h"
#include "VerticalInfo.h"
#include "DependencyCandidate.h"

namespace util {

template <class Value>
std::shared_ptr<Value> VerticalMap<Value>::SetTrie::associate(
        bitset const& key, size_t nextBit, std::shared_ptr<Value> value) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == boost::dynamic_bitset<>::npos) {
        std::swap(value, value_);
        return value;
    }
    return getOrCreateSubTrie(nextBit)->associate(key, nextBit + 1, std::move(value));
}

template <class Value>
std::shared_ptr<Value const> VerticalMap<Value>::SetTrie::get(bitset const&key, size_t nextBit) const {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == boost::dynamic_bitset<>::npos) {
        return value_;
    }

    auto subtrie = getSubtrie(nextBit);
    if (subtrie == nullptr) return nullptr;
    return subtrie->get(key, nextBit + 1);
}

template <class Value>
std::shared_ptr<Value> VerticalMap<Value>::SetTrie::remove(bitset const&key, size_t nextBit) {
    nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    if (nextBit == bitset::npos) {
        auto removedValue = value_;
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
    return subtries_[index - offset_].get();
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
                                                  std::function<void(bitset const&, std::shared_ptr<Value const>)> collector) const {
    if (value_ != nullptr) {
        collector(bitset(subsetKey), value_);
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
        std::function<bool (bitset const&, std::shared_ptr<Value const>)> const& collector) const {
    if (value_ != nullptr) {
        if (!collector(bitset(subsetKey), value_)) return false;
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
        std::function<bool (bitset const&, std::shared_ptr<Value const>)> const& collector) const {
    if (nextBit != bitset::npos) {
        nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    }
    if (nextBit == bitset::npos) {
        if (value_ != nullptr) {
            if (!collector(bitset(supersetKey), value_)) return false;
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
        std::function<void (bitset const&, std::shared_ptr<Value const>)> const& collector) const {
    if (nextBit != bitset::npos) {
        nextBit = (nextBit == 0 ? key.find_first() : key.find_next(nextBit - 1));
    }
    if (nextBit == bitset::npos) {
        if (value_ != nullptr) {
            collector(bitset(supersetKey), value_);
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
            [&subsetKeys, this](auto& indices, [[maybe_unused]] auto value)
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
            [&entries, this](auto& indices, auto value) {
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
            [&entry, this](auto& indices, auto value) {
                entry = {relation_->getVertical(indices), value};
                return false;
            }
    );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::getAnySubsetEntry(
        const Vertical &vertical,std::function<bool(Vertical const*, std::shared_ptr<Value const>)> const& condition) const {
    typename VerticalMap<Value>::Entry entry;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.collectSubsetKeys(
            vertical.getColumnIndices(),
            0,
            subsetKey,
            [&entry, this, &condition](auto& indices, auto value) {
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
            [&entries, this](auto& indices, auto value) {
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
            [&entry, this](auto& indices, auto value) {
                entry = {relation_->getVertical(indices), value};
                return false;
            }
    );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry
VerticalMap<Value>::getAnySupersetEntry(
        Vertical const &vertical, std::function<bool(Vertical const*, std::shared_ptr<Value const>)> condition) const {
    typename VerticalMap<Value>::Entry entry;
    bitset supersetKey(relation_->getNumColumns());
    setTrie_.collectSupersetKeys(
            vertical.getColumnIndices(),
            0,
            supersetKey,
            [&entry, this, &condition](auto& indices, auto value) {
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
            [&entries, this](auto& indices, auto value) {
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
    setTrie_.traverseEntries(subsetKey, [&keySet, this](auto& k, [[maybe_unused]] auto v) {
        keySet.insert(relation_->getVertical(k));
    });
    return keySet;
}

template<class Value>
std::vector<std::shared_ptr<Value const>> VerticalMap<Value>::values() {
    std::vector<std::shared_ptr<Value const>> values;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&values]([[maybe_unused]] auto& k, auto v) -> void { values.push_back(v); }
    );
    return values;
}

template<class Value>
std::unordered_set<typename VerticalMap<Value>::Entry> VerticalMap<Value>::entrySet() {
    std::unordered_set<typename VerticalMap<Value>::Entry> entrySet;
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&entrySet, this](auto& k, auto v) -> void { entrySet.emplace(relation_->getVertical(k), v); }
    );
    return entrySet;
}

template<class Value>
unsigned int VerticalMap<Value>::removeFromUsageCounter(
        std::unordered_map<Vertical, unsigned int>& usageCounter, const Vertical& key) {
    return usageCounter.erase(key);
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::remove(Vertical const &key) {
    auto removedValue = setTrie_.remove(key.getColumnIndices(), 0);
    if (removedValue != nullptr) size_--;
    return removedValue;
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::remove(const VerticalMap::bitset &key)  {
    auto removedValue = setTrie_.remove(key, 0);
    if (removedValue != nullptr) size_--;
    return removedValue;
}

//comparator is of Compare type - check ascending/descending issues
template<class Value>
void VerticalMap<Value>::shrink(double factor, std::function<bool(Entry, Entry)> const &compare,
                                std::function<bool(Entry)> const &canRemove,
                                [[maybe_unused]] ProfilingContext::ObjectToCache cacheObject) {
    //some logging

    std::priority_queue<Entry, std::vector<Entry>, std::function<bool(Entry, Entry)>> keyQueue(
            compare, std::vector<Entry>(size_));
    bitset subsetKey(relation_->getNumColumns());
    setTrie_.traverseEntries(
            subsetKey,
            [&keyQueue, this, &canRemove](auto& k, auto v) {
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
            [&keyQueue, this, &canRemove, &usageCounter, medianOfUsage](auto& k, auto v) -> void {
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
std::shared_ptr<Value> VerticalMap<Value>::put(Vertical const &key, std::shared_ptr<Value> value) {
    auto oldValue = setTrie_.associate(key.getColumnIndices(), 0, std::move(value));
    if (oldValue == nullptr) size_++;

    return oldValue;
}

template<class Value>
std::shared_ptr<Value const> VerticalMap<Value>::get(Vertical const &key) const {
    return setTrie_.get(key.getColumnIndices(), 0); ;
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::get(Vertical const &key) {
    return std::const_pointer_cast<Value>(setTrie_.get(key.getColumnIndices(), 0)); ;
}

template<class Value>
std::shared_ptr<Value const> VerticalMap<Value>::get(bitset const &key) const {
    return setTrie_.get(key, 0); ;
}

//explicitly instantiate to solve template implementation linking issues
template class VerticalMap<PositionListIndex>;
template class VerticalMap<AgreeSetSample>;
template class VerticalMap<DependencyCandidate>;
template class VerticalMap<VerticalInfo>;
template class VerticalMap<Vertical>;

template<class V>
size_t BlockingVerticalMap<V>::getSize() const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getSize();
}

template<class V>
bool BlockingVerticalMap<V>::isEmpty() const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::isEmpty();
}

template<class V>
std::shared_ptr<V const> BlockingVerticalMap<V>::get(Vertical const& key) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::get(key);
}

template<class V>
std::shared_ptr<V const> BlockingVerticalMap<V>::get(bitset const& key) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::get(key);
}

template<class V>
bool BlockingVerticalMap<V>::containsKey(const Vertical &key) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::containsKey(key);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::put(const Vertical &key, std::shared_ptr<V> value) {
    std::scoped_lock writeLock(readWriteMutex_);
    return VerticalMap<V>::put(key, value);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::remove(const Vertical &key) {
    std::scoped_lock writeLock(readWriteMutex_);
    return VerticalMap<V>::remove(key);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::remove(const bitset &key) {
    std::scoped_lock writeLock(readWriteMutex_);
    return VerticalMap<V>::remove(key);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::get(const Vertical &key) {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::get(key);
}

template<class V>
std::unordered_set<Vertical> BlockingVerticalMap<V>::keySet() {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::keySet();
}

template<class V>
std::vector<std::shared_ptr<V const>> BlockingVerticalMap<V>::values() {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::values();
}

template<class V>
std::unordered_set<typename BlockingVerticalMap<V>::Entry> BlockingVerticalMap<V>::entrySet() {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::entrySet();
}

template<class V>
std::vector<Vertical> BlockingVerticalMap<V>::getSubsetKeys(const Vertical &vertical) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getSubsetKeys(vertical);
}

template<class V>
std::vector<typename BlockingVerticalMap<V>::Entry> BlockingVerticalMap<V>::getSubsetEntries(
        const Vertical &vertical) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getSubsetEntries(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::getAnySubsetEntry(const Vertical &vertical) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getAnySubsetEntry(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::getAnySubsetEntry(const Vertical &vertical,
                                                             const std::function<bool(const Vertical *,
                                                                                      std::shared_ptr<V const>)> &condition) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getAnySubsetEntry(vertical, condition);
}

template<class V>
std::vector<typename BlockingVerticalMap<V>::Entry> BlockingVerticalMap<V>::getSupersetEntries(const Vertical &vertical) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getSupersetEntries(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::getAnySupersetEntry(const Vertical &vertical) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getAnySupersetEntry(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::getAnySupersetEntry(const Vertical &vertical,
                                                               std::function<bool(const Vertical *,
                                                                                  std::shared_ptr<V const>)> condition) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getAnySupersetEntry(vertical, condition);
}

template<class V>
std::vector<typename BlockingVerticalMap<V>::Entry>
BlockingVerticalMap<V>::getRestrictedSupersetEntries(const Vertical &vertical, const Vertical &exclusion) const {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getRestrictedSupersetEntries(vertical, exclusion);
}

template<class V>
bool BlockingVerticalMap<V>::removeSupersetEntries(const Vertical &key) {
    std::scoped_lock writeLock(readWriteMutex_);
    return VerticalMap<V>::removeSupersetEntries(key);
}

template<class V>
bool BlockingVerticalMap<V>::removeSubsetEntries(const Vertical &key) {
    std::scoped_lock writeLock(readWriteMutex_);
    return VerticalMap<V>::removeSubsetEntries(key);
}

template<class V>
void BlockingVerticalMap<V>::shrink(double factor, const std::function<bool(Entry, Entry)> &compare,
                                    const std::function<bool(Entry)> &canRemove,
                                    ProfilingContext::ObjectToCache cacheObject) {
    std::scoped_lock writeLock(readWriteMutex_);
    VerticalMap<V>::shrink(factor, compare, canRemove, cacheObject);
}

template<class V>
void BlockingVerticalMap<V>::shrink(std::unordered_map<Vertical, unsigned int> &usageCounter,
                                    const std::function<bool(Entry)> &canRemove) {
    std::scoped_lock writeLock(readWriteMutex_);
    VerticalMap<V>::shrink(usageCounter, canRemove);
}

template<class V>
long long BlockingVerticalMap<V>::getShrinkInvocations() {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getShrinkInvocations();
}

template<class V>
long long BlockingVerticalMap<V>::getTimeSpentOnShrinking() {
    std::shared_lock readLock(readWriteMutex_);
    return VerticalMap<V>::getTimeSpentOnShrinking();
}

template class BlockingVerticalMap<PositionListIndex>;
template class BlockingVerticalMap<AgreeSetSample>;
template class BlockingVerticalMap<DependencyCandidate>;
template class BlockingVerticalMap<VerticalInfo>;
template class BlockingVerticalMap<Vertical>;

} // namespace util

