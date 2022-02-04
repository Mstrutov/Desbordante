#include <exception>
#include <queue>
#include <unordered_set>
#include "VerticalMap.h"
#include "VerticalInfo.h"
#include "DependencyCandidate.h"

namespace util {

template <class Value>
std::shared_ptr<Value> VerticalMap<Value>::SetTrie::Associate(
    Bitset const& key, size_t next_bit, std::shared_ptr<Value> value) {
  next_bit = (next_bit == 0 ? key.find_first() : key.find_next(next_bit - 1));
    if (next_bit == boost::dynamic_bitset<>::npos) {
        std::swap(value, value_);
        return value;
    }
    return GetOrCreateSubTrie(next_bit)->Associate(key, next_bit + 1, std::move(value));
}

template <class Value>
std::shared_ptr<Value const> VerticalMap<Value>::SetTrie::Get(Bitset const&key, size_t next_bit) const {
  next_bit = (next_bit == 0 ? key.find_first() : key.find_next(next_bit - 1));
    if (next_bit == boost::dynamic_bitset<>::npos) {
        return value_;
    }

    auto subtrie = GetSubtrie(next_bit);
    if (subtrie == nullptr) return nullptr;
    return subtrie->Get(key, next_bit + 1);
}

template <class Value>
std::shared_ptr<Value> VerticalMap<Value>::SetTrie::Remove(Bitset const&key, size_t next_bit) {
  next_bit = (next_bit == 0 ? key.find_first() : key.find_next(next_bit - 1));
    if (next_bit == Bitset::npos) {
        auto removed_value = value_;
        value_ = nullptr;
        return removed_value;
    }

    auto subtrie = GetSubtrie(next_bit);
    if (subtrie == nullptr) return nullptr;
    auto removed_value = subtrie->Remove(key, next_bit + 1);
    if (subtrie->IsEmpty()) {
        subtries_[next_bit - offset_] = nullptr;
    }

    return removed_value;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::IsEmpty() const {
    if (value_ == nullptr) return false;
    return std::all_of(subtries_.begin(), subtries_.end(),
                       [](auto& subtrie_ptr) { return subtrie_ptr == nullptr; });
}

template <class Value>
typename VerticalMap<Value>::SetTrie* VerticalMap<Value>::SetTrie::GetOrCreateSubTrie(size_t index) {
    if (subtries_.empty()) {
        subtries_ = std::vector<std::unique_ptr<SetTrie>>(dimension_ - offset_);
    }
    SetTrie* subtrie = GetSubtrie(index);
    if (subtrie == nullptr) {
        subtries_[index - offset_] = std::make_unique<SetTrie>(index + 1, dimension_);
    }
    return subtries_[index - offset_].get();
}

template <class Value>
typename VerticalMap<Value>::SetTrie* VerticalMap<Value>::SetTrie::GetSubtrie(size_t index) {
    return const_cast<VerticalMap<Value>::SetTrie*>(
        (const_cast<VerticalMap<Value>::SetTrie const *>(this))->GetSubtrie(index)
            );
}

template <class Value>
typename VerticalMap<Value>::SetTrie const* VerticalMap<Value>::SetTrie::GetSubtrie(size_t index) const {
    if (subtries_.empty()) {
        return nullptr;
    }
    if (index < offset_ || index >= dimension_) {
        throw std::runtime_error("Error in GetSubtrie: index must be in [offset_; dimension_)");
    }
    return subtries_[index - offset_].get();
}

template <class Value>
void VerticalMap<Value>::SetTrie::TraverseEntries(Bitset &subset_key,
                                                  std::function<void(Bitset const&, std::shared_ptr<Value const>)> collector) const {
    if (value_ != nullptr) {
        collector(Bitset(subset_key), value_);
    }   
    for (size_t i = offset_; i < dimension_; i++) {
        auto subtrie = GetSubtrie(i);
        if (subtrie != nullptr) {
            subset_key.set(i);
            subtrie->TraverseEntries(subset_key, collector);
            subset_key.reset(i);
        }
    }
}

template <class Value>
bool VerticalMap<Value>::SetTrie::CollectSubsetKeys(
    Bitset const& key, size_t next_bit, Bitset & subset_key,
    std::function<bool (Bitset const&, std::shared_ptr<Value const>)> const& collector) const {
    if (value_ != nullptr) {
        if (!collector(Bitset(subset_key), value_)) return false;
    }
    
    for (next_bit = (next_bit == 0 ? key.find_first() : key.find_next(next_bit - 1));
         next_bit != Bitset::npos;
         next_bit = key.find_next(next_bit)) {
        auto subtrie = GetSubtrie(next_bit);
        if (subtrie != nullptr) {
            subset_key.set(next_bit);
            if (!subtrie->CollectSubsetKeys(key, next_bit + 1, subset_key, collector)) return false;
            subset_key.reset(next_bit);
        }
    }
    return true;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::CollectSupersetKeys(
    Bitset const& key, size_t next_bit, Bitset &superset_key,
    std::function<bool (Bitset const&, std::shared_ptr<Value const>)> const& collector) const {
    if (next_bit != Bitset::npos) {
      next_bit = (next_bit == 0 ? key.find_first() : key.find_next(next_bit - 1));
    }
    if (next_bit == Bitset::npos) {
        if (value_ != nullptr) {
            if (!collector(Bitset(superset_key), value_)) return false;
        }
        for (size_t i = offset_; i < dimension_; i++) {
            auto subtrie = GetSubtrie(i);
            if (subtrie != nullptr) {
                superset_key.set(i);
                if (!subtrie->CollectSupersetKeys(key, next_bit, superset_key, collector)) return false;
                superset_key.reset(i);
            }
        }
    }
    else {
        for (size_t i = offset_; i < next_bit; i++) {
            auto subtrie = GetSubtrie(i);
            if (subtrie != nullptr) {
                superset_key.set(i);
                if (!subtrie->CollectSupersetKeys(key, next_bit, superset_key, collector)) return false;
                superset_key.reset(i);
            }
        }
        
        auto subtrie = GetSubtrie(next_bit);
        if (subtrie != nullptr) {
            superset_key.set(next_bit);
            if (!subtrie->CollectSupersetKeys(key, next_bit + 1, superset_key, collector)) return false;
            superset_key.reset(next_bit);
        }        
    }
    return true;
}

template <class Value>
bool VerticalMap<Value>::SetTrie::CollectRestrictedSupersetKeys(
    Bitset const&key, Bitset const& blacklist, size_t next_bit, Bitset &superset_key,
    std::function<void (Bitset const&, std::shared_ptr<Value const>)> const& collector) const {
    if (next_bit != Bitset::npos) {
      next_bit = (next_bit == 0 ? key.find_first() : key.find_next(next_bit - 1));
    }
    if (next_bit == Bitset::npos) {
        if (value_ != nullptr) {
            collector(Bitset(superset_key), value_);
        }
        for (size_t i = offset_; i < dimension_; i++) {
            if (blacklist.test(i)) continue;
            auto subtrie = GetSubtrie(i);
            if (subtrie != nullptr) {
                superset_key.set(i);
                if (!subtrie->CollectRestrictedSupersetKeys(key, blacklist, next_bit, superset_key, collector)) {
                    return false;
                }
                superset_key.reset(i);
            }
        }
    }
    else {
        for (size_t i = offset_; i < next_bit; i++) {
            if (blacklist.test(i)) continue;
            auto subtrie = GetSubtrie(i);
            if (subtrie != nullptr) {
                superset_key.set(i);
                if (!subtrie->CollectRestrictedSupersetKeys(key, blacklist, next_bit, superset_key, collector)) {
                    return false;
                }
                superset_key.reset(i);
            }
        }

        auto subtrie = GetSubtrie(next_bit);
        if (subtrie != nullptr) {
            superset_key.set(next_bit);
            if (!subtrie->CollectRestrictedSupersetKeys(key, blacklist, next_bit + 1, superset_key, collector)) {
                return false;
            }
            superset_key.reset(next_bit);
        }
    }
    return true;
}

template<class Value>
std::vector<Vertical> VerticalMap<Value>::GetSubsetKeys(Vertical const& vertical) const {
    std::vector<Vertical> subset_keys;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.CollectSubsetKeys(
      vertical.GetColumnIndices(),
      0,
      subset_key,
      [&subset_keys, this](auto &indices, [[maybe_unused]] auto value) {
        subset_keys.push_back(relation_->GetVertical(indices));
        return true;
      }
  );
    return subset_keys;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::GetSubsetEntries(const Vertical &vertical) const {
    std::vector<typename VerticalMap<Value>::Entry> entries;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.CollectSubsetKeys(
      vertical.GetColumnIndices(),
      0,
      subset_key,
      [&entries, this](auto &indices, auto value) {
        entries.emplace_back(relation_->GetVertical(indices), value);
        return true;
      }
  );
    return entries;
}

//returns an empty pair if no entry is found
template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::GetAnySubsetEntry(Vertical const& vertical) const {
    typename VerticalMap<Value>::Entry entry;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.CollectSubsetKeys(
      vertical.GetColumnIndices(),
      0,
      subset_key,
      [&entry, this](auto &indices, auto value) {
        entry = {relation_->GetVertical(indices), value};
        return false;
      }
  );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::GetAnySubsetEntry(
        const Vertical &vertical,std::function<bool(Vertical const*, std::shared_ptr<Value const>)> const& condition) const {
    typename VerticalMap<Value>::Entry entry;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.CollectSubsetKeys(
      vertical.GetColumnIndices(),
      0,
      subset_key,
      [&entry, this, &condition](auto &indices, auto value) {
        auto kv = relation_->GetVertical(indices);
        if (condition(&kv, value)) {
          entry = {kv, value};
          return false;
        } else {
          return true;
        }
      }
  );
    return entry;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::GetSupersetEntries(Vertical const &vertical) const {
    std::vector<typename VerticalMap<Value>::Entry> entries;
    Bitset superset_key(relation_->GetNumColumns());
  set_trie_.CollectSupersetKeys(
      vertical.GetColumnIndices(),
      0,
      superset_key,
      [&entries, this](auto &indices, auto value) {
        entries.emplace_back(relation_->GetVertical(indices), value);
        return true;
      }
  );
    return entries;
}

template<class Value>
typename VerticalMap<Value>::Entry VerticalMap<Value>::GetAnySupersetEntry(Vertical const &vertical) const {
    typename VerticalMap<Value>::Entry entry;
    Bitset superset_key(relation_->GetNumColumns());
  set_trie_.CollectSupersetKeys(
      vertical.GetColumnIndices(),
      0,
      superset_key,
      [&entry, this](auto &indices, auto value) {
        entry = {relation_->GetVertical(indices), value};
        return false;
      }
  );
    return entry;
}

template<class Value>
typename VerticalMap<Value>::Entry
VerticalMap<Value>::GetAnySupersetEntry(
        Vertical const &vertical, std::function<bool(Vertical const*, std::shared_ptr<Value const>)> condition) const {
    typename VerticalMap<Value>::Entry entry;
    Bitset superset_key(relation_->GetNumColumns());
  set_trie_.CollectSupersetKeys(
      vertical.GetColumnIndices(),
      0,
      superset_key,
      [&entry, this, &condition](auto &indices, auto value) {
        auto kv = relation_->GetVertical(indices);
        if (condition(&kv, value)) {
          entry = {kv, value};
          return false;
        } else {
          return true;
        }
      }
  );
    return entry;
}

template<class Value>
std::vector<typename VerticalMap<Value>::Entry> VerticalMap<Value>::GetRestrictedSupersetEntries(
        Vertical const &vertical, Vertical const &exclusion) const {
    if (vertical.GetColumnIndices().intersects(exclusion.GetColumnIndices()))
        throw std::runtime_error("Error in GetRestrictedSupersetEntries: a vertical shouldn't intersect with a restriction");

    std::vector<typename VerticalMap<Value>::Entry> entries;
    Bitset superset_key(relation_->GetNumColumns());
  set_trie_.CollectRestrictedSupersetKeys(
      vertical.GetColumnIndices(),
      exclusion.GetColumnIndices(),
      0,
      superset_key,
      [&entries, this](auto &indices, auto value) {
        entries.emplace_back(relation_->GetVertical(indices), value);
        return true;
      }
  );
    return entries;
}

template<class Value>
bool VerticalMap<Value>::RemoveSupersetEntries(Vertical const& key) {
    std::vector<typename VerticalMap<Value>::Entry> superset_entries = GetSupersetEntries(key);
    for (auto superset_entry : superset_entries) {
      Remove(superset_entry.first);
    }
    return !superset_entries.empty();
}

template<class Value>
bool VerticalMap<Value>::RemoveSubsetEntries(Vertical const& key) {
    std::vector<typename VerticalMap<Value>::Entry> subset_entries = GetSubsetEntries(key);
    for (auto subset_entry : subset_entries) {
      Remove(subset_entry.first);
    }
    return !subset_entries.empty();
}

template<class Value>
std::unordered_set<Vertical> VerticalMap<Value>::KeySet() {
    std::unordered_set<Vertical> key_set;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.TraverseEntries(subset_key, [&key_set, this](auto &k, [[maybe_unused]] auto v) {
    key_set.insert(relation_->GetVertical(k));
  });
    return key_set;
}

template<class Value>
std::vector<std::shared_ptr<Value const>> VerticalMap<Value>::Values() {
    std::vector<std::shared_ptr<Value const>> values;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.TraverseEntries(
      subset_key,
      [&values]([[maybe_unused]] auto &k, auto v) -> void { values.push_back(v); }
  );
    return values;
}

template<class Value>
std::unordered_set<typename VerticalMap<Value>::Entry> VerticalMap<Value>::EntrySet() {
    std::unordered_set<typename VerticalMap<Value>::Entry> entry_set;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.TraverseEntries(
      subset_key,
      [&entry_set, this](auto &k, auto v) -> void { entry_set.emplace(relation_->GetVertical(k), v); }
  );
    return entry_set;
}

template<class Value>
unsigned int VerticalMap<Value>::RemoveFromUsageCounter(
    std::unordered_map<Vertical, unsigned int>& usage_counter, const Vertical& key) {
    return usage_counter.erase(key);
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::Remove(Vertical const &key) {
    auto removed_value = set_trie_.Remove(key.GetColumnIndices(), 0);
    if (removed_value != nullptr) size_--;
    return removed_value;
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::Remove(const VerticalMap::Bitset &key)  {
    auto removed_value = set_trie_.Remove(key, 0);
    if (removed_value != nullptr) size_--;
    return removed_value;
}

//comparator is of Compare type - check ascending/descending issues
template<class Value>
void VerticalMap<Value>::Shrink(double factor, std::function<bool(Entry, Entry)> const &compare,
                                std::function<bool(Entry)> const &can_remove,
                                [[maybe_unused]] ProfilingContext::ObjectToCache cache_object) {
    //some logging

    std::priority_queue<Entry, std::vector<Entry>, std::function<bool(Entry, Entry)>> key_queue(
            compare, std::vector<Entry>(size_));
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.TraverseEntries(
      subset_key,
      [&key_queue, this, &can_remove](auto &k, auto v) {
        if (Entry entry(relation_->GetVertical(k), v); can_remove(entry)) {
          key_queue.push(entry);
        }
      }
  );
    unsigned int num_of_removed = 0;
    unsigned int target_size = size_ * factor;
    while (!key_queue.empty() && size_ > target_size) {
        auto key = key_queue.top().first;
        key_queue.pop();

        //insert additional logging

        num_of_removed++;
      Remove(key);
    }
    shrink_invocations_++;
  time_spent_on_shrinking_+= 1; //haven't implemented time measuring yet
}

template<class Value>
void VerticalMap<Value>::Shrink(std::unordered_map<Vertical, unsigned int> &usage_counter,
                                std::function<bool(Entry)> const &can_remove) {
    //some logging

    std::vector<int> usage_counters(usage_counter.size());
    for (auto& [first, second] : usage_counter) {
        usage_counters.push_back(second);
    }
    std::sort(usage_counters.begin(), usage_counters.end());
    unsigned int median_of_usage = usage_counters.size() % 2 == 0 ?
                                   (usage_counters[usage_counters.size() / 2 + 1] + usage_counters[usage_counters.size() / 2]) / 2 :
                                   usage_counters[usage_counters.size() / 2];

    std::queue<Entry> key_queue;
    Bitset subset_key(relation_->GetNumColumns());
  set_trie_.TraverseEntries(
      subset_key,
      [&key_queue, this, &can_remove, &usage_counter, median_of_usage](auto &k, auto v) -> void {
        if (Entry entry(relation_->GetVertical(k), v); can_remove(entry)
            && usage_counter.at(entry.first) <= median_of_usage) {
          key_queue.push(entry);
        }
      }
  );
    unsigned int num_of_removed = 0;
    while (!key_queue.empty()) {
        auto key = key_queue.front().first;
        key_queue.pop();

        //insert additional logging

        num_of_removed++;
      Remove(key);
      RemoveFromUsageCounter(usage_counter, key);
    }

    //TODO: what do we want to accomplish here? - looks ok btw
    for (auto& [first, second] : usage_counter) {
        second = 0;
    }

    shrink_invocations_++;
  time_spent_on_shrinking_+= 1; //haven't implemented time measuring yet
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::Put(Vertical const &key, std::shared_ptr<Value> value) {
    auto old_value = set_trie_.Associate(key.GetColumnIndices(), 0, std::move(value));
    if (old_value == nullptr) size_++;

    return old_value;
}

template<class Value>
std::shared_ptr<Value const> VerticalMap<Value>::Get(Vertical const &key) const {
    return set_trie_.Get(key.GetColumnIndices(), 0); ;
}

template<class Value>
std::shared_ptr<Value> VerticalMap<Value>::Get(Vertical const &key) {
    return std::const_pointer_cast<Value>(set_trie_.Get(key.GetColumnIndices(), 0)); ;
}

template<class Value>
std::shared_ptr<Value const> VerticalMap<Value>::Get(Bitset const &key) const {
    return set_trie_.Get(key, 0); ;
}

//explicitly instantiate to solve template implementation linking issues
template class VerticalMap<PositionListIndex>;
template class VerticalMap<AgreeSetSample>;
template class VerticalMap<DependencyCandidate>;
template class VerticalMap<VerticalInfo>;
template class VerticalMap<Vertical>;

template<class V>
size_t BlockingVerticalMap<V>::GetSize() const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetSize();
}

template<class V>
bool BlockingVerticalMap<V>::IsEmpty() const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::IsEmpty();
}

template<class V>
std::shared_ptr<V const> BlockingVerticalMap<V>::Get(Vertical const& key) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::Get(key);
}

template<class V>
std::shared_ptr<V const> BlockingVerticalMap<V>::Get(Bitset const& key) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::Get(key);
}

template<class V>
bool BlockingVerticalMap<V>::ContainsKey(const Vertical &key) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::ContainsKey(key);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::Put(const Vertical &key, std::shared_ptr<V> value) {
    std::scoped_lock write_lock(read_write_mutex_);
    return VerticalMap<V>::Put(key, value);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::Remove(const Vertical &key) {
    std::scoped_lock write_lock(read_write_mutex_);
    return VerticalMap<V>::Remove(key);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::Remove(const Bitset &key) {
    std::scoped_lock write_lock(read_write_mutex_);
    return VerticalMap<V>::Remove(key);
}

template<class V>
std::shared_ptr<V> BlockingVerticalMap<V>::Get(const Vertical &key) {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::Get(key);
}

template<class V>
std::unordered_set<Vertical> BlockingVerticalMap<V>::KeySet() {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::KeySet();
}

template<class V>
std::vector<std::shared_ptr<V const>> BlockingVerticalMap<V>::Values() {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::Values();
}

template<class V>
std::unordered_set<typename BlockingVerticalMap<V>::Entry> BlockingVerticalMap<V>::EntrySet() {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::EntrySet();
}

template<class V>
std::vector<Vertical> BlockingVerticalMap<V>::GetSubsetKeys(const Vertical &vertical) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetSubsetKeys(vertical);
}

template<class V>
std::vector<typename BlockingVerticalMap<V>::Entry> BlockingVerticalMap<V>::GetSubsetEntries(
        const Vertical &vertical) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetSubsetEntries(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::GetAnySubsetEntry(const Vertical &vertical) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetAnySubsetEntry(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::GetAnySubsetEntry(const Vertical &vertical,
                                                                                 const std::function<bool(const Vertical *,
                                                                                      std::shared_ptr<V const>)> &condition) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetAnySubsetEntry(vertical, condition);
}

template<class V>
std::vector<typename BlockingVerticalMap<V>::Entry> BlockingVerticalMap<V>::GetSupersetEntries(const Vertical &vertical) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetSupersetEntries(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::GetAnySupersetEntry(const Vertical &vertical) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetAnySupersetEntry(vertical);
}

template<class V>
typename BlockingVerticalMap<V>::Entry BlockingVerticalMap<V>::GetAnySupersetEntry(const Vertical &vertical,
                                                                                   std::function<bool(const Vertical *,
                                                                                  std::shared_ptr<V const>)> condition) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetAnySupersetEntry(vertical, condition);
}

template<class V>
std::vector<typename BlockingVerticalMap<V>::Entry>
BlockingVerticalMap<V>::GetRestrictedSupersetEntries(const Vertical &vertical, const Vertical &exclusion) const {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetRestrictedSupersetEntries(vertical, exclusion);
}

template<class V>
bool BlockingVerticalMap<V>::RemoveSupersetEntries(const Vertical &key) {
    std::scoped_lock write_lock(read_write_mutex_);
    return VerticalMap<V>::RemoveSupersetEntries(key);
}

template<class V>
bool BlockingVerticalMap<V>::RemoveSubsetEntries(const Vertical &key) {
    std::scoped_lock write_lock(read_write_mutex_);
    return VerticalMap<V>::RemoveSubsetEntries(key);
}

template<class V>
void BlockingVerticalMap<V>::Shrink(double factor, const std::function<bool(Entry, Entry)> &compare,
                                    const std::function<bool(Entry)> &can_remove,
                                    ProfilingContext::ObjectToCache cache_object) {
    std::scoped_lock write_lock(read_write_mutex_);
  VerticalMap<V>::Shrink(factor, compare, can_remove, cache_object);
}

template<class V>
void BlockingVerticalMap<V>::Shrink(std::unordered_map<Vertical, unsigned int> &usage_counter,
                                    const std::function<bool(Entry)> &can_remove) {
    std::scoped_lock write_lock(read_write_mutex_);
  VerticalMap<V>::Shrink(usage_counter, can_remove);
}

template<class V>
long long BlockingVerticalMap<V>::GetShrinkInvocations() {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetShrinkInvocations();
}

template<class V>
long long BlockingVerticalMap<V>::GetTimeSpentOnShrinking() {
    std::shared_lock read_lock(read_write_mutex_);
    return VerticalMap<V>::GetTimeSpentOnShrinking();
}

template class BlockingVerticalMap<PositionListIndex>;
template class BlockingVerticalMap<AgreeSetSample>;
template class BlockingVerticalMap<DependencyCandidate>;
template class BlockingVerticalMap<VerticalInfo>;
template class BlockingVerticalMap<Vertical>;

} // namespace util

