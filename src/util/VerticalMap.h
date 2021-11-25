#pragma once
#include <functional>
#include <memory>
#include <shared_mutex>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "custom/CustomHashes.h"
#include "ProfilingContext.h"

namespace util {

//difficulties with const methods

//Value: PLI, AgreeSetSample, VerticalInfo, DependencyCandidate, Vertical <- all of these are shared_ptrs?
//Use template specialization for shared_ptr<Value>?
template <class Value>
class VerticalMap {
protected:
    using bitset = boost::dynamic_bitset<>;
    //typename std::shared_ptr<Value> shared_ptr<Value>;

    // Each node corresponds to a bit in a bitset. Each node also has a vector of the possible consequent set bits.
    class SetTrie {
    private:
        size_t offset_;
        size_t dimension_;
        std::vector<std::unique_ptr<SetTrie>> subtries_; //unique_ptr?
        std::shared_ptr<Value> value_;

        bool isEmpty() const;

    public:
        explicit SetTrie(size_t dimension) : SetTrie(0, dimension) {}
        SetTrie(size_t offset, size_t dimension) : offset_(offset), dimension_(dimension) {}

        // Sets given key to a given value
        // Returns the old value with ownership
        // Not a const method as SetTrie gets changed
        std::shared_ptr<Value> associate(bitset const& key, size_t nextBit, std::shared_ptr<Value> value);

        // Returns a pointer to the value mapped by the given key
        std::shared_ptr<Value const> get(bitset const& key, size_t nextBit) const;

        // Erases an entry with the given key
        // Returns the old value with ownership
        std::shared_ptr<Value> remove(bitset const& key, size_t nextBit);

        // Gets the subtrie with the given index. If such a subtrie does not exist, creates one
        // Not a const method as a SetTrie may be created
        SetTrie* getOrCreateSubTrie(size_t index);

        // Gets the subtrie with the given index
        // Non-const version of the method -- needed in getOrCreateSubtrie to manipulate SetTrie*
        SetTrie* getSubtrie(size_t index);

        // Gets the subtrie with the given index
        // Const version of the method
        SetTrie const* getSubtrie(size_t index) const;

        // Calls collector on every trie that is a subset of the given subsetKey
        bool collectSubsetKeys(bitset const& key, size_t nextBit, bitset& subsetKey,
                               std::function<bool(bitset const&, std::shared_ptr<Value const>)> const& collector) const;

        // Calls collector on every trie that is a superset of the given subsetKey
        bool collectSupersetKeys(bitset const& key, size_t nextBit, bitset& supersetKey,
                                 std::function<bool(bitset const&, std::shared_ptr<Value const>)> const& collector) const;

        // Calls collector on every trie that is a superset of the given subsetKey with no bits from the blacklist
        bool collectRestrictedSupersetKeys(bitset const& key, bitset const& blacklist, size_t nextBit,
                                           bitset& supersetKey,
                                           std::function<void(bitset const&, std::shared_ptr<Value const>)> const& collector) const;

        // Calls collector on every entry
        void traverseEntries(bitset& subsetKey, std::function<void(bitset const&, std::shared_ptr<Value const>)> collector) const;
    };

    RelationalSchema const* relation_;
    size_t size_ = 0;
    long long shrinkInvocations_ = 0;
    long long timeSpentOnShrinking_ = 0;
    SetTrie setTrie_;

    unsigned int removeFromUsageCounter (std::unordered_map<Vertical, unsigned int>& usageCounter, const Vertical& key);
public:
    using Entry = std::pair<Vertical, std::shared_ptr<Value const>>;
    explicit VerticalMap(RelationalSchema const* relation) :
        relation_(relation), setTrie_(relation->getNumColumns()) {}
    virtual size_t getSize() const { return size_; }
    virtual bool isEmpty() const { return size_ == 0; }

    // basic get-check-insert-remove operations
    virtual std::shared_ptr<Value const> get(Vertical const &key) const;
    virtual std::shared_ptr<Value const> get(bitset const &key) const;
    virtual bool containsKey(Vertical const& key) const { return get(key) != nullptr; }
    virtual std::shared_ptr<Value> put(Vertical const& key, std::shared_ptr<Value> value);
    virtual std::shared_ptr<Value> remove(Vertical const& key);
    virtual std::shared_ptr<Value> remove(bitset const& key);

    // non-const version of get() for костыль purposes
    virtual std::shared_ptr<Value> get(Vertical const& key);

    // get all keys/values/entries for traversing
    virtual std::unordered_set<Vertical> keySet();
    virtual std::vector<std::shared_ptr<Value const>> values();
    virtual std::unordered_set<Entry> entrySet();

    // get specific entries for traversing
    virtual std::vector<Vertical> getSubsetKeys(Vertical const& vertical) const;
    virtual std::vector<Entry> getSubsetEntries(Vertical const& vertical) const;
    virtual Entry getAnySubsetEntry(Vertical const& vertical) const;
    virtual Entry getAnySubsetEntry(Vertical const& vertical, std::function<bool(Vertical const*, std::shared_ptr<Value const>)> const& condition) const;
    virtual std::vector<Entry> getSupersetEntries(Vertical const& vertical) const;
    virtual Entry getAnySupersetEntry(Vertical const& vertical) const;
    virtual Entry getAnySupersetEntry(Vertical const& vertical, std::function<bool(Vertical const*, std::shared_ptr<Value const>)> condition) const;
    virtual std::vector<Entry> getRestrictedSupersetEntries(Vertical const& vertical, Vertical const& exclusion) const;
    virtual bool removeSupersetEntries(Vertical const& key);
    virtual bool removeSubsetEntries(Vertical const& key);

    /* methods to shrink the map by deleting removable entries
     * !!! Untested yet - use carefully
     * */
    virtual void shrink(double factor, std::function<bool(Entry, Entry)> const& compare,
                        std::function<bool (Entry)> const& canRemove, ProfilingContext::ObjectToCache cacheObject);
    virtual void shrink(std::unordered_map<Vertical, unsigned int>& usageCounter, std::function<bool (Entry)> const& canRemove);

    virtual long long getShrinkInvocations() { return shrinkInvocations_; }
    virtual long long getTimeSpentOnShrinking() { return timeSpentOnShrinking_; }

    virtual ~VerticalMap() = default;
};

/*
 * A version of VerticalMap for parallel processing. Uses reader-writer mutex for blocking.
 * */
template <class V>
class BlockingVerticalMap : public VerticalMap<V> {
private:
    // mutable нужен, чтобы const методы могли его лочить (т.к. меняется состояние мьютекса)
    mutable std::shared_mutex readWriteMutex_;
    // shared_mutex, scoped_lock, shared_lock
public:
    using typename VerticalMap<V>::Entry;
    using typename VerticalMap<V>::bitset;

    explicit BlockingVerticalMap(RelationalSchema const* relation) : VerticalMap<V>(relation) {}
    virtual size_t getSize() const override;
    virtual bool isEmpty() const override;

    virtual std::shared_ptr<V const> get(Vertical const &key) const override;
    virtual std::shared_ptr<V const> get(bitset const &key) const override;
    virtual bool containsKey(Vertical const& key) const override;
    virtual std::shared_ptr<V> put(Vertical const& key, std::shared_ptr<V> value) override;
    virtual std::shared_ptr<V> remove(Vertical const& key) override;
    virtual std::shared_ptr<V> remove(bitset const& key) override;

    virtual std::shared_ptr<V> get(Vertical const& key) override;

    virtual std::unordered_set<Vertical> keySet() override;
    virtual std::vector<std::shared_ptr<V const>> values() override;
    virtual std::unordered_set<Entry> entrySet() override;

    virtual std::vector<Vertical> getSubsetKeys(Vertical const& vertical) const override;
    virtual std::vector<Entry> getSubsetEntries(Vertical const& vertical) const override;
    virtual Entry getAnySubsetEntry(Vertical const& vertical) const override;
    virtual Entry getAnySubsetEntry(
            Vertical const& vertical, std::function<bool(Vertical const*, std::shared_ptr<V const>)> const& condition) const override;
    virtual std::vector<Entry> getSupersetEntries(Vertical const& vertical) const override;
    virtual Entry getAnySupersetEntry(Vertical const& vertical) const override;
    virtual Entry getAnySupersetEntry(
            Vertical const& vertical, std::function<bool(Vertical const*, std::shared_ptr<V const>)> condition) const override;
    virtual std::vector<Entry> getRestrictedSupersetEntries(
            Vertical const& vertical, Vertical const& exclusion) const override;
    virtual bool removeSupersetEntries(Vertical const& key) override;
    virtual bool removeSubsetEntries(Vertical const& key) override;

    virtual void shrink(double factor, std::function<bool(Entry, Entry)> const& compare,
                        std::function<bool (Entry)> const& canRemove, ProfilingContext::ObjectToCache cacheObject)
                        override;
    virtual void shrink(std::unordered_map<Vertical, unsigned int>& usageCounter,
                        std::function<bool (Entry)> const& canRemove) override;

    virtual long long getShrinkInvocations() override;
    virtual long long getTimeSpentOnShrinking() override;

    virtual ~BlockingVerticalMap() = default;
};

} // namespace util

