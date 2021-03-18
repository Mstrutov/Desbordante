#pragma once
#include <functional>
#include <memory>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "custom/CustomHashes.h"
#include "ProfilingContext.h"



//difficulties with const methods

//Value: PLI, AgreeSetSample, VerticalInfo, DependencyCandidate, Vertical <- all of these are shared_ptrs?
//Use template specialization for shared_ptr<Value>?
template <class Value>
class VerticalMap {
private:
    using bitset = boost::dynamic_bitset<>;

    // Each node corresponds to a bit in a bitset. Each node also has a vector of the possible consequent set bits.
    class SetTrie {
    private:
        size_t offset_;
        size_t dimension_;
        std::vector<std::unique_ptr<SetTrie>> subtries_; //unique_ptr?
        std::unique_ptr<Value> value_;

        bool isEmpty() const;

    public:
        explicit SetTrie(size_t dimension) : SetTrie(0, dimension) {}
        SetTrie(size_t offset, size_t dimension) : offset_(offset), dimension_(dimension) {}

        // Sets given key to a given value
        // Returns the old value with ownership
        // Not a const method as SetTrie gets changed
        std::unique_ptr<Value> associate(bitset const& key, size_t nextBit, std::unique_ptr<Value> value);

        // Returns a pointer to the value mapped by the given key
        Value const* get(bitset const& key, size_t nextBit) const;

        // Erases an entry with the given key
        // Returns the old value with ownership
        std::unique_ptr<Value> remove(bitset const& key, size_t nextBit);

        // Gets the subtrie with the given index. If such a subtrie does not exist, creates one
        // Not a const method as a SetTrie may be created
        SetTrie* getOrCreateSubTrie(size_t index);

        // Gets the subtrie with the given index
        // Non-const version of the method
        SetTrie* getSubtrie(size_t index);

        // Gets the subtrie with the given index
        // Const version of the method
        SetTrie const* getSubtrie(size_t index) const;

        // Calls collector on every trie that is a subset of the given subsetKey
        bool collectSubsetKeys(bitset const& key, size_t nextBit, bitset& subsetKey,
                               std::function<bool(bitset&&, Value const*)> const& collector) const;

        // Calls collector on every trie that is a superset of the given subsetKey
        bool collectSupersetKeys(bitset const& key, size_t nextBit, bitset& supersetKey,
                                 std::function<bool(bitset&&, Value const*)> const& collector) const;

        // Calls collector on every trie that is a superset of the given subsetKey with no bits from the blacklist
        bool collectRestrictedSupersetKeys(bitset const& key, bitset const& blacklist, size_t nextBit,
                                           bitset& supersetKey,
                                           std::function<void(bitset&&, Value const*)> const& collector) const;

        // Calls collector on every entry
        void traverseEntries(bitset& subsetKey, std::function<void(bitset const&, Value const*)> collector) const;
    };

    RelationalSchema const* relation_;
    size_t size_ = 0;
    long long shrinkInvocations_ = 0;
    long long timeSpentOnShrinking_ = 0;
    SetTrie setTrie_;

    unsigned int removeFromUsageCounter (std::unordered_map<Vertical, unsigned int>& usageCounter, const Vertical& key);
public:
    using Entry = std::pair<Vertical, Value const*>;
    explicit VerticalMap(RelationalSchema const* relation) :
        relation_(relation), setTrie_(relation->getNumColumns()) {}
    size_t getSize() const { return size_; }
    bool isEmpty() const { return size_ == 0; }

    // basic get-check-insert-remove operations
    Value const* get(Vertical const &key) const;
    Value const* get(bitset const &key) const;
    bool containsKey(Vertical const& key) { return get(key) != nullptr; }
    std::unique_ptr<Value> put(Vertical const& key, std::unique_ptr<Value> value);
    std::unique_ptr<Value> remove(Vertical const& key);
    std::unique_ptr<Value> remove(bitset const& key);

    // non-const version of get() for костыль purposes
    Value* get(Vertical const& key);

    // get all keys/values/entries for traversing
    std::unordered_set<Vertical> keySet();
    std::vector<Value const*> values();
    std::unordered_set<Entry> entrySet();

    // get specific entries for traversing
    std::vector<Vertical> getSubsetKeys(Vertical const& vertical) const;
    std::vector<Entry> getSubsetEntries(Vertical const& vertical) const;
    Entry getAnySubsetEntry(Vertical const& vertical) const;
    Entry getAnySubsetEntry(Vertical const& vertical, std::function<bool(Vertical const*, Value const*)> const& condition) const;
    std::vector<Entry> getSupersetEntries(Vertical const& vertical) const;
    Entry getAnySupersetEntry(Vertical const& vertical) const;
    Entry getAnySupersetEntry(Vertical const& vertical, std::function<bool(Vertical const*, Value const*)> condition) const;
    std::vector<Entry> getRestrictedSupersetEntries(Vertical const& vertical, Vertical const& exclusion) const;
    bool removeSupersetEntries(Vertical const& key);
    bool removeSubsetEntries(Vertical const& key);

    // methods to shrink the map by deleting removable entries
    void shrink(double factor, std::function<bool(Entry, Entry)> const& compare,
                std::function<bool (Entry)> const& canRemove, ProfilingContext::ObjectToCache cacheObject);
    void shrink(std::unordered_map<Vertical, unsigned int>& usageCounter, std::function<bool (Entry)> const& canRemove);

    long long getShrinkInvocations() { return shrinkInvocations_; }
    long long getTimeSpentOnShrinking() { return timeSpentOnShrinking_; }
};

