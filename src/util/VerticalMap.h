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



//TODO: Expect linking problems. Implement function constructing specific objects to give the compiler necessary info - done at the end of .cpp

//difficulties with const methods

//Value: PLI, AgreeSetSample, VerticalInfo, DependencyCandidate, Vertical <- all of these are shared_ptrs?
//Use template specialization for shared_ptr<Value>?
template <class Value>
class VerticalMap {
private:
    //in Java nested class has implicitly a reference to a parent - is this used? - looks like no
    using bitset = boost::dynamic_bitset<>;
    class SetTrie {
    private:
        size_t offset_;
        size_t dimension_;
        std::vector<std::shared_ptr<SetTrie>> subtries_; //unique_ptr?
        Value value_;

        bool isEmpty() const;

    public:
        explicit SetTrie(size_t dimension) : SetTrie(0, dimension) {}
        SetTrie(size_t offset, size_t dimension) : offset_(offset), dimension_(dimension) {}

        Value associate(bitset const& key, size_t nextBit, Value value);
        Value get(bitset const& key, size_t nextBit);
        Value remove(bitset const& key, size_t nextBit);
        //unchecked// - TODO: shared_ptr<SetTrie>
        std::shared_ptr<SetTrie> getOrCreateSubTrie(size_t index);
        //unchecked//
        std::shared_ptr<SetTrie> getSubtrie(size_t index);
        bool collectSubsetKeys(bitset const& key, size_t nextBit, bitset& subsetKey, std::function<bool(bitset&&, Value)> const& collector);
        bool collectSupersetKeys(bitset const& key, size_t nextBit, bitset& supersetKey, std::function<bool(bitset&&, Value)> const& collector);
        bool collectRestrictedSupersetKeys(bitset const& key, bitset const& blacklist, size_t nextBit, bitset& supersetKey, std::function<void(bitset&&, Value)> const& collector);

        void traverseEntries(bitset& subsetKey, std::function<void(bitset const&, Value)> collector);
    };

    std::weak_ptr<RelationalSchema> relation_;
    size_t size_ = 0;
    long long shrinkInvocations_ = 0;
    long long timeSpentOnShrinking_ = 0;
    //std::unordered_map<Vertical, Value> map_;
    SetTrie setTrie_;

    unsigned int removeFromUsageCounter (std::unordered_map<Vertical, unsigned int>& usageCounter, Vertical key);
public:
    //key = shared_ptr<const Vertical>?? Attempt to achieve const Vertical?
    using Entry = std::pair<std::shared_ptr<Vertical>, Value>;
    explicit VerticalMap(std::shared_ptr<RelationalSchema> relation) : relation_(relation), setTrie_(relation->getNumColumns()) {}
    size_t getSize() const { return size_; }
    bool isEmpty() const { return size_ == 0; }

    Value get(Vertical const &key);
    Value get(bitset const &key);

    bool containsKey(Vertical const& key) { return get(key) != nullptr; }
    Value put(Vertical const& key, Value value);

    Value remove(Vertical const& key);
    Value remove(bitset const& key);

    std::vector<std::shared_ptr<Vertical>> getSubsetKeys(Vertical const& vertical);
    std::vector<Entry> getSubsetEntries(Vertical const& vertical);
    Entry getAnySubsetEntry(Vertical const& vertical);
    Entry getAnySubsetEntry(Vertical const& vertical, std::function<bool(Vertical*, Value)> const& condition);
    std::vector<Entry> getSupersetEntries(Vertical const& vertical);
    Entry getAnySupersetEntry(Vertical const& vertical);
    Entry getAnySupersetEntry(Vertical const& vertical, std::function<bool(Vertical*, Value)> condition);
    std::vector<Entry> getRestrictedSupersetEntries(Vertical const& vertical, Vertical const& exclusion);
    bool removeSupersetEntries(Vertical const& key);
    bool removeSubsetEntries(Vertical const& key);
    std::unordered_set<std::shared_ptr<Vertical>> keySet();
    std::vector<Value> values();
    std::unordered_set<Entry> entrySet();
    void shrink(double factor, std::function<bool(Entry, Entry)> const& compare, std::function<bool (Entry)> const& canRemove, ProfilingContext::ObjectToCache cacheObject);
    void shrink(std::unordered_map<Vertical, unsigned int>& usageCounter, std::function<bool (Entry)> const& canRemove);
    void shrink() = delete; //concurrent shrink?

    long long getShrinkInvocations() { return shrinkInvocations_; }
    long long getTimeSpentOnShrinking() { return timeSpentOnShrinking_; }
};

