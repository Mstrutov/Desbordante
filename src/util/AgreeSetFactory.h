#pragma once

#include <set>
#include <deque>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include <boost/functional/hash.hpp>

#include "Vertical.h"
#include "ColumnLayoutRelationData.h"
#include "custom/CustomHashes.h"
#include "FDAlgorithm.h"

namespace util {

using AgreeSet = Vertical;

enum class AgreeSetsGenMethod {
    kUsingVectorOfIDSets = 0, /*< Generates agree sets using identifier sets.
                               *  Algorithm works as follows:
                               *  1. Generates maximal representation
                               *     (check out the `kUsingMCAndGetAgreeSet` description).
                               *  2. Fills vector<IdentifierSet> with identifier sets by
                               *     iterating over each cluster in max representation.
                               *     In order to avoid adding identifier set of the same tuple
                               *     twice (the same tuple can be in different clusters), the
                               *     set of ids of the already added tuples is used.
                               *  3. Iterates over all pairs of identifier sets from vector.
                               *  4. Gets agree set for current pair by intersecting.
                               */
    kUsingMapOfIDSets       , /*< Metanome approach.
                               *  Generates agree sets using identifier sets.
                               *  Algorithm works as follows:
                               *  1. Generates maximal representation
                               *     (check out the `kUsingMCAndGetAgreeSet` description).
                               *  2. Fills map<int, IdentifierSet> with <tuple index, tuple idset>
                               *     pairs by iterating over each cluster in max representation.
                               *  3. Iterates over all pairs of tuples from each cluster
                               *     of maximal representation.
                               *  4. Gets agree set for current pair of tuples by intersecting
                               *     their identifier sets.
                               */
    kUsingGetAgreeSet       , /*< The most naive (so the slowest) way to generate agree sets.
                               *  Generates agree set for all pairs of tuples that
                               *  can form agree set. Tuples can form agree set if
                               *  they have the same value in at least one attribute.
                               *  Such tuples are in the same cluster of pli.
                               *  So algorithm works as follows:
                               *  1. Iterates over all table attributes.
                               *  2. Iterates over all pairs of tuples from each cluster
                               *     of current attribute pli.
                               *  3. Gets agree set for current pair using getAgreeSet.
                               */
    kUsingMCAndGetAgreeSet  , /*< In `kUsingGetAgreeSet` method, the same pair of tuples
                               *  can be processed (passed to getAgreeSet) multiple times.
                               *  This method similar to `kUsingGetAgreeSet`,
                               *  but avoids the above problem.
                               *  It reduces the number of redundant tuple comparisons by
                               *  transforming partitions into its maximal representation.
                               *  And processes tuple pairs from maximal representation clusters.
                               *  Algorithm of maximal representation generation was taken from
                               *  Metanome. For more information about maximal representation
                               *  check out http://www.vldb.org/pvldb/vol8/p1082-papenbrock.pdf
                               */
};

/* Max representation generation method */
enum class MCGenMethod {
    kUsingHandleEqvClass = 0, /*< Naive way to generate maximal representation.
                               *  It sorts all equivalence classes (represented by vector<int>)
                               *  from all partitions in ascending order. One equivalence class
                               *  is less than another iff its size is less or sizes are equal
                               *  and one is lexicographically smaller. After equivalence classes
                               *  are sorted by insertion into set, algorithm proceeds to maximize
                               *  partitions.
                               *  It maintains map<eqv_class_size, Set<eqv_class>>:
                               *  1. Adds set of equivalence classes with minimum size to the map.
                               *  2. Iterates over sorted equivalence classes with non minimal size,
                               *     adds each class to the map and checks for every subset of
                               *     current class if map contains it. If it is, deletes subset
                               *     from the map.
                               *  After all equivalence classes have been processed, map stores the
                               *  maximum representation.
                               */
    kUsingHandlePartition,     /*< Sorts all equivalence classes from all partitions in descending
                                *  order using std::set. Definition of comparison for equivalence
                                *  classes is same as for kUsingHandleEqvClass. Iterates over
                                *  sorted_eqv_classes and maintains index as std::map<int, set<int>>.
                                *  A more detailed description of index is in the implementation.
                                *  'handlePartition':
                                *  For current eqv_class checks if it has superset in the index using
                                *  isSubset(), if not, adds it to the index and to
                                *  max_representation.
                                */
    kUsingCalculateSupersets,  /*< Fills max_representation with equivalence classes from first not
                                *  empty partition. Then iterates over the remaining partitions
                                *  and edits max_representation on the fly using calculateSupersets:
                                *  1. Iterates over sets from max_representation.
                                *  2. If current set is a subset of equivalence class from the
                                *     current partition, then deletes (after max_representation
                                *     has been fully processed) it from the max_representation.
                                *     And adds (also delayed) appropriate equivalence class to
                                *     max_representation.
                                */
    kParallel                  /*< Algorithm is the same as in kUsingHandlePartition method.
                                *  Uses thread pool of config_.threads_num threads to perform
                                *  'handlePartition' part on equivalence classes of the same size.
                                */
};

class AgreeSetFactory {
public:
    struct Configuration {
        AgreeSetsGenMethod as_gen_method = AgreeSetsGenMethod::kUsingMapOfIDSets;
        MCGenMethod mc_gen_method = MCGenMethod::kUsingCalculateSupersets;
        ushort threads_num = 1;

        /* Not using default keyword because of gcc bug:
         * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=88165
         */
        Configuration() noexcept {}
        explicit Configuration(AgreeSetsGenMethod as_gen_m,
                               MCGenMethod mc_gen_m,
                               ushort threads_num) noexcept
            : as_gen_method(as_gen_m), mc_gen_method(mc_gen_m), threads_num(threads_num) {}
        explicit Configuration(AgreeSetsGenMethod as_gen_m) noexcept
            : as_gen_method(as_gen_m) {}
        explicit Configuration(MCGenMethod mc_gen_m) noexcept
            : mc_gen_method(mc_gen_m) {}
        explicit Configuration(ushort threads_num) noexcept
            : threads_num(threads_num) {}
    };
    using SetOfVectors = std::unordered_set<std::vector<int>,
                                            boost::hash<std::vector<int>>>;
    using SetOfAgreeSets = std::unordered_set<AgreeSet>;

    explicit AgreeSetFactory(ColumnLayoutRelationData const* const rel,
                             Configuration const& c = Configuration(),
                             FDAlgorithm* algo = nullptr)
        : relation_(rel), config_(c), algo_(algo) {}

    ColumnLayoutRelationData const* getRelation() const { return relation_; }
    void SetConfiguration(Configuration const& c) { config_ = c; }

    // Computes all agree sets of `relation_` using specified method
    SetOfAgreeSets genAgreeSets() const;

    SetOfVectors genPLIMaxRepresentation() const;

    AgreeSet getAgreeSet(int const tuple1_index, int const tuple2_index) const;
private:
    /* Implementations of generation agree sets algorithms */
    SetOfAgreeSets genASUsingVectorOfIDSets() const;
    SetOfAgreeSets genASUsingMapOfIDSets() const;
    SetOfAgreeSets genASUsingGetAgreeSets() const;
    SetOfAgreeSets genASUsingMCAndGetAgreeSets() const;

    /* Implementations of generation MC algorithms */
    SetOfVectors genMCUsingHandleEqvClass() const;
    SetOfVectors genMCUsingHandlePartition() const;
    SetOfVectors genMCUsingCalculateSupersets() const;
    SetOfVectors genMCParallel() const;

    void calculateSupersets(SetOfVectors& max_representation,
                            std::deque<std::vector<int>> const& partition) const;
    /* From Metanome: `handleList`.
     * Extremely slow for anything big eqv_class,
     * I think it is not usable at all
     */
    void handleEqvClass(std::vector<int>& eqv_class,
                        std::unordered_map<size_t, SetOfVectors>& max_sets,
                        bool const first_step) const;
    /* From Metanome.
     * Checks if index 'contains' eqvivalence class which is superset of eqv_class.
     */
    bool isSubset(std::vector<int> const& eqv_class,
                  std::unordered_map<int, std::unordered_set<size_t>> const& index) const;
    using VectorComp = std::function<bool (std::vector<int> const&, std::vector<int> const&)>;
    std::set<std::vector<int>, VectorComp> genSortedEqvClasses(VectorComp comp) const;

    void addProgress(double const val) const noexcept {
        if (algo_ != nullptr) {
            algo_->addProgress(val);
        }
    }


    ColumnLayoutRelationData const* const relation_;

    Configuration config_;
    FDAlgorithm* algo_;
};

} // namespace util

