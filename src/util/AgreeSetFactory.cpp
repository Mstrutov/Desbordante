#include "AgreeSetFactory.h"

#include <unordered_set>

#include "IdentifierSet.h"
#include "logging/easylogging++.h"

using std::set, std::vector, std::unordered_set;

template<AgreeSetsGenMethod method>
AgreeSetFactory::SetOfAgreeSets AgreeSetFactory::genAgreeSets() const {
    auto start_time = std::chrono::system_clock::now();
    std::string method_str;
    SetOfAgreeSets agree_sets;

    if constexpr (method == AgreeSetsGenMethod::kUsingVectorOfIDSets) {
        method_str = "`kUsingVectorOfIDSets`";
        vector<IdentifierSet> identifier_sets;
        SetOfVectors const max_representation = genPLIMaxRepresentation();

        auto start_time = std::chrono::system_clock::now();

        // compute identifier sets
        // identifier_sets is vector
        std::unordered_set<int> cache;
        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                if (!cache.insert(*p).second) {
                    continue;
                }
                identifier_sets.emplace_back(IdentifierSet(relation_, *p));
            }
        }

        auto elapsed_mills_to_gen_id_sets =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
        LOG(INFO) << "TIME TO IDENTIFIER SETS GENERATION: "
                  << elapsed_mills_to_gen_id_sets.count();

        LOG(DEBUG) << "Identifier sets:";
        for (auto const& id_set : identifier_sets) {
            LOG(DEBUG) << id_set.toString();
        }

        // compute agree sets using identifier sets
        // using vector of identifier sets
        if (!identifier_sets.empty()) {
            auto back_it = std::prev(identifier_sets.end());
            for (auto p = identifier_sets.begin(); p != back_it; ++p) {
                for (auto q = std::next(p); q != identifier_sets.end(); ++q) {
                    agree_sets.insert(p->intersect(*q));
                }
            }
        }
    } else if constexpr (method == AgreeSetsGenMethod::kUsingMapOfIDSets) {
        method_str = "`kUsingMapOfIDSets`";
        std::unordered_map<int, IdentifierSet> identifier_sets;
        SetOfVectors const max_representation = genPLIMaxRepresentation();

        auto start_time = std::chrono::system_clock::now();

        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                identifier_sets.emplace(*p, IdentifierSet(relation_, *p));
            }
        }

        auto elapsed_mills_to_gen_id_sets =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
        LOG(INFO) << "TIME TO IDENTIFIER SETS GENERATION: "
                  << elapsed_mills_to_gen_id_sets.count();


        LOG(DEBUG) << "Identifier sets:";
        for (auto const& [index, id_set] : identifier_sets) {
            LOG(DEBUG) << id_set.toString();
        }


        // compute agree sets using identifier sets
        // metanome approach (using map of identifier sets)
        for (auto const &cluster : max_representation) {
            auto back_it = std::prev(cluster.end());
            for (auto p = cluster.begin(); p != back_it; ++p) {
                for (auto q = std::next(p); q != cluster.end(); ++q) {
                    IdentifierSet const& id_set1 = identifier_sets.at(*p);
                    IdentifierSet const& id_set2 = identifier_sets.at(*q);
                    agree_sets.insert(id_set1.intersect(id_set2));
                }
            }
        }
    } else if constexpr (method == AgreeSetsGenMethod::kUsingMCAndGetAgreeSet) {
        method_str = "`kUsingMCAndGetAgreeSet`";
        SetOfVectors const max_representation = genPLIMaxRepresentation();

        // Compute agree sets from maximal representation using getAgreeSet()
        // ~3300 ms on CIPublicHighway700 (Debug build), ~250 ms (Release)
        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                for (auto q = std::next(p); q != cluster.end(); ++q) {
                    agree_sets.insert(getAgreeSet(*p, *q));
                }
            }
        }
    } else if constexpr (method == AgreeSetsGenMethod::kUsingGetAgreeSet) {
        method_str = "`kUsingGetAgreeSet`";
        vector<ColumnData> const& columns_data = relation_->getColumnData();

        // Compute agree sets from stripped partitions (simplest method by Wyss)
        // ~40436 ms on CIPublicHighway700 (Debug build)
        for (ColumnData const& column_data : columns_data) {
            PositionListIndex const* const pli = column_data.getPositionListIndex();
            for (vector<int> const& cluster : pli->getIndex()) {
                for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                    for (auto q = std::next(p); q != cluster.end(); ++q) {
                        agree_sets.insert(getAgreeSet(*p, *q));
                    }
                }
            }
        }
    }

    // metanome kostil, doesn't work properly in general
    agree_sets.insert(*relation_->getSchema()->emptyVertical);

    auto elapsed_mills_to_gen_agree_sets =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
    LOG(INFO) << "TIME TO GENERATE AGREE SETS WITH METHOD "
              << method_str << ": "
              << elapsed_mills_to_gen_agree_sets.count();

    return agree_sets;
}

AgreeSet AgreeSetFactory::getAgreeSet(int const tuple1_index,
                                      int const tuple2_index) const {
    std::vector<int> const tuple1 = relation_->getTuple(tuple1_index);
    std::vector<int> const tuple2 = relation_->getTuple(tuple2_index);
    boost::dynamic_bitset<> agree_set_indices(relation_->getNumColumns());

    for (size_t i = 0; i < agree_set_indices.size(); ++i) {
        if (tuple1[i] != 0 && tuple1[i] == tuple2[i]) {
            agree_set_indices.set(i);
        }
    }

    return relation_->getSchema()->getVertical(agree_set_indices);
}

/* It seems very cumbersome (so is the genAgreeSet method),
 * maybe need to revise the algorithm for choosing the desired method.
 */
template<MCGenMethod method>
AgreeSetFactory::SetOfVectors AgreeSetFactory::genPLIMaxRepresentation() const {
    auto start_time = std::chrono::system_clock::now();
    vector<ColumnData> const& columns_data = relation_->getColumnData();
    std::string method_str;
    SetOfVectors max_representation;

    if constexpr (method == MCGenMethod::kUsingCalculateSupersets) {
        method_str = "`kUsingCalculateSupersets`";
        auto not_empty_pli =
            std::find_if(columns_data.begin(), columns_data.end(),
                         [](ColumnData const& c) {
                             return c.getPositionListIndex()->getSize() != 0;
                         }
            );

        if (not_empty_pli == columns_data.end()) {
            goto OUT;
        }

        max_representation.insert(not_empty_pli->getPositionListIndex()->getIndex().begin(),
                                  not_empty_pli->getPositionListIndex()->getIndex().end());

        for (auto p = std::next(not_empty_pli); p != columns_data.end(); ++p) {
            PositionListIndex const* pli = p->getPositionListIndex();
            if (pli->getSize() != 0) {
                calculateSupersets(max_representation, pli->getIndex());
            }
        }
    } else if constexpr (method == MCGenMethod::kUsingHandleEqvClass) {
        method_str = "`kUsingHandleEqvClass`";
        // set of all equivalence classes of all paritions
        auto comp = [](vector<int> const& lhs, vector<int> const& rhs) {
            if (lhs.size() != rhs.size()) {
                return lhs.size() < rhs.size();
            }
            return std::lexicographical_compare(lhs.begin(), lhs.end(),
                                                rhs.begin(), rhs.end());
        };
        set<vector<int>, decltype(comp)> sorted_eqv_classes(comp);

        // Fill sorted_partitions
        for (ColumnData const& data : columns_data) {
            std::deque<vector<int>> const& index = data.getPositionListIndex()->getIndex();
            sorted_eqv_classes.insert(index.begin(), index.end());
        }

        if (sorted_eqv_classes.empty()) {
            goto OUT;
        }


        // Maximize partitions
        size_t const min_size = sorted_eqv_classes.begin()->size();
        std::unordered_map<size_t, SetOfVectors> max_sets;
        SetOfVectors min_set;

        auto const first_not_min
            = std::find_if_not(sorted_eqv_classes.begin(), sorted_eqv_classes.end(),
                               [min_size](auto const& p) { return p.size() == min_size; });

        min_set.insert(sorted_eqv_classes.begin(), first_not_min);
        max_sets.emplace(min_size, std::move(min_set));

        for (auto it = first_not_min; it != sorted_eqv_classes.end();) {
            // So that the eqv_class can be modified
            vector<int> eqv_class = sorted_eqv_classes.extract(it++).value();
            handleEqvClass(eqv_class, max_sets, true);
        }


        // Metanome `mergeResult`
        max_representation.reserve(max_sets.size());
        for (auto it = max_sets.begin(); it != max_sets.end();) {
            SetOfVectors set = max_sets.extract(it++).mapped();
            max_representation.insert(std::make_move_iterator(set.begin()),
                                      std::make_move_iterator(set.end()));
        }
        assert(max_sets.empty());
    }

OUT:
    auto elapsed_mills_to_gen_max_representation =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );
    LOG(INFO) << "TIME TO GENERATE MAX REPRESENTATION WITH METHOD "
              << method_str << ": "
              << elapsed_mills_to_gen_max_representation.count();

    return max_representation;
}

void AgreeSetFactory::handleEqvClass(vector<int>& eqv_class,
                                     std::unordered_map<size_t, SetOfVectors>& max_sets,
                                     bool const first_step) const {
    for (auto it = eqv_class.begin(); it != eqv_class.end(); ++it) {
        vector<int> copy(eqv_class.begin(), it);
        copy.insert(copy.end(), std::next(it), eqv_class.end());

        size_t const size = copy.size();
        assert(size == eqv_class.size() - 1);

        /* Need to check if an element with copy.size() key exists before accessing it
         * to avoid new SetOfVectors() creation if it doesn't.
         */
        if (max_sets.count(size) != 0 &&
            max_sets[size].find(copy) != max_sets[size].end()) {
            max_sets[size].erase(copy);
        } else {
            if (size > 2) {
                handleEqvClass(copy, max_sets, false);
            }
        }
    }

    if (first_step)
        max_sets[eqv_class.size()].insert(std::move(eqv_class));
}

void AgreeSetFactory::calculateSupersets(SetOfVectors& max_representation,
                                         std::deque<vector<int>> const& partition) const {
    SetOfVectors to_add_to_mc;
    auto hash = [beg = max_representation.begin()](SetOfVectors::const_iterator it) {
        return std::distance<SetOfVectors::const_iterator>(beg, it);
    };
    unordered_set<SetOfVectors::const_iterator, decltype(hash)> to_delete_from_mc(1, hash);
    set<std::deque<vector<int>>::const_iterator> to_exclude_from_partition;

    for (auto it = max_representation.begin(); it != max_representation.end(); ++it) {
        for (auto p = partition.begin();
             to_exclude_from_partition.size() != partition.size() && p != partition.end();
             ++p) {
            if (to_exclude_from_partition.find(p) != to_exclude_from_partition.end()) {
                continue;
            }

            if (it->size() >= p->size() &&
                std::includes(it->begin(), it->end(), p->begin(), p->end())) {
                to_add_to_mc.erase(*p);
                to_exclude_from_partition.insert(p);
                break;
            }

            if (p->size() >= it->size() &&
                std::includes(p->begin(), p->end(), it->begin(), it->end())) {
                to_delete_from_mc.insert(it);
            }

            to_add_to_mc.insert(*p);
        }
    }

    for (auto& cluster : to_add_to_mc) {
        max_representation.insert(std::move(cluster));
    }
    for (auto it : to_delete_from_mc) {
        max_representation.erase(it);
    }
}


template AgreeSetFactory::SetOfAgreeSets
AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingVectorOfIDSets>() const;
template AgreeSetFactory::SetOfAgreeSets
AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingMapOfIDSets>() const;
template AgreeSetFactory::SetOfAgreeSets
AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingGetAgreeSet>() const;
template AgreeSetFactory::SetOfAgreeSets
AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingMCAndGetAgreeSet>() const;

template AgreeSetFactory::SetOfVectors
AgreeSetFactory::genPLIMaxRepresentation<MCGenMethod::kUsingCalculateSupersets>() const;
template AgreeSetFactory::SetOfVectors
AgreeSetFactory::genPLIMaxRepresentation<MCGenMethod::kUsingHandleEqvClass>() const;

