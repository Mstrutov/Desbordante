#include "AgreeSetFactory.h"

#include <unordered_set>

#include "IdentifierSet.h"

#ifndef NDEBUG
#define AGREESET_DEBUG
#endif

#ifdef AGREESET_DEBUG
#define DEBUG_AGREESET(fmt, ...) \
            do { fprintf(stderr, fmt, ##__VA_ARGS__); } while (0)
#else
#define DEBUG_AGREESET(fmt, ...)
#endif

using std::set, std::vector;

template<AgreeSetsGenMethod method>
set<AgreeSet> AgreeSetFactory::genAgreeSets() const {
    auto start_time = std::chrono::system_clock::now();
    std::string method_str;
    set<AgreeSet> agree_sets;

    if constexpr (method == AgreeSetsGenMethod::kUsingVectorOfIDSets) {
        method_str = "`kUsingVectorOfIDSets`";
        vector<IdentifierSet> identifier_sets;
        set<vector<int>> const max_representation = genPLIMaxRepresentation();

        auto start_time = std::chrono::system_clock::now();

        // compute identifier sets
        // identifier_sets is vector
        std::unordered_set<int> cache;
        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                if (cache.find(*p) != cache.end())
                    continue;
                cache.insert(*p);
                identifier_sets.emplace_back(IdentifierSet(relation_, *p));
            }
        }

        auto elapsed_mills_to_gen_id_sets =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
        std::cout << "TIME TO IDENTIFIER SETS GENERATION: "
                  << elapsed_mills_to_gen_id_sets.count() << '\n';

        DEBUG_AGREESET("Identifier sets:\n");
        for (auto const& id_set : identifier_sets) {
            DEBUG_AGREESET("%s\n", id_set.toString().c_str());
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
        set<vector<int>> const max_representation = genPLIMaxRepresentation();

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
        std::cout << "TIME TO IDENTIFIER SETS GENERATION: "
                  << elapsed_mills_to_gen_id_sets.count() << '\n';


        DEBUG_AGREESET("Identifier sets:\n");
        for (auto const& [index, id_set] : identifier_sets) {
            DEBUG_AGREESET("%s\n", id_set.toString().c_str());
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
        set<vector<int>> const max_representation = genPLIMaxRepresentation();

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
    std::cout << "TIME TO AGREE SETS GENERATION WITH METHOD "
              << method_str << ": "
              << elapsed_mills_to_gen_agree_sets.count() << '\n';

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

set<vector<int>> AgreeSetFactory::genPLIMaxRepresentation() const {
    auto start_time = std::chrono::system_clock::now();
    vector<ColumnData> const& columns_data = relation_->getColumnData();
    auto not_empty_pli = std::find_if(columns_data.begin(), columns_data.end(),
                                      [](ColumnData const& c) {
        return c.getPositionListIndex()->getSize() != 0;
    });

    if (not_empty_pli == columns_data.end()) {
        return {};
    }

    //inefficient, metanome uses HashSet(std::unordered_set, need hash for vector, boost?)
    set<vector<int>> max_representation(not_empty_pli->getPositionListIndex()->getIndex().begin(),
                                        not_empty_pli->getPositionListIndex()->getIndex().end());

    for (auto p = columns_data.begin(); p != columns_data.end(); ++p) {
        if (p == not_empty_pli) //already examined
            continue;

        PositionListIndex const* pli = p->getPositionListIndex();
        if (pli->getSize() != 0)
            calculateSupersets(max_representation, pli->getIndex());
    }

    auto elapsed_mills_to_gen_max_representation =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );
    std::cout << "TIME TO MAX REPRESENTATION GENERATION: "
              << elapsed_mills_to_gen_max_representation.count() << '\n';

    return max_representation;
}

void AgreeSetFactory::calculateSupersets(set<vector<int>>& max_representation,
                                         std::deque<vector<int>> partition) const {
    set<vector<int>> to_add;
    set<vector<int>> to_delete;
    auto erase_from_partition = partition.end();
    for (auto const& max_set : max_representation) {
        for (auto p = partition.begin(); p != partition.end(); ++p) {
            if (max_set.size() >= p->size() &&
                std::includes(max_set.begin(), max_set.end(), p->begin(), p->end())) {
                to_add.erase(*p);
                erase_from_partition = p;
                break;
            }
            if (p->size() >= max_set.size() &&
                std::includes(p->begin(), p->end(), max_set.begin(), max_set.end())) {
                to_delete.insert(max_set);
            }
            to_add.insert(*p);
        }

        if (erase_from_partition != partition.end()) {
            partition.erase(erase_from_partition);
            erase_from_partition = partition.end();
        }
    }

    for (auto& cluster : to_add)
        max_representation.insert(std::move(cluster));
    for (auto& cluster : to_delete)
        max_representation.erase(cluster);
}

template
set<AgreeSet> AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingVectorOfIDSets>() const;
template
set<AgreeSet> AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingMapOfIDSets>() const;
template
set<AgreeSet> AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingGetAgreeSet>() const;
template
set<AgreeSet> AgreeSetFactory::genAgreeSets<AgreeSetsGenMethod::kUsingMCAndGetAgreeSet>() const;
