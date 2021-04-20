#include "FastFDs.h"

#include <algorithm>

#include <boost/dynamic_bitset.hpp>

#include "IdentifierSet.h"

#ifndef NDEBUG
    #define FASTFDS_DEBUG
#endif

using std::vector, std::shared_ptr, std::set;

namespace std {
    template<>
    struct less<shared_ptr<Vertical>> {
        bool operator()(shared_ptr<Vertical> const& l, shared_ptr<Vertical> const& r) const {
            return *l < *r;
        }
    };
}

unsigned long long FastFDs::execute() {
    relation_ = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    schema_ = relation_->getSchema();

    if (schema_->getNumColumns() == 0)
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");

    auto start_time = std::chrono::system_clock::now();

    genDiffSets();

    auto elapsed_mills_to_gen_diff_sets =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
    std::cout << "TIME TO DIFF SETS GENERATION: " << elapsed_mills_to_gen_diff_sets.count() << '\n';

    if (diff_sets_.size() == 1 && *diff_sets_.back() == *schema_->emptyVertical) { 
        auto elapsed_milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
        return elapsed_milliseconds.count();
    }

    for (auto const& column : schema_->getColumns()) {
        if (columnContainsOnlyEqualValues(*column)) {
            std::cout << "Registered FD: " << schema_->emptyVertical->toString() << "->" << column->toString() << '\n';
            registerFD(Vertical(), *column);
            continue;
        }

        vector<shared_ptr<Vertical>> diff_sets_mod = getDiffSetsMod(*column);
        assert(!diff_sets_mod.empty());
        if (!(diff_sets_mod.size() == 1 && *diff_sets_mod.back() == *schema_->emptyVertical)) {
            // use vector instead of set?
            set<Column, OrderingComparator> init_ordering = getInitOrdering(diff_sets_mod, *column);
            findCovers(*column, diff_sets_mod, diff_sets_mod, *schema_->emptyVertical, init_ordering);
        }
    }

    auto elapsed_milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);

    return elapsed_milliseconds.count();
}

bool FastFDs::columnContainsOnlyEqualValues(Column const& column) const {
    auto pli = relation_->getColumnData(column.getIndex())->getPositionListIndex();
    bool column_contains_only_equal_values = pli->getNumNonSingletonCluster() == 1 &&
                                             static_cast<unsigned int>(pli->getSize()) == relation_->getNumRows();
    return column_contains_only_equal_values;
}

void FastFDs::findCovers(Column const& attribute, vector<shared_ptr<Vertical>> const& diff_sets_mod,
                         vector<shared_ptr<Vertical>> const& cur_diff_sets, Vertical const& path,
                         set<Column, OrderingComparator> const& ordering) {
    if (ordering.size() == 0 && !cur_diff_sets.empty())
        return; // no FDs here

    if (cur_diff_sets.empty()) {
        if (coverMinimal(path, diff_sets_mod)) {
            std::cout << "Registered FD: " << path.toString() << "->" << attribute.toString() << '\n';
            registerFD(path, attribute);
            return;
        }
        return; // wasted effort, non-minimal result
    }

    for (auto const& column : ordering) {
        vector<shared_ptr<Vertical>> next_diff_sets;
        for (auto const& diff_set : cur_diff_sets) {
            if (!diff_set->contains(column))
                next_diff_sets.push_back(diff_set);
        }

        auto next_ordering = getNextOrdering(next_diff_sets, column, ordering);
        findCovers(attribute, diff_sets_mod, next_diff_sets, *path.Union(column), next_ordering);
    }
}

bool FastFDs::isCover(Vertical const& candidate, vector<shared_ptr<Vertical>> const& sets) const {
    bool covers = true;
    for (auto const& set: sets) {
        if (!set->intersects(candidate)) {
            covers = false;
            break;
        }
    }
    return covers;
}

bool FastFDs::coverMinimal(Vertical const& cover, vector<shared_ptr<Vertical>> const& diff_sets_mod) const {
    for (auto const& column : cover.getColumns()) {
        shared_ptr<Vertical> subset = cover.without(*column);
        bool subset_covers = isCover(*subset, diff_sets_mod);
        if (subset_covers)
            return false; // cover is not minimal
    }
    return true; // cover is minimal
}

bool FastFDs::orderingComp(vector<shared_ptr<Vertical>> const& diff_sets,
                           Column const& l_col, Column const& r_col) const {
        unsigned cov_l = 0;
        unsigned cov_r = 0;

        for (auto& diff_set : diff_sets) {
            if (diff_set->contains(l_col))
                ++cov_l;
            if (diff_set->contains(r_col))
                ++cov_r;
        }
        if (cov_l != cov_r)
            return cov_l > cov_r;
        // l_index < r_index => l > r. Or maybe need to compare columns using their names?
        return l_col.getIndex() < r_col.getIndex();
}

set<Column, FastFDs::OrderingComparator>
FastFDs::getInitOrdering(vector<shared_ptr<Vertical>> const& diff_sets, Column const& attribute) const {
    auto ordering_comp = [&diff_sets, this](Column const& l_col, Column const& r_col) {
        return orderingComp(diff_sets, l_col, r_col);
    };
    set<Column, OrderingComparator> ordering(ordering_comp);
    for (auto const& col : schema_->getColumns()) {
        if (*col != attribute)
            ordering.insert(*col);
    }
    return ordering;
}

set<Column, FastFDs::OrderingComparator>
FastFDs::getNextOrdering(vector<shared_ptr<Vertical>> const& diff_sets, Column const& attribute,
                         set<Column, OrderingComparator> const& cur_ordering) const {
    auto ordering_comp = [&diff_sets, this](Column const& l_col, Column const& r_col) {
        return orderingComp(diff_sets, l_col, r_col);
    };
    set<Column, OrderingComparator> ordering(ordering_comp);

    auto p = cur_ordering.find(attribute);
    assert(p != cur_ordering.end());
    for (++p; p != cur_ordering.end(); ++p) {
        //awful kostil
        for (auto const& diff_set : diff_sets) {
            if (diff_set->contains(*p)) {
                ordering.insert(*p);
                break;
            }
        }
    }
    return ordering;
}

vector<shared_ptr<Vertical>> FastFDs::getDiffSetsMod(Column const& col) const {
    vector<shared_ptr<Vertical>> diff_sets_mod;

    /* diff_sets_ is sorted, before adding next diff_set to
     * diff_sets_mod need to check if diff_sets_mod contains
     * a subset of diff_set, that means that diff_set
     * is not minimal.
     */
    for (auto const& diff_set : diff_sets_) {
        if (diff_set->contains(col)) {
            bool is_miminal = true;
            for (auto const& min_diff_set : diff_sets_mod) {
                if (diff_set->contains(*min_diff_set)) {
                    is_miminal = false;
                    break;
                }
            }
            if (is_miminal)
                diff_sets_mod.push_back(diff_set->without(col));
        }
    }

    #ifdef FASTFDS_DEBUG
        std::cout << "Compute minimal difference sets modulo " << col.toString() << ":\n";
        for (auto& item : diff_sets_mod) {
            std::cout << item->toString() << '\n';
        }
    #endif

    return diff_sets_mod;
}

void FastFDs::genDiffSets() {
    // std::set to get rid of repeating agree sets during inserting
    set<shared_ptr<Vertical>> const agree_sets = genAgreeSets();

    #ifdef FASTFDS_DEBUG
        std::cout << "Agree sets:\n";
        for (auto const& agree_set : agree_sets) {
            std::cout << agree_set->toString() << '\n';
        }
    #endif

    // Complement agree sets to get difference sets
    diff_sets_.reserve(agree_sets.size());
    for (auto const& agree_set : agree_sets) {
        diff_sets_.push_back(agree_set->invert());
    }
    // sort diff_sets_, it will be used further to find minimal difference sets modulo column
    std::sort(diff_sets_.begin(), diff_sets_.end());

    #ifdef FASTFDS_DEBUG
        std::cout << "Compute difference sets:\n";
        for (auto const& diff_set : diff_sets_)
            std::cout << diff_set->toString() << '\n';
    #endif
}

set<vector<int>> FastFDs::getPLIMaxRepresentation() const {
    auto start_time = std::chrono::system_clock::now();
    vector<shared_ptr<ColumnData>> const& columns_data = relation_->getColumnData();
    //inefficient, metanome uses HashSet(std::unordered_set, need hash for vector, boost?)
    auto not_empty_pli = std::find_if(columns_data.begin(), columns_data.end(),
                                      [](shared_ptr<ColumnData> const& c) {
        return c->getPositionListIndex()->getSize() != 0;
    });

    if (not_empty_pli == columns_data.end()) {
        return {};
    }

    set<vector<int>> max_representation((*not_empty_pli)->getPositionListIndex()->getIndex().begin(),
                                        (*not_empty_pli)->getPositionListIndex()->getIndex().end());

    for (auto p = columns_data.begin(); p != columns_data.end(); ++p) {
        if (p == not_empty_pli) //already examined
            continue;

        shared_ptr<PositionListIndex> pli = (*p)->getPositionListIndex();
        if (pli->getSize() != 0)
            calculateSupersets(max_representation, pli->getIndex());
    }

    auto elapsed_mills_to_gen_max_representation =
        std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
    std::cout << "TIME TO MAX REPRESENTATION GENERATION: " << elapsed_mills_to_gen_max_representation.count() << '\n';

    return max_representation;
}

template<FastFDs::AgreeSetsGenMethod method>
set<shared_ptr<Vertical>> FastFDs::genAgreeSets() const {
    set<shared_ptr<Vertical>> agree_sets;

    if constexpr (method == AgreeSetsGenMethod::kUsingVectorOfIDSets) {
        vector<IdentifierSet> identifier_sets;
        set<vector<int>> const max_representation = getPLIMaxRepresentation();

        auto start_time = std::chrono::system_clock::now();

        // compute identifier sets
        // identifier_sets is vector
        std::unordered_set<int> cache;
        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                if (cache.find(*p) != cache.end())
                    continue;
                cache.insert(*p);
                identifier_sets.push_back(IdentifierSet(*relation_, *p));
            }
        }

        auto elapsed_mills_to_gen_id_sets =
            std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
        std::cout << "TIME TO IDENTIFIER SETS GENERATION: " << elapsed_mills_to_gen_id_sets.count() << '\n';

        #ifdef FASTFDS_DEBUG
            std::cout << "Identifier sets:\n";
            for (auto const& id_set : identifier_sets) {
                std::cout << id_set.toString() << '\n';
            }
        #endif

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
        std::unordered_map<int, IdentifierSet> identifier_sets;
        set<vector<int>> const max_representation = getPLIMaxRepresentation();

        auto start_time = std::chrono::system_clock::now();

        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                identifier_sets.emplace(*p, IdentifierSet(*relation_, *p));
            }
        }

        auto elapsed_mills_to_gen_id_sets =
            std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start_time);
        std::cout << "TIME TO IDENTIFIER SETS GENERATION: " << elapsed_mills_to_gen_id_sets.count() << '\n';


        #ifdef FASTFDS_DEBUG
            std::cout << "Identifier sets:\n";
            for (auto const& [index, id_set] : identifier_sets) {
                std::cout << id_set.toString() << '\n';
            }
        #endif


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
        set<vector<int>> const max_representation = getPLIMaxRepresentation();

        // Compute agree sets from maximal representation using getAgreeSet()
        // ~3300 ms on CIPublicHighway700 (Debug build), ~250 ms (Release)
        for (auto const& cluster : max_representation) {
            for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                for (auto q = std::next(p); q != cluster.end(); ++q) {
                    agree_sets.insert(getAgreeSet(relation_->getTuple(*p), relation_->getTuple(*q)));
                }
            }
        }
    } else if constexpr (method == AgreeSetsGenMethod::kUsingGetAgreeSet) {
        vector<shared_ptr<ColumnData>> const columns_data = relation_->getColumnData();

        // Compute agree sets from stripped partitions (simplest method by Wyss)
        // ~40436 ms on CIPublicHighway700 (Debug build)
        for (auto const& column_data : columns_data) {
            shared_ptr<PositionListIndex> pli = column_data->getPositionListIndex();
            for (vector<int> const& cluster : pli->getIndex()) {
                for (auto p = cluster.begin(); p != cluster.end(); ++p) {
                    for (auto q = std::next(p); q != cluster.end(); ++q) {
                        agree_sets.insert(getAgreeSet(relation_->getTuple(*p), relation_->getTuple(*q)));
                    }
                }
            }
        }
    }

    // metanome kostil, doesn't work properly in general
    agree_sets.insert(schema_->emptyVertical);

    return agree_sets;
}

void FastFDs::calculateSupersets(set<vector<int>>& max_representation,
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

shared_ptr<Vertical> FastFDs::getAgreeSet(vector<int> const& tuple1, vector<int> const& tuple2) const {
    assert(tuple1.size() == tuple2.size() && tuple1.size() == relation_->getNumColumns());
    boost::dynamic_bitset<> agree_set_indices(relation_->getNumColumns());

    for (size_t i = 0; i < agree_set_indices.size(); ++i) {
        if (tuple1[i] != 0 && tuple1[i] == tuple2[i]) {
            agree_set_indices.set(i);
        }
    }
    return schema_->getVertical(agree_set_indices);
}