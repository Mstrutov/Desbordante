#include "FastFDs.h"

#include <algorithm>

#include <boost/dynamic_bitset.hpp>

#include "AgreeSetFactory.h"

#ifndef NDEBUG
    #define FASTFDS_DEBUG
#endif

#ifdef FASTFDS_DEBUG
#define DEBUG_FASTFDS(fmt, ...) \
            do { fprintf(stderr, fmt, ##__VA_ARGS__); } while (0)
#else
#define DEBUG_FASTFDS(fmt, ...)
#endif

using std::vector, std::set;

unsigned long long FastFDs::execute() {
    relation_ = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    schema_ = relation_->getSchema();

    if (schema_->getNumColumns() == 0) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }

    auto start_time = std::chrono::system_clock::now();

    genDiffSets();

    auto elapsed_mills_to_gen_diff_sets =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );
    std::cout << "TIME TO DIFF SETS GENERATION: "
              << elapsed_mills_to_gen_diff_sets.count() << '\n';

    if (diff_sets_.size() == 1 && diff_sets_.back() == *schema_->emptyVertical) {
        auto elapsed_milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
        return elapsed_milliseconds.count();
    }

    for (auto const& column : schema_->getColumns()) {
        if (columnContainsOnlyEqualValues(*column)) {
            std::cout << "Registered FD: " << schema_->emptyVertical->toString()
                      << "->" << column->toString() << '\n';
            registerFD(Vertical(), *column);
            continue;
        }

        vector<DiffSet> diff_sets_mod = getDiffSetsMod(*column);
        assert(!diff_sets_mod.empty());
        if (!(diff_sets_mod.size() == 1 && diff_sets_mod.back() == *schema_->emptyVertical)) {
            // use vector instead of set?
            set<Column, OrderingComparator> init_ordering = getInitOrdering(diff_sets_mod, *column);
            findCovers(*column, diff_sets_mod, diff_sets_mod, *schema_->emptyVertical, init_ordering);
        }
    }

    auto elapsed_milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );

    return elapsed_milliseconds.count();
}

bool FastFDs::columnContainsOnlyEqualValues(Column const& column) const {
    auto pli = relation_->getColumnData(column.getIndex()).getPositionListIndex();
    bool column_contains_only_equal_values =
        pli->getNumNonSingletonCluster() == 1 && pli->getSize() == relation_->getNumRows();

    return column_contains_only_equal_values;
}

void FastFDs::findCovers(Column const& attribute, vector<DiffSet> const& diff_sets_mod,
                         vector<DiffSet> const& cur_diff_sets, Vertical const& path,
                         set<Column, OrderingComparator> const& ordering) {
    if (ordering.size() == 0 && !cur_diff_sets.empty()) {
        return; // no FDs here
    }

    if (cur_diff_sets.empty()) {
        if (coverMinimal(path, diff_sets_mod)) {
            std::cout << "Registered FD: " << path.toString()
                      << "->" << attribute.toString() << '\n';
            registerFD(path, attribute);
            return;
        }
        return; // wasted effort, non-minimal result
    }

    for (Column const& column : ordering) {
        vector<DiffSet> next_diff_sets;
        for (DiffSet const& diff_set : cur_diff_sets) {
            if (!diff_set.contains(column)) {
                next_diff_sets.push_back(diff_set);
            }
        }

        auto next_ordering = getNextOrdering(next_diff_sets, column, ordering);
        findCovers(attribute, diff_sets_mod, next_diff_sets, path.Union(column), next_ordering);
    }
}

bool FastFDs::isCover(Vertical const& candidate, vector<Vertical> const& sets) const {
    bool covers = true;

    for (Vertical const& set: sets) {
        if (!set.intersects(candidate)) {
            covers = false;
            break;
        }
    }

    return covers;
}

bool FastFDs::coverMinimal(Vertical const& cover,
                           vector<DiffSet> const& diff_sets_mod) const {
    for (Column const* column : cover.getColumns()) {
        Vertical subset = cover.without(*column);
        bool subset_covers = isCover(subset, diff_sets_mod);
        if (subset_covers) {
            return false; // cover is not minimal
        }
    }
    return true; // cover is minimal
}

bool FastFDs::orderingComp(vector<DiffSet> const& diff_sets,
                           Column const& l_col, Column const& r_col) const {
    unsigned cov_l = 0;
    unsigned cov_r = 0;

    for (DiffSet const& diff_set : diff_sets) {
        if (diff_set.contains(l_col)) {
            ++cov_l;
        }
        if (diff_set.contains(r_col)) {
            ++cov_r;
        }
    }

    if (cov_l != cov_r) {
        return cov_l > cov_r;
    }

    return l_col > r_col;
}

set<Column, FastFDs::OrderingComparator>
FastFDs::getInitOrdering(vector<DiffSet> const& diff_sets, Column const& attribute) const {
    auto ordering_comp = [&diff_sets, this](Column const& l_col, Column const& r_col) {
        return orderingComp(diff_sets, l_col, r_col);
    };
    set<Column, OrderingComparator> ordering(ordering_comp);

    for (auto const& col : schema_->getColumns()) {
        if (*col != attribute) {
            ordering.insert(*col);
        }
    }

    return ordering;
}

set<Column, FastFDs::OrderingComparator>
FastFDs::getNextOrdering(vector<DiffSet> const& diff_sets, Column const& attribute,
                         set<Column, OrderingComparator> const& cur_ordering) const {
    auto ordering_comp = [&diff_sets, this](Column const& l_col, Column const& r_col) {
        return orderingComp(diff_sets, l_col, r_col);
    };
    set<Column, OrderingComparator> ordering(ordering_comp);
    // columns that are contained in at least one diff set
    std::unordered_set<Column> diff_sets_cols;

    for (DiffSet const& diff_set : diff_sets) {
        for (Column const* const col : diff_set.getColumns()) {
            diff_sets_cols.insert(*col);
        }
    }

    auto p = cur_ordering.find(attribute);
    assert(p != cur_ordering.end());
    for (++p; p != cur_ordering.end(); ++p) {
        if (diff_sets_cols.find(*p) != diff_sets_cols.end()) {
            ordering.insert(*p);
        }
    }
    return ordering;
}

vector<FastFDs::DiffSet> FastFDs::getDiffSetsMod(Column const& col) const {
    vector<DiffSet> diff_sets_mod;

    /* diff_sets_ is sorted, before adding next diff_set to
     * diff_sets_mod need to check if diff_sets_mod contains
     * a subset of diff_set, that means that diff_set
     * is not minimal.
     */
    for (DiffSet const& diff_set : diff_sets_) {
        if (diff_set.contains(col)) {
            bool is_miminal = true;

            for (DiffSet const& min_diff_set : diff_sets_mod) {
                if (diff_set.contains(min_diff_set)) {
                    is_miminal = false;
                    break;
                }
            }

            if (is_miminal) {
                diff_sets_mod.push_back(diff_set.without(col));
            }
        }
    }

    DEBUG_FASTFDS("Compute minimal difference sets modulo %s:\n", col.toString().c_str());
    for (auto& item : diff_sets_mod) {
         DEBUG_FASTFDS("%s\n", item.toString().c_str());
    }

    return diff_sets_mod;
}

void FastFDs::genDiffSets() {
    AgreeSetFactory factory(relation_.get());
    AgreeSetFactory::SetOfAgreeSets agree_sets = factory.genAgreeSets();

    DEBUG_FASTFDS("Agree sets:\n");
    for (auto const& agree_set : agree_sets) {
        DEBUG_FASTFDS("%s\n", agree_set.toString().c_str());
    }

    // Complement agree sets to get difference sets
    diff_sets_.reserve(agree_sets.size());
    for (AgreeSet const& agree_set : agree_sets) {
        diff_sets_.push_back(agree_set.invert());
    }
    // sort diff_sets_, it will be used further to find minimal difference sets modulo column
    std::sort(diff_sets_.begin(), diff_sets_.end());

    DEBUG_FASTFDS("Compute difference sets:\n");
    for (auto const& diff_set : diff_sets_)
        DEBUG_FASTFDS("%s\n", diff_set.toString().c_str());
}
