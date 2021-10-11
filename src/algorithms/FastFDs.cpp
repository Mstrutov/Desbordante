#include "FastFDs.h"

#include <algorithm>
#include <mutex>
#include <thread>

#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/thread.hpp>

#include "AgreeSetFactory.h"
#include "logging/easylogging++.h"

using std::vector, std::set;

FastFDs::FastFDs(std::filesystem::path const& path,
                 char separator, bool hasHeader,
                 unsigned int max_lhs, ushort parallelism) :
    FDAlgorithm(path, separator, hasHeader,
                { "Agree sets generation", "Finding minimal covers" }), max_lhs_(max_lhs) {
    if (parallelism == 0) {
        threads_num_ = std::thread::hardware_concurrency();
        if (threads_num_ == 0) {
            throw std::runtime_error("Unable to detect number of concurrent"
                                     " threads supported. Specify it manually.");
        }
    } else {
        threads_num_ = parallelism;
    }
}


// Should be FDAlgorithm method I think
void FastFDs::registerFD(Vertical lhs, Column rhs) {
    if (threads_num_ > 1) {
        boost::lock_guard<boost::mutex> lock(register_mutex_);
        FDAlgorithm::registerFD(std::move(lhs), std::move(rhs));
    } else {
        FDAlgorithm::registerFD(std::move(lhs), std::move(rhs));
    }
}

unsigned long long FastFDs::execute() {
    relation_ = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    schema_ = relation_->getSchema();
    percent_per_col_ = kTotalProgressPercent / schema_->getNumColumns();

    if (schema_->getNumColumns() == 0) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }

    auto start_time = std::chrono::system_clock::now();

    genDiffSets();
    setProgress(kTotalProgressPercent);
    toNextProgressPhase();

    auto elapsed_mills_to_gen_diff_sets =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );
    LOG(INFO) << "TIME TO DIFF SETS GENERATION: "
              << elapsed_mills_to_gen_diff_sets.count();

    if (diff_sets_.size() == 1 && diff_sets_.back() == *schema_->emptyVertical) {
        auto elapsed_milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
        return elapsed_milliseconds.count();
    }

    auto task = [this](std::unique_ptr<Column> const& column) {
        if (columnContainsOnlyEqualValues(*column)) {
            LOG(INFO) << "Registered FD: " << schema_->emptyVertical->toString()
                      << "->" << column->toString();
            registerFD(Vertical(), *column);
            return;
        }

        vector<DiffSet> diff_sets_mod = getDiffSetsMod(*column);
        assert(!diff_sets_mod.empty());
        if (!(diff_sets_mod.size() == 1 && diff_sets_mod.back() == *schema_->emptyVertical)) {
            set<Column, OrderingComparator> init_ordering = getInitOrdering(diff_sets_mod, *column);
            findCovers(*column, diff_sets_mod, diff_sets_mod,
                       *schema_->emptyVertical, init_ordering);
        } else {
            addProgress(percent_per_col_);
        }
    };

    if (threads_num_ > 1) {
        boost::asio::thread_pool pool(threads_num_);

        for (std::unique_ptr<Column> const& column : schema_->getColumns()) {
            boost::asio::post(pool, [&column, task](){ return task(column); });
        }

        pool.join();
    } else {
        for (std::unique_ptr<Column> const& column : schema_->getColumns()) {
            task(column);
        }
    }

    setProgress(kTotalProgressPercent);

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
    if (path.getArity() > max_lhs_) {
        return;
    }

    if (ordering.size() == 0 && !cur_diff_sets.empty()) {
        return; // no FDs here
    }

    if (cur_diff_sets.empty()) {
        if (coverMinimal(path, diff_sets_mod)) {
            LOG(INFO) << "Registered FD: " << path.toString()
                      << "->" << attribute.toString();
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

        // First findCovers call, calculate progress
        if (path.getArity() == 0) {
            addProgress(percent_per_col_ / ordering.size());
        }
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

    LOG(DEBUG) << "Compute minimal difference sets modulo "
               << col.toString() << ":";
    for (auto& item : diff_sets_mod) {
         LOG(DEBUG) << item.toString();
    }

    return diff_sets_mod;
}

void FastFDs::genDiffSets() {
    AgreeSetFactory::Configuration c;
    c.threads_num = threads_num_;
    if (threads_num_ > 1) {
        // Not implemented properly, check the description of AgreeSetFactory::genMCParallel()
        //c.mc_gen_method = MCGenMethod::kParallel;
    }
    AgreeSetFactory factory(relation_.get(), c, this);
    AgreeSetFactory::SetOfAgreeSets agree_sets = factory.genAgreeSets();

    LOG(DEBUG) << "Agree sets:";
    for (auto const& agree_set : agree_sets) {
        LOG(DEBUG) << agree_set.toString();
    }

    // Complement agree sets to get difference sets
    diff_sets_.reserve(agree_sets.size());
    if (threads_num_ > 1) {
        size_t agree_sets_per_thread = agree_sets.size() / threads_num_;
        std::mutex m;
        std::vector<std::thread> threads;
        auto task = [&m, this](AgreeSetFactory::SetOfAgreeSets::const_iterator first,
                               AgreeSetFactory::SetOfAgreeSets::const_iterator last) {
            for (; first != last; ++first) {
                DiffSet diff_set = first->invert();
                std::lock_guard lock(m);
                diff_sets_.push_back(std::move(diff_set));
            }
        };

        threads.reserve(threads_num_);

        auto p = agree_sets.begin();
        auto q = agree_sets.end();
        for (ushort i = 0; i < threads_num_; ++i) {
            auto prev = p;

            if (i != threads_num_ - 1) {
                std::advance(p, agree_sets_per_thread);
            } else {
                p = q;
            }
            threads.emplace_back(task, prev, p);
        }

        for (auto& thread : threads) {
            thread.join();
        }
    } else {
        for (AgreeSet const& agree_set : agree_sets) {
            diff_sets_.push_back(agree_set.invert());
        }
    }

    // sort diff_sets_, it will be used further to find minimal difference sets modulo column
    std::sort(diff_sets_.begin(), diff_sets_.end());

    LOG(DEBUG) << "Compute difference sets:";
    for (auto const& diff_set : diff_sets_) {
        LOG(DEBUG) << diff_set.toString();
    }
}

