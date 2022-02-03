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
#include "ParallelFor.h"

using std::vector, std::set;

FastFDs::FastFDs(std::filesystem::path const& path, char separator, bool has_header,
                 unsigned int max_lhs, ushort parallelism)
    : PliBasedFDAlgorithm(path, separator, has_header, true,
                          {"Agree sets generation", "Finding minimal covers"}),
      max_lhs_(max_lhs) {
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

unsigned long long FastFDs::ExecuteInternal() {
    schema_ = relation_->GetSchema();
    percent_per_col_ = kTotalProgressPercent / schema_->GetNumColumns();

    auto start_time = std::chrono::system_clock::now();

    GenDiffSets();
    SetProgress(kTotalProgressPercent);
    ToNextProgressPhase();

    auto elapsed_mills_to_gen_diff_sets =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );
    LOG(INFO) << "TIME TO DIFF SETS GENERATION: "
              << elapsed_mills_to_gen_diff_sets.count();

    if (diff_sets_.size() == 1 && diff_sets_.back() == *schema_->empty_vertical_) {
        auto elapsed_milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now() - start_time
            );
        return elapsed_milliseconds.count();
    }

    auto task = [this](std::unique_ptr<Column> const& column) {
        if (ColumnContainsOnlyEqualValues(*column)) {
            LOG(DEBUG) << "Registered FD: " << schema_->empty_vertical_->ToString()
                      << "->" << column->ToString();
            RegisterFd(Vertical(), *column);
            return;
        }

        vector<DiffSet> diff_sets_mod = GetDiffSetsMod(*column);
        assert(!diff_sets_mod.empty());
        if (!(diff_sets_mod.size() == 1 && diff_sets_mod.back() == *schema_->empty_vertical_)) {
            set<Column, OrderingComparator> init_ordering = GetInitOrdering(diff_sets_mod, *column);
            FindCovers(*column, diff_sets_mod, diff_sets_mod,
                       *schema_->empty_vertical_, init_ordering);
        } else {
            AddProgress(percent_per_col_);
        }
    };

    if (threads_num_ > 1) {
        boost::asio::thread_pool pool(threads_num_);

        for (std::unique_ptr<Column> const& column : schema_->GetColumns()) {
            boost::asio::post(pool, [&column, task]() { return task(column); });
        }

        pool.join();
    } else {
        for (std::unique_ptr<Column> const& column : schema_->GetColumns()) {
            task(column);
        }
    }

    SetProgress(kTotalProgressPercent);

    auto elapsed_milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time
        );

    return elapsed_milliseconds.count();
}

bool FastFDs::ColumnContainsOnlyEqualValues(Column const& column) const {
    auto pli = relation_->GetColumnData(column.GetIndex()).GetPositionListIndex();
    bool column_contains_only_equal_values =
        pli->GetNumNonSingletonCluster() == 1 && pli->GetSize() == relation_->GetNumRows();

    return column_contains_only_equal_values;
}

void FastFDs::FindCovers(Column const& attribute, vector<DiffSet> const& diff_sets_mod,
                         vector<DiffSet> const& cur_diff_sets, Vertical const& path,
                         set<Column, OrderingComparator> const& ordering) {
    if (path.GetArity() > max_lhs_) {
        return;
    }

    if (ordering.size() == 0 && !cur_diff_sets.empty()) {
        return; // no FDs here
    }

    if (cur_diff_sets.empty()) {
        if (CoverMinimal(path, diff_sets_mod)) {
            LOG(DEBUG) << "Registered FD: " << path.ToString()
                      << "->" << attribute.ToString();
            RegisterFd(path, attribute);
            return;
        }
        return; // wasted effort, non-minimal result
    }

    for (Column const& column : ordering) {
        vector<DiffSet> next_diff_sets;
        for (DiffSet const& diff_set : cur_diff_sets) {
            if (!diff_set.Contains(column)) {
                next_diff_sets.push_back(diff_set);
            }
        }

        auto next_ordering = GetNextOrdering(next_diff_sets, column, ordering);
        FindCovers(attribute, diff_sets_mod, next_diff_sets, path.Union(column), next_ordering);

        // First FindCovers call, calculate progress
        if (path.GetArity() == 0) {
            AddProgress(percent_per_col_ / ordering.size());
        }
    }
}

bool FastFDs::IsCover(Vertical const& candidate, vector<Vertical> const& sets) const {
    bool covers = true;

    for (Vertical const& set: sets) {
        if (!set.Intersects(candidate)) {
            covers = false;
            break;
        }
    }

    return covers;
}

bool FastFDs::CoverMinimal(Vertical const& cover,
                           vector<DiffSet> const& diff_sets_mod) const {
    for (Column const* column : cover.GetColumns()) {
        Vertical subset = cover.Without(*column);
        bool subset_covers = IsCover(subset, diff_sets_mod);
        if (subset_covers) {
            return false; // cover is not minimal
        }
    }
    return true; // cover is minimal
}

bool FastFDs::OrderingComp(vector<DiffSet> const& diff_sets,
                           Column const& l_col, Column const& r_col) const {
    unsigned cov_l = 0;
    unsigned cov_r = 0;

    for (DiffSet const& diff_set : diff_sets) {
        if (diff_set.Contains(l_col)) {
            ++cov_l;
        }
        if (diff_set.Contains(r_col)) {
            ++cov_r;
        }
    }

    if (cov_l != cov_r) {
        return cov_l > cov_r;
    }

    return l_col > r_col;
}

set<Column, FastFDs::OrderingComparator>
FastFDs::GetInitOrdering(vector<DiffSet> const& diff_sets, Column const& attribute) const {
    auto ordering_comp = [&diff_sets, this](Column const& l_col, Column const& r_col) {
        return OrderingComp(diff_sets, l_col, r_col);
    };
    set<Column, OrderingComparator> ordering(ordering_comp);

    for (auto const& col : schema_->GetColumns()) {
        if (*col != attribute) {
            ordering.insert(*col);
        }
    }

    return ordering;
}

set<Column, FastFDs::OrderingComparator>
FastFDs::GetNextOrdering(vector<DiffSet> const& diff_sets, Column const& attribute,
                         set<Column, OrderingComparator> const& cur_ordering) const {
    auto ordering_comp = [&diff_sets, this](Column const& l_col, Column const& r_col) {
        return OrderingComp(diff_sets, l_col, r_col);
    };
    set<Column, OrderingComparator> ordering(ordering_comp);
    // columns that are contained in at least one diff set
    std::unordered_set<Column> diff_sets_cols;

    for (DiffSet const& diff_set : diff_sets) {
        for (Column const* const col : diff_set.GetColumns()) {
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

/* Metanome uses thread pool here. No need for it because main loop over columns in
 * execute() is parallelized, this approach should be much better.
 */
vector<FastFDs::DiffSet> FastFDs::GetDiffSetsMod(Column const& col) const {
    vector<DiffSet> diff_sets_mod;

    /* diff_sets_ is sorted, before adding next diff_set to
     * diff_sets_mod need to check if diff_sets_mod contains
     * a subset of diff_set, that means that diff_set
     * is not minimal.
     */
    for (DiffSet const& diff_set : diff_sets_) {
        if (diff_set.Contains(col)) {
            bool is_miminal = true;

            for (DiffSet const& min_diff_set : diff_sets_mod) {
                if (diff_set.Contains(min_diff_set)) {
                    is_miminal = false;
                    break;
                }
            }

            if (is_miminal) {
                diff_sets_mod.push_back(diff_set.Without(col));
            }
        }
    }

    LOG(DEBUG) << "Compute minimal difference sets modulo "
               << col.ToString() << ":";
    for (auto& item : diff_sets_mod) {
        LOG(DEBUG) << item.ToString();
    }

    return diff_sets_mod;
}

void FastFDs::GenDiffSets() {
    util::AgreeSetFactory::Configuration c;
    c.threads_num = threads_num_;
    if (threads_num_ > 1) {
        // Not implemented properly, check the description of AgreeSetFactory::GenMcParallel()
        //c.mc_gen_method = MCGenMethod::kParallel;
    }
    util::AgreeSetFactory factory(relation_.get(), c, this);
    util::AgreeSetFactory::SetOfAgreeSets agree_sets = factory.GenAgreeSets();

    LOG(DEBUG) << "Agree sets:";
    for (auto const& agree_set : agree_sets) {
        LOG(DEBUG) << agree_set.ToString();
    }

    // Complement agree sets to get difference sets
    diff_sets_.reserve(agree_sets.size());
    if (threads_num_ > 1) {
        std::mutex m;
        auto const task = [&m, this](util::AgreeSet const& as) {
            DiffSet diff_set = as.Invert();
            std::lock_guard lock(m);
            diff_sets_.push_back(std::move(diff_set));
        };

        util::parallel_foreach(agree_sets.begin(), agree_sets.end(), threads_num_, task);
    } else {
        for (util::AgreeSet const& agree_set : agree_sets) {
            diff_sets_.push_back(agree_set.Invert());
        }
    }

    // sort diff_sets_, it will be used further to find minimal difference sets modulo column
    std::sort(diff_sets_.begin(), diff_sets_.end());

    LOG(DEBUG) << "Compute difference sets:";
    for (auto const& diff_set : diff_sets_) {
        LOG(DEBUG) << diff_set.ToString();
    }
}

