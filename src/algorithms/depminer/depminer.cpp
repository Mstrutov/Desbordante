#include "algorithms/depminer/depminer.h"

#include <chrono>
#include <list>
#include <memory>

#include <easylogging++.h>

#include "model/column_combination.h"
#include "model/relational_schema.h"
#include "util/agree_set_factory.h"

namespace algos {

Depminer::Depminer(bool request_prepared_data)
    : PliBasedFDAlgorithm({"AgreeSets generation", "Finding CMAXSets", "Finding LHS"},
                          request_prepared_data) {}

using boost::dynamic_bitset, std::make_shared, std::shared_ptr, std::setw, std::vector, std::list,
        std::dynamic_pointer_cast;

unsigned long long Depminer::ExecuteInternal() {
    const auto start_time = std::chrono::system_clock::now();

    schema_ = relation_->GetSchema();

    progress_step_ = kTotalProgressPercent / schema_->GetNumColumns();

    //Agree sets
    const util::AgreeSetFactory agree_set_factory =
        util::AgreeSetFactory(relation_.get(), util::AgreeSetFactory::Configuration(), this);
    const auto agree_sets = agree_set_factory.GenAgreeSets();
    ToNextProgressPhase();

    //maximal sets
    const std::vector<CMAXSet> c_max_cets = GenerateCmaxSets(agree_sets);
    ToNextProgressPhase();

    //LHS
    const auto lhs_time = std::chrono::system_clock::now();
    // 1
    for (auto const& column : schema_->GetColumns()) {
        LhsForColumn(column, c_max_cets);
        AddProgress(progress_step_);
    }

    const auto lhs_elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - lhs_time);
    LOG(INFO) << "> LHS FIND TIME: " << lhs_elapsed_milliseconds.count();
    LOG(INFO) << "> FD COUNT: " << this->fd_collection_.Size();
    const auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

std::vector<CMAXSet> Depminer::GenerateCmaxSets(std::unordered_set<Vertical> const& agree_sets) {
    const auto start_time = std::chrono::system_clock::now();

    std::vector<CMAXSet> c_max_cets;

    for (auto const& column : this->schema_->GetColumns()) {
        CMAXSet result(*column);

        // finding all sets, which doesn't contain column
        for (auto const& ag : agree_sets) {
            if (!ag.Contains(*column)) {
                result.AddCombination(ag);
            }
        }

        // finding max sets
        std::unordered_set<Vertical> super_sets;
        std::unordered_set<Vertical> sets_delete;
        bool to_add = true;

        for (auto const& set : result.GetCombinations()) {
            for (auto const& super_set : super_sets) {
                if (set.Contains(super_set)) {
                    sets_delete.insert(super_set);
                }
                if (to_add) {
                    to_add = !super_set.Contains(set);
                }
            }
            for (auto const& to_delete : sets_delete) {
                super_sets.erase(to_delete);
            }
            if (to_add) {
                super_sets.insert(set);
            } else {
                to_add = true;
            }
            sets_delete.clear();
        }

        // Inverting MaxSet
        std::unordered_set<Vertical> result_super_sets;
        for (auto const& combination : super_sets) {
            result_super_sets.insert(combination.Invert());
        }
        result.MakeNewCombinations(std::move(result_super_sets));
        c_max_cets.push_back(result);
        AddProgress(progress_step_);
    }

    const auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    LOG(INFO) << "> CMAX GENERATION TIME: " << elapsed_milliseconds.count();
    LOG(INFO) << "> CMAX SETS COUNT: " << c_max_cets.size();

    return c_max_cets;
}

void Depminer::LhsForColumn(std::unique_ptr<Column> const& column,
                            std::vector<CMAXSet> const& c_max_cets) {
    std::unordered_set<Vertical> level;
    // 3
    CMAXSet correct = GenFirstLevel(c_max_cets, *column, level);

    const auto pli = relation_->GetColumnData(column->GetIndex()).GetPositionListIndex();
    bool column_contains_only_equal_values = pli->IsConstant();
    if (column_contains_only_equal_values) {
        RegisterFd(Vertical(), *column);
        return;
    }

    // 4
    while (!level.empty()) {
        std::unordered_set<Vertical> level_copy = level;
        // 5
        for (auto const& l : level) {
            bool is_fd = true;
            for (auto const& combination : correct.GetCombinations()) {
                if (!l.Intersects(combination)) {
                    is_fd = false;
                    break;
                }
            }
            // 6
            if (is_fd) {
                if (!l.Contains(*column)) {
                    this->RegisterFd(l, *column);
                }
                level_copy.erase(l);
            }
            if (level_copy.empty()) {
                break;
            }
        }
        // 7
        level = GenNextLevel(level_copy);
    }
}

CMAXSet Depminer::GenFirstLevel(std::vector<CMAXSet> const& c_max_cets, Column const& attribute,
                                std::unordered_set<Vertical>& level) {
    CMAXSet correct_set(attribute);
    for (auto const& set : c_max_cets) {
        if (!(set.GetColumn() == attribute)) {
            continue;
        }
        correct_set = set;
        for (auto const& combination : correct_set.GetCombinations()) {
            for (auto const& column : combination.GetColumns()) {
                if (level.count(Vertical(*column)) == 0)
                    level.insert(Vertical(*column));
            }
        }
        break;
    }
    return correct_set;
}

//Apriori-gen function
std::unordered_set<Vertical> Depminer::GenNextLevel(std::unordered_set<Vertical> const& prev_level) {
    std::unordered_set<Vertical> candidates;
    for (auto const& p : prev_level) {
        for (auto const& q : prev_level) {
            if (!CheckJoin(p, q)) {
                continue;
            }
            Vertical candidate(p);
            candidate = candidate.Union(q);
            candidates.insert(candidate);
        }
    }
    std::unordered_set<Vertical> result;
    for (Vertical candidate : candidates) {
        bool prune = false;
        for (auto const& column : candidate.GetColumns()) {
            candidate = candidate.Invert(Vertical(*column));
            if (prev_level.count(candidate) == 0) {
                prune = true;
                break;
            }
            candidate = candidate.Invert(Vertical(*column));
        }
        if (!prune) {
            result.insert(candidate);
        }
    }
    return result;
}

bool Depminer::CheckJoin(Vertical const& _p, Vertical const& _q) {
    dynamic_bitset<> p = _p.GetColumnIndices();
    dynamic_bitset<> q = _q.GetColumnIndices();

    size_t p_last = -1, q_last = -1;

    for (size_t i = 0; i < p.size(); i++) {
        p_last = p[i] ? i : p_last;
        q_last = q[i] ? i : q_last;
    }
    if (p_last >= q_last) return false;
    dynamic_bitset<> intersection = p;
    intersection.intersects(q);
    return p.count() == intersection.count() && q.count() == intersection.count();
}

}  // namespace algos
