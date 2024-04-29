#include "fun.h"

#include <easylogging++.h>

namespace algos {

FunQuadruple FunQuadruple::Union(Column const& that) const {
    return FunQuadruple(candidate_.Union(that));
}

FunQuadruple FunQuadruple::Union(Vertical const& that) const {
    return FunQuadruple(candidate_.Union(that));
}

bool FunQuadruple::Contains(FunQuadruple const& that) const {
    return candidate_.Contains(that.candidate_);
}

bool FunQuadruple::Contains(Vertical const& that) const {
    return candidate_.Contains(that);
}

FUN::FUN(std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : PliBasedFDAlgorithm({kDefaultPhaseName}, relation_manager) {}

void FUN::ResetStateFd() {
    fds_.clear();
}

bool FUN::IsKey(FunQuadruple const& l) const {
    return l.GetCount() == relation_->GetNumRows();
}

void FUN::DisplayFD(Level const& l_k_minus_1) {
    for (FunQuadruple const& l : l_k_minus_1) {
        /*  our other algorithms mine l.candidate.GetArity() == 0,
         *  while Metanome's FUN explicitly ignores
         */
        for (Column const* rhs : l.GetClosure().Without(l.GetQuasiclosure()).GetColumns()) {
            fds_.try_emplace(*rhs, std::set<Vertical>());
            bool subset_exists_already = false;
            for (Vertical const& lhs : fds_.at(*rhs)) {
                if (l.Contains(lhs)) {
                    subset_exists_already = true;
                    break;
                }
            }
            if (!subset_exists_already) {
                fds_.at(*rhs).emplace(l.GetCandidate());
            }
        }
    }
}

void FUN::PurePrune(Level const& l_k_minus_1, Level& l_k) const {
    for (auto l = l_k.begin(); l != l_k.end();) {
        bool erased = false;
        for (FunQuadruple const& s : l_k_minus_1) {
            if (l->GetCount() == s.GetCount() && l->Contains(s)) {
                l = l_k.erase(l);
                erased = true;
                break;
            }
        }
        if (!erased) {
            l++;
        }
    }
}

void FUN::ComputeClosure(Level& l_k_minus_1, Level const& l_k) const {
    for (FunQuadruple& l : l_k_minus_1) {
        if (IsKey(l)) {
            continue;
        }
        l.SetClosure(l.GetQuasiclosure());
        for (Column const* a : r_prime_.Without(l.GetQuasiclosure()).GetColumns()) {
            if (FastCount(l_k_minus_1, l_k, l.Union(*a)) == l.GetCount()) {
                l.SetClosure(l.GetClosure().Union(*a));
            }
        }
    }
}

void FUN::ComputeQuasiClosure(Level const& l_k_minus_1, Level& l_k) const {
    for (FunQuadruple& l : l_k) {
        if (IsKey(l)) {
            l.SetClosure(r_);
        }
        l.SetQuasiclosure(l.GetCandidate());
        for (FunQuadruple const& s : l_k_minus_1) {
            if (l.Contains(s)) {
                l.SetQuasiclosure(l.GetQuasiclosure().Union(s.GetClosure()));
            }
        }
    }
}

unsigned long FUN::Count(Vertical const& l) const {
    std::vector<ColumnData> const& column_data = relation_->GetColumnData();
    size_t first_column_index = l.GetColumnIndices().find_first();

    model::PositionListIndex const* pli = column_data.at(first_column_index).GetPositionListIndex();

    if (l.GetColumnIndices().count() == 1) {
        return pli->GetNumCluster();
    }

    //  workaround to avoid auto-destruction of plis
    std::unique_ptr<model::PositionListIndex> holder;

    for (size_t i = l.GetColumnIndices().find_next(first_column_index);
         i != boost::dynamic_bitset<>::npos; i = l.GetColumnIndices().find_next(i)) {
        pli->Intersect(column_data.at(i).GetPositionListIndex()).swap(holder);
        pli = holder.get();
    }

    return pli->GetNumCluster();
}

unsigned long FUN::FastCount(Level const& l_k_minus_1, Level const& l_k,
                             FunQuadruple const& l) const {
    auto position_at_l_k = std::find(l_k.begin(), l_k.end(), l);
    if (position_at_l_k != l_k.end()) {
        return position_at_l_k->GetCount();
    }
    unsigned long max = 0;
    for (FunQuadruple const& l_prime : l_k_minus_1) {
        if (l.Contains(l_prime)) {
            max = std::max(max, l_prime.GetCount());
        }
    }
    return max;
}

std::list<FunQuadruple> FUN::GenerateCandidate(Level const& l_k) const {
    std::set<FunQuadruple> l_k_plus_1;
    for (FunQuadruple const& l_prime : l_k) {
        if (IsKey(l_prime)) {
            continue;
        }
        for (Column const* a : r_prime_.Without(l_prime.GetCandidate()).GetColumns()) {
            FunQuadruple l = l_prime.Union(*a);
            if (l_k_plus_1.find(l) == l_k_plus_1.end()) {
                l.SetCount(Count(l.GetCandidate()));
                l_k_plus_1.emplace(l);
            }
        }
    }
    return {l_k_plus_1.begin(), l_k_plus_1.end()};
}

unsigned long long FUN::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();
    schema_ = relation_->GetSchema();
    double progress_step = kTotalProgressPercent / (schema_->GetNumColumns() + 1);
    AddProgress(progress_step);
    Vertical empty_vertical = *schema_->empty_vertical_;

    r_ = empty_vertical;
    r_prime_ = empty_vertical;
    Level l_k_minus_1{FunQuadruple(empty_vertical)};
    Level l_k;
    for (Column const& a : schema_->GetColumns()) {
        FunQuadruple attribute(a);
        attribute.SetCount(Count(attribute.GetCandidate()));
        l_k.push_back(attribute);
        r_ = r_.Union(a);
        if (!IsKey(attribute)) {
            r_prime_ = r_prime_.Union(a);
        }
        if (attribute.GetCount() == 1) {
            fds_.emplace(a, std::set<Vertical>{empty_vertical});
        }
    }

    while (!l_k.empty()) {
        ComputeClosure(l_k_minus_1, l_k);
        ComputeQuasiClosure(l_k_minus_1, l_k);
        DisplayFD(l_k_minus_1);
        PurePrune(l_k_minus_1, l_k);
        l_k_minus_1 = l_k;
        l_k = GenerateCandidate(l_k);
        AddProgress(progress_step);
    }
    DisplayFD(l_k_minus_1);

    int total_fds = 0;
    for (auto const& [rhs, lverticals] : fds_) {
        for (Vertical const& lhs : lverticals) {
            RegisterFd(lhs, rhs);
            total_fds++;
        }
    }

    SetProgress(kTotalProgressPercent);
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    LOG(INFO) << "Total FD count: " << total_fds;
    LOG(INFO) << "HASH: " << Fletcher16();

    return elapsed_milliseconds.count();
}

}  // namespace algos
