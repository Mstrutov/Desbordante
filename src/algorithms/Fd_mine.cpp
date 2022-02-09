#include "Fd_mine.h"

#include <boost/unordered_map.hpp>
#include <queue>
#include <vector>

unsigned long long Fd_mine::ExecuteInternal() {
    // 1
    schema_ = relation_->GetSchema();
    auto start_time = std::chrono::system_clock::now();

    relation_indices_ = dynamic_bitset<>(schema_->GetNumColumns());

    for (size_t column_index = 0; column_index < schema_->GetNumColumns(); column_index++) {
        dynamic_bitset<> tmp(schema_->GetNumColumns());
        tmp[column_index] = 1;
        relation_indices_[column_index] = 1;
        candidate_set_.insert(std::move(tmp));
    }

    for (auto const& candidate : candidate_set_) {
        closure_[candidate] = dynamic_bitset<>(schema_->GetNumColumns());
    }

    // 2
    while (!candidate_set_.empty()) {
        for (auto const& candidate : candidate_set_) {
            ComputeNonTrivialClosure(candidate);
            ObtainFDandKey(candidate);
        }
        ObtainEqSet();
        PruneCandidates();
        GenerateNextLevelCandidates();
    }

    // 3
    Reconstruct();
    Display();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

void Fd_mine::ComputeNonTrivialClosure(dynamic_bitset<> const& xi) {
    if (!closure_.count(xi)) {
        closure_[xi] = dynamic_bitset<>(xi.size());
    }
    for (size_t column_index = 0; column_index < schema_->GetNumColumns(); column_index++) {
        if ((relation_indices_ - xi - closure_[xi])[column_index]) {
            dynamic_bitset<> candidate_xy = xi;
            dynamic_bitset<> candidate_y(schema_->GetNumColumns());
            candidate_xy[column_index] = 1;
            candidate_y[column_index] = 1;

            if (xi.count() == 1) {
                auto candidate_x_pli =
                    relation_->GetColumnData(xi.find_first()).GetPositionListIndex();
                auto candidate_y_pli =
                    relation_->GetColumnData(column_index).GetPositionListIndex();

                plis_[candidate_xy] = candidate_x_pli->Intersect(candidate_y_pli);

                if (candidate_x_pli->GetNumCluster() == plis_[candidate_xy]->GetNumCluster()) {
                    closure_[xi][column_index] = 1;
                }

                continue;
            }

            if (!plis_.count(candidate_xy)) {
                auto candidate_y_pli =
                    relation_->GetColumnData(candidate_y.find_first()).GetPositionListIndex();
                plis_[candidate_xy] = plis_[xi]->Intersect(candidate_y_pli);
            }

            if (plis_[xi]->GetNumCluster() == plis_[candidate_xy]->GetNumCluster()) {
                closure_[xi][column_index] = 1;
            }
        }
    }
}

void Fd_mine::ObtainFDandKey(dynamic_bitset<> const& xi) {
    fd_set_[xi] = closure_[xi];
    if (relation_indices_ == (xi | closure_[xi])) {
        key_set_.insert(xi);
    }
}

void Fd_mine::ObtainEqSet() {
    for (auto const& candidate: candidate_set_) {
        for (auto &[lhs, closure] : fd_set_) {
            auto common_atrs = candidate & lhs;
            if ((candidate - common_atrs).is_subset_of(closure) &&
                (lhs - common_atrs).is_subset_of(closure_[candidate])) {
                if (lhs != candidate) {
                    eq_set_[lhs].insert(candidate);
                    eq_set_[candidate].insert(lhs);
                }
            }
        }
    }
}

void Fd_mine::PruneCandidates() {
    auto it = candidate_set_.begin();
    while (it != candidate_set_.end()) {
        bool found = false;
        auto const& xi = *it;

        for (auto const& xj : eq_set_[xi]) {
            if (candidate_set_.find(xj) != candidate_set_.end()) {
                it = candidate_set_.erase(it);
                found = true;
                break;
            }
        }

        if (found) continue;

        if (key_set_.find(xi) != key_set_.end()) {
            it = candidate_set_.erase(it);
            continue;
        }
        it++;
    }
}

void Fd_mine::GenerateNextLevelCandidates() {
    std::vector<dynamic_bitset<>> candidates(candidate_set_.begin(), candidate_set_.end());

    dynamic_bitset<> candidate_i;
    dynamic_bitset<> candidate_j;
    dynamic_bitset<> candidate_ij;

    for (size_t i = 0; i < candidates.size(); i++) {
        candidate_i = candidates[i];

        for (size_t j = i + 1; j < candidates.size(); j++) {
            candidate_j = candidates[j];

            // apriori-gen
            bool similar = true;
            size_t set_bits = 0;

            assert(candidate_i.count() != 0);
            for (size_t k = 0; set_bits < candidate_i.count() - 1; k++) {
                if (candidate_i[k] == candidate_j[k]) {
                    if (candidate_i[k]) {
                        set_bits++;
                    }
                } else {
                    similar = false;
                    break;
                }
            }
            //

            if (similar) {
                candidate_ij = candidate_i | candidate_j;

                if (!(candidate_j).is_subset_of(fd_set_[candidate_i]) &&
                    !(candidate_i).is_subset_of(fd_set_[candidate_j])) {
                    if (candidate_i.count() == 1) {
                        auto candidate_i_pli = relation_->GetColumnData(candidate_i.find_first())
                                                   .GetPositionListIndex();
                        auto candidate_j_pli = relation_->GetColumnData(candidate_j.find_first())
                                                   .GetPositionListIndex();
                        plis_[candidate_ij] = candidate_i_pli->Intersect(candidate_j_pli);
                    } else {
                        plis_[candidate_ij] =
                            plis_[candidate_i]->Intersect(plis_[candidate_j].get());
                    }

                    auto closure_ij = closure_[candidate_i] | closure_[candidate_j];
                    if (relation_indices_ == (candidate_ij | closure_ij)) {
                        key_set_.insert(candidate_ij);
                    } else {
                        candidate_set_.insert(candidate_ij);
                    }
                }
            }
        }

        candidate_set_.erase(candidate_i);
    }
}

void Fd_mine::Reconstruct() {
    std::queue<dynamic_bitset<>> queue;
    dynamic_bitset<> generated_lhs(relation_indices_.size());
    dynamic_bitset<> generated_lhs_tmp(relation_indices_.size());

    for (const auto &[lhs, rhs] : fd_set_) {
        std::unordered_map<dynamic_bitset<>, bool> observed;

        observed[lhs] = true;
        auto rhs_copy = rhs;
        queue.push(lhs);

        for (const auto &[eq, eqset] : eq_set_) {
            if (eq.is_subset_of(rhs_copy)) {
                for (const auto& eq_rhs : eqset) {
                    rhs_copy |= eq_rhs;
                }
            }
        }
        bool rhs_will_not_change = false;

        while (!queue.empty()) {
            dynamic_bitset<> current_lhs = queue.front();
            queue.pop();
            size_t rhs_count = rhs_copy.count();
            for (const auto &[eq, eqset] : eq_set_) {
                if (!rhs_will_not_change && eq.is_subset_of(rhs_copy)) {
                    for (const auto& eq_rhs : eqset) {
                        rhs_copy |= eq_rhs;
                    }
                }

                if (eq.is_subset_of(current_lhs)) {
                    generated_lhs_tmp = current_lhs - eq;
                    for (const auto& new_eq : eqset) {
                        generated_lhs = generated_lhs_tmp;
                        generated_lhs |= new_eq;

                        if (!observed[generated_lhs]) {
                            queue.push(generated_lhs);
                            observed[generated_lhs] = true;
                        }
                    }
                }
            }
            if (rhs_count == rhs_copy.count()) {
                rhs_will_not_change = true;
            }
        }

        for (auto &[lhs, rbool] : observed) {
            if (final_fd_set_.count(lhs)) {
                final_fd_set_[lhs] |= rhs_copy;
            } else {
                final_fd_set_[lhs] = rhs_copy;
            }
        }
    }
}

void Fd_mine::Display() {
    unsigned int fd_counter = 0;

    for (auto const& [lhs, rhs] : final_fd_set_) {
        for (size_t j = 0; j < rhs.size(); j++) {
            if (!rhs[j] || (rhs[j] && lhs[j])) {
                continue;
            }
            std::cout << "Discovered FD: ";
            for (size_t i = 0; i < lhs.size(); i++) {
                if (lhs[i]) {
                    std::cout << schema_->GetColumn(i)->GetName() << " ";
                }
            }
            std::cout << "-> " << schema_->GetColumn(j)->GetName() << "\n";
            RegisterFd(Vertical(schema_, lhs), *schema_->GetColumn(j));
            fd_counter++;
        }
    }
    std::cout << "TOTAL FDs " << fd_counter << "\n";
}

