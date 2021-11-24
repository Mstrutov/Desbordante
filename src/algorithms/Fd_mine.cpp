#include "Fd_mine.h"

#include <boost/unordered_map.hpp>
#include <queue>
#include <vector>

unsigned long long Fd_mine::executeInternal() {
    // 1
    schema = relation_->getSchema();
    auto startTime = std::chrono::system_clock::now();

    relationIndices = dynamic_bitset<>(schema->getNumColumns());

    for (size_t columnIndex = 0; columnIndex < schema->getNumColumns(); columnIndex++) {
        dynamic_bitset<> tmp(schema->getNumColumns());
        tmp[columnIndex] = 1;
        relationIndices[columnIndex] = 1;
        candidateSet.insert(std::move(tmp));
    }

    for (auto const& candidate : candidateSet) {
        closure[candidate] = dynamic_bitset<>(schema->getNumColumns());
    }

    // 2
    while (!candidateSet.empty()) {
        for (auto const& candidate : candidateSet) {
            computeNonTrivialClosure(candidate);
            obtainFDandKey(candidate);
        }
        obtainEQSet();
        pruneCandidates();
        generateNextLevelCandidates();
    }

    // 3
    reconstruct();
    display();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    return elapsed_milliseconds.count();
}

void Fd_mine::computeNonTrivialClosure(dynamic_bitset<> const& candidateX) {
    if (!closure.count(candidateX)) {
        closure[candidateX] = dynamic_bitset<>(candidateX.size());
    }
    for (int columnIndex = 0; columnIndex < schema->getNumColumns(); columnIndex++) {
        if ((relationIndices - candidateX - closure[candidateX])[columnIndex]) {
            dynamic_bitset<> candidateXY = candidateX;
            dynamic_bitset<> candidateY(schema->getNumColumns());
            candidateXY[columnIndex] = 1;
            candidateY[columnIndex] = 1;

            if (candidateX.count() == 1) {
                auto candidateXPli = relation_->getColumnData(candidateX.find_first()).getPositionListIndex();
                auto candidateYPli = relation_->getColumnData(columnIndex).getPositionListIndex();

                plis[candidateXY] = candidateXPli->intersect(candidateYPli);

                if (candidateXPli->getNumCluster() == plis[candidateXY]->getNumCluster()) {
                    closure[candidateX][columnIndex] = 1;
                }

                continue;
            }

            if (!plis.count(candidateXY)) {
                auto candidateYPli = relation_->getColumnData(candidateY.find_first()).getPositionListIndex();
                plis[candidateXY] = plis[candidateX]->intersect(candidateYPli);
            }

            if (plis[candidateX]->getNumCluster() == plis[candidateXY]->getNumCluster()) {
                closure[candidateX][columnIndex] = 1;
            }
        }
    }
}

void Fd_mine::obtainFDandKey(dynamic_bitset<> const& candidate) {
    fdSet[candidate] = closure[candidate];
    if (relationIndices == (candidate | closure[candidate])) {
        keySet.insert(candidate);
    }
}

void Fd_mine::obtainEQSet() {
    for (auto const& candidate: candidateSet) {
        for (auto &[lhs, Closure] : fdSet) {
            auto commonAtrs = candidate & lhs;
            if ((candidate - commonAtrs).is_subset_of(Closure) && (lhs - commonAtrs).is_subset_of(closure[candidate])) {
                if (lhs != candidate) {
                    eqSet[lhs].insert(candidate);
                    eqSet[candidate].insert(lhs);
                }
            }
        }
    }
}

void Fd_mine::pruneCandidates() {
    auto it = candidateSet.begin();
    while (it != candidateSet.end()) {
        bool found = false;
        auto const& xi = *it;

        for (auto const& xj : eqSet[xi]) {
            if (candidateSet.find(xj) != candidateSet.end()) {
                it = candidateSet.erase(it);
                found = true;
                break;
            }
        }

        if (found) continue;

        if (keySet.find(xi) != keySet.end()) {
            it = candidateSet.erase(it);
            continue;
        }
        it++;
    }
}

void Fd_mine::generateNextLevelCandidates() {
    std::vector<dynamic_bitset<>> candidates(candidateSet.begin(), candidateSet.end());

    dynamic_bitset<> candidateI;
    dynamic_bitset<> candidateJ;
    dynamic_bitset<> candidateIJ;

    for (size_t i = 0; i < candidates.size(); i++) {
        candidateI = candidates[i];

        for (size_t j = i + 1; j < candidates.size(); j++) {
            candidateJ = candidates[j];

            // apriori-gen
            bool similar = true;
            int setBits = 0;

            for (int k = 0; setBits < candidateI.count() - 1; k++) {
                if (candidateI[k] == candidateJ[k]) {
                    if (candidateI[k]) {
                        setBits++;
                    }
                } else {
                    similar = false;
                    break;
                }
            }
            //

            if (similar) {
                candidateIJ = candidateI | candidateJ;

                if (!(candidateJ).is_subset_of(fdSet[candidateI]) && !(candidateI).is_subset_of(fdSet[candidateJ])) {
                    if (candidateI.count() == 1) {
                        auto candidateIPli = relation_->getColumnData(candidateI.find_first()).getPositionListIndex();
                        auto candidateJPli = relation_->getColumnData(candidateJ.find_first()).getPositionListIndex();
                        plis[candidateIJ] = candidateIPli->intersect(candidateJPli);
                    } else {
                        plis[candidateIJ] = plis[candidateI]->intersect(plis[candidateJ].get());
                    }

                    auto closureIJ = closure[candidateI] | closure[candidateJ];
                    if (relationIndices == (candidateIJ | closureIJ)) {
                        keySet.insert(candidateIJ);
                    } else {
                        candidateSet.insert(candidateIJ);
                    }
                }
            }
        }

        candidateSet.erase(candidateI);
    }
}

void Fd_mine::reconstruct() {
    std::queue<dynamic_bitset<>> queue;
    dynamic_bitset<> generatedLhs(relationIndices.size());
    dynamic_bitset<> generatedLhs_tmp(relationIndices.size());

    for (const auto &[lhs, rhs] : fdSet) {
        std::unordered_map<dynamic_bitset<>, bool> observed;

        observed[lhs] = true;
        auto Rhs = rhs;
        queue.push(lhs);

        for (const auto &[eq, eqset] : eqSet) {
            if (eq.is_subset_of(Rhs)) {
                for (const auto &eqRhs : eqset) {
                    Rhs |= eqRhs;
                }
            }
        }
        bool rhsWillNotChange = false;

        while (!queue.empty()) {
            dynamic_bitset<> currentLhs = queue.front();
            queue.pop();
            int Rhs_count = Rhs.count();
            for (const auto &[eq, eqset] : eqSet) {
                if (!rhsWillNotChange && eq.is_subset_of(Rhs)) {
                    for (const auto &eqRhs : eqset) {
                        Rhs |= eqRhs;
                    }
                }

                if (eq.is_subset_of(currentLhs)) {
                    generatedLhs_tmp = currentLhs - eq;
                    for (const auto &newEq : eqset) {
                        generatedLhs = generatedLhs_tmp;
                        generatedLhs |= newEq;

                        if (!observed[generatedLhs]) {
                            queue.push(generatedLhs);
                            observed[generatedLhs] = true;
                        }
                    }
                }
            }
            if (Rhs_count == Rhs.count()) {
                rhsWillNotChange = true;
            }
        }

        for (auto &[lhs, rbool] : observed) {
            if (final_fdSet.count(lhs)) {
                final_fdSet[lhs] |= Rhs;
            } else {
                final_fdSet[lhs] = Rhs;
            }
        }
    }
}

void Fd_mine::display() {
    unsigned int fd_counter = 0;

    for (auto const& [lhs, rhs] : final_fdSet) {
        for (size_t j = 0; j < rhs.size(); j++) {
            if (!rhs[j] || rhs[j] && lhs[j]) continue;
            std::cout << "Discovered FD: ";
            for (size_t i = 0; i < lhs.size(); i++) {
                if (lhs[i]) std::cout << schema->getColumn(i)->getName() << " ";
            }
            std::cout << "-> " << schema->getColumn(j)->getName() << "\n";
            registerFD(Vertical(schema, lhs), *schema->getColumn(j));
            fd_counter++;
        }
    }
    std::cout << "TOTAL FDs " << fd_counter << "\n";
}
