#include "Fd_mine.h"

#include <map>
#include <set>
#include <vector>

#include "ColumnLayoutRelationData.h"
#include "LatticeVertex.h"
#include "Vertical.h"

void Fd_mine::execute() {
    // 1
    relation = ColumnLayoutRelationData::createFrom(inputGenerator, true);
    schema = relation->getSchema();

    r = dynamic_bitset<>(schema->getNumColumns());

    for (size_t columnIndex = 0; columnIndex < schema->getNumColumns();
         columnIndex++) {
        dynamic_bitset<> tmp(schema->getNumColumns());
        tmp[columnIndex] = 1;
        r[columnIndex] = 1;
        candidateSet.insert(std::move(tmp));
    }

    for (auto &xi : candidateSet) {
        closure[xi] = dynamic_bitset<>(schema->getNumColumns());
    }

    // 2

    while (!candidateSet.empty()) {
        for (auto &xi : candidateSet) {
            computeNonTrivialClosure(xi);
            obtainFDandKey(xi);
        }
        obtainEQSet(candidateSet);
        pruneCandidates(candidateSet);
        generateCandidates(candidateSet);
    }

    // 3
    display();
}

void Fd_mine::computeNonTrivialClosure(dynamic_bitset<> xi) {
    closure[xi] = dynamic_bitset<>(xi.size());
    for (int columnIndex = 0; columnIndex < schema->getNumColumns();
         columnIndex++) {
        if ((r - xi - closure[xi])[columnIndex]) {
            dynamic_bitset<> xiy = xi;
            xiy[columnIndex] = 1;

            if (plis[xi] == nullptr) {
                shared_ptr<ColumnData> columnData =
                    relation->getColumnData(xi.find_first());
                plis[xi] = columnData->getPositionListIndex();
            }

            dynamic_bitset<> y(schema->getNumColumns());
            y[columnIndex] = 1;
            if (plis[y] == nullptr) {
                shared_ptr<ColumnData> columnData =
                    relation->getColumnData(columnIndex);
                plis[y] = columnData->getPositionListIndex();
            }

            plis[xiy] = plis[xi]->intersect(plis[y]);

            if (plis[xi]->getSize() == plis[xiy]->getSize()) {
                closure[xi][columnIndex] = 1;
            }
        }
    }
}

void Fd_mine::obtainFDandKey(dynamic_bitset<> xi) {
    fdSet[xi] = closure[xi];
    if (r == (xi | closure[xi])) {
        keySet.insert(xi);
    }
}

void Fd_mine::obtainEQSet(std::set<dynamic_bitset<>> &candidateSet) {
    for (auto &xi : candidateSet) {
        for (auto &[x, xClosure] : fdSet) {
            dynamic_bitset<> z = xi & x;
            if ((xi - z).is_subset_of(xClosure) &&
                (x - z).is_subset_of(closure[xi])) {
                if (x != xi) {
                    eqSet[x].insert(xi);
                    eqSet[xi].insert(x);
                }
            }
        }
    }
}

void Fd_mine::pruneCandidates(std::set<dynamic_bitset<>> &candidateSet) {
    std::set<dynamic_bitset<>>::iterator it = candidateSet.begin();
    while (it != candidateSet.end()) {
        bool found = false;
        const dynamic_bitset<> &xi = *it;

        for (auto &xj : eqSet[xi]) {
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

void Fd_mine::generateCandidates(std::set<dynamic_bitset<>> &candidateSet) {
    std::vector<dynamic_bitset<>> candidates(candidateSet.begin(),
                                             candidateSet.end());

    dynamic_bitset<> xi;
    dynamic_bitset<> xj;
    dynamic_bitset<> xij;

    for (size_t i = 0; i < candidates.size(); i++) {
        xi = candidates[i];

        for (size_t j = i + 1; j < candidates.size(); j++) {
            xj = candidates[j];
            // sets different in two bits (in two arguments)
            if ((xi ^ xj).count() == 2) {
                xij = xi | xj;

                if (!(xj).is_subset_of(fdSet[xi]) &&
                    !(xi).is_subset_of(fdSet[xj])) {
                    if (plis[xi] == nullptr) {
                        shared_ptr<ColumnData> columnData =
                            relation->getColumnData(xi.find_first());
                        plis[xi] = columnData->getPositionListIndex();
                    }
                    if (plis[xj] == nullptr) {
                        shared_ptr<ColumnData> columnData =
                            relation->getColumnData(xj.find_first());
                        plis[xj] = columnData->getPositionListIndex();
                    }

                    plis[xij] = plis[xi]->intersect(plis[xj]);

                    dynamic_bitset<> closureXij = closure[xi] | closure[xj];
                    if (r == (xij | closureXij)) {
                        keySet.insert(xij);
                    } else {
                        candidateSet.insert(xij);
                    }
                }
            }
        }

        candidateSet.erase(xi);
    }
}

void Fd_mine::display() {
    int count_fd = 0;
    for (auto [l, r] : fdSet) {
        for (size_t j = 0; j < r.size(); j++) {
            if (!r[j]) continue;
            std::cout << "Discovered FD: ";
            for (size_t i = 0; i < l.size(); i++) {
                if (l[i]) std::cout << schema->getColumn(i)->getName() << " ";
            }
            std::cout << " -> " << schema->getColumn(j)->getName() << std::endl;
            count_fd++;
        }
    }
    std::cout << "TOTAL FDs " << count_fd << std::endl;
}