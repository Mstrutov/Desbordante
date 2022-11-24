#include "partition_table.h"

// see ./LICENSE
#include "algorithms/cfd/miner_node.h"
namespace algos {

[[maybe_unused]] int PartitionTable::database_row_number_;

// Computes intersection
[[maybe_unused]] std::vector<PartitionTidList> PartitionTable::Intersection(
        const PartitionTidList& lhs, const std::vector<PartitionTidList*> rhses) {
    std::unordered_map<int, int> eqIndices(support(lhs.tids));
    std::vector<std::vector<int> > eqClasses(lhs.sets_number);
    // Construct a lookup from tid to equivalence class
    int eix = 0;
    int count = 0;
    for (unsigned ix = 0; ix <= lhs.tids.size(); ix++) {
        count++;
        if (ix == lhs.tids.size() || lhs.tids[ix] == PartitionTidList::SEP) {
            eqClasses[eix++].reserve(count);
            count = 0;
        } else {
            eqIndices[lhs.tids[ix]] = eix + 1;
        }
    }
    std::vector<PartitionTidList> res;
    for (const PartitionTidList* rhs : rhses) {
        res.emplace_back();
        res.back().sets_number = 0;
        res.back().tids.reserve(lhs.tids.size());
        for (unsigned ix = 0; ix <= rhs->tids.size(); ix++) {
            if (ix == rhs->tids.size() || rhs->tids[ix] == PartitionTidList::SEP) {
                for (auto& eqcl : eqClasses) {
                    if (eqcl.size()) {
                        res.back().tids.insert(res.back().tids.end(), eqcl.begin(), eqcl.end());
                        res.back().tids.push_back(PartitionTidList::SEP);
                        res.back().sets_number++;
                        eqcl.clear();
                    }
                }
            } else {
                const int jt = rhs->tids[ix];
                if (eqIndices.count(jt)) {
                    eqClasses[eqIndices[jt] - 1].push_back(jt);
                }
            }
        }
        if (res.back().tids.size() && res.back().tids.back() == PartitionTidList::SEP) {
            res.back().tids.pop_back();
        }
    }
    return res;
}

std::vector<PartitionTidList> PartitionTable::Intersection(
        const PartitionTidList& lhs, const std::vector<const PartitionTidList*> rhses) {
    std::unordered_map<int, int> eqIndices(support(lhs.tids));
    std::vector<std::vector<int> > eqClasses(lhs.sets_number);
    // Construct a lookup from tid to equivalence class
    int eix = 0;
    int count = 0;
    for (unsigned ix = 0; ix <= lhs.tids.size(); ix++) {
        count++;
        if (ix == lhs.tids.size() || lhs.tids[ix] == PartitionTidList::SEP) {
            eqClasses[eix++].reserve(count);
            count = 0;
        } else {
            eqIndices[lhs.tids[ix]] = eix + 1;
        }
    }
    std::vector<PartitionTidList> res;
    for (const PartitionTidList* rhs : rhses) {
        res.emplace_back();
        res.back().sets_number = 0;
        res.back().tids.reserve(lhs.tids.size());
        for (unsigned ix = 0; ix <= rhs->tids.size(); ix++) {
            if (ix == rhs->tids.size() || rhs->tids[ix] == PartitionTidList::SEP) {
                for (auto& eqcl : eqClasses) {
                    if (eqcl.size()) {
                        res.back().tids.insert(res.back().tids.end(), eqcl.begin(), eqcl.end());
                        res.back().tids.push_back(PartitionTidList::SEP);
                        res.back().sets_number++;
                        eqcl.clear();
                    }
                }
            } else {
                const int jt = rhs->tids[ix];
                if (eqIndices.count(jt)) {
                    eqClasses[eqIndices[jt] - 1].push_back(jt);
                }
            }
        }
        if (res.back().tids.size() && res.back().tids.back() == PartitionTidList::SEP) {
            res.back().tids.pop_back();
        }
    }
    return res;
}

[[maybe_unused]] PartitionTidList PartitionTable::Intersection(const PartitionTidList& lhs,
                                                               const PartitionTidList& rhs) {
    std::unordered_map<int, int> eqIndices;
    eqIndices.reserve(lhs.tids.size() + 1 - lhs.sets_number);
    std::vector<std::vector<int> > eqClasses(lhs.sets_number);

    // Construct a lookup from tid to equivalence class
    int eix = 0;
    int count = 0;
    for (unsigned ix = 0; ix <= lhs.tids.size(); ix++) {
        count++;
        if (ix == lhs.tids.size() || lhs.tids[ix] == PartitionTidList::SEP) {
            eqClasses[eix].reserve(count);
            count = 0;
            eix++;
        } else {
            eqIndices[lhs.tids[ix]] = eix + 1;
        }
    }
    // For each rhs partition, for each eq class: spread all tids over eq classes in lhs
    PartitionTidList res;
    res.sets_number = 0;
    for (unsigned ix = 0; ix <= rhs.tids.size(); ix++) {
        if (ix == rhs.tids.size() || rhs.tids[ix] == PartitionTidList::SEP) {
            for (auto& eqcl : eqClasses) {
                if (eqcl.size()) {
                    res.tids.insert(res.tids.end(), eqcl.begin(), eqcl.end());
                    res.tids.push_back(PartitionTidList::SEP);
                    res.sets_number++;
                    eqcl.clear();
                }
            }
        } else {
            const int jt = rhs.tids[ix];
            if (eqIndices[jt]) {
                eqClasses[eqIndices[jt] - 1].push_back(jt);
            }
        }
    }
    if (res.tids.size() && res.tids.back() == PartitionTidList::SEP) {
        res.tids.pop_back();
    }
    return res;
}

[[maybe_unused]] std::vector<int> PartitionTable::OtherPartitionMap(const PartitionTidList& x,
                                                                    const PartitionTidList& xa) {
    std::vector<int> res;

    // Map all partitions in xa to their size; use a random tid as identifier
    std::unordered_map<int, int> bigt;
    bigt.reserve(xa.sets_number);
    int count = 0;
    for (unsigned pi = 0; pi <= xa.tids.size(); pi++) {
        if (xa.tids[pi] == PartitionTidList::SEP || pi == xa.tids.size()) {
            bigt[xa.tids[pi - 1]] = count;
            count = 0;
        } else {
            count++;
        }
    }

    int rep = -1;
    std::vector<int> bins;
    for (unsigned cix = 0; cix < x.tids.size(); cix++) {
        if (x.tids[cix] == PartitionTidList::SEP) {
            // res[eix] = std::make_pair(rep, bins);
            res.push_back(PartitionTidList::SEP);
            rep = -1;
            bins.clear();
        } else {
            int t = x.tids[cix];
            if (bigt.count(t)) {
                if (rep < 0) {
                    rep = t;
                    res.push_back(t);
                }
                res.push_back(bigt[t]);
                // res[eix].first = t;
                // res[eix].second.push_back(bigt[t]);
            }
        }
    }
    return res;
}

[[maybe_unused]] std::vector<std::pair<int, std::vector<int> > > PartitionTable::PartitionMap(
        const PartitionTidList& x, const PartitionTidList& xa) {
    std::vector<std::pair<int, std::vector<int> > > res(x.sets_number);

    // Map all partitions in xa to their size; use a random tid as identifier
    std::unordered_map<int, int> bigt;
    bigt.reserve(xa.sets_number);
    int count = 0;
    for (unsigned pi = 0; pi <= xa.tids.size(); pi++) {
        if (xa.tids[pi] == PartitionTidList::SEP || pi == xa.tids.size()) {
            bigt[xa.tids[pi - 1]] = count;
            count = 0;
        } else {
            count++;
        }
    }

    int eix = 0;
    int rep = -1;
    std::vector<int> bins;
    for (unsigned cix = 0; cix < x.tids.size(); cix++) {
        if (x.tids[cix] == PartitionTidList::SEP) {
            res[eix] = std::make_pair(rep, bins);
            rep = -1;
            bins.clear();
            eix++;
        } else {
            int t = x.tids[cix];
            if (bigt.count(t)) {
                rep = t;
                bins.push_back(bigt[t]);
                // res[eix].first = t;
                // res[eix].second.push_back(bigt[t]);
            }
        }
    }
    return std::vector<std::pair<int, std::vector<int> > >();
}

int PartitionTable::PartitionError(const PartitionTidList& x, const PartitionTidList& xa) {
    int e = 0;

    std::map<int, int> bigt;
    // bigt.reserve(xa.sets_number);
    int count = 0;
    for (unsigned pi = 0; pi <= xa.tids.size(); pi++) {
        if (pi == xa.tids.size() || xa.tids[pi] == PartitionTidList::SEP) {
            bigt[xa.tids[pi - 1]] = count;
            count = 0;
        } else {
            count++;
        }
    }

    count = 0;
    int m = 0;
    for (unsigned cix = 0; cix <= x.tids.size(); cix++) {
        if (cix == x.tids.size() || x.tids[cix] == PartitionTidList::SEP) {
            e += count - m;
            m = 0;
            count = 0;
        } else {
            count++;
            int t = x.tids[cix];
            if (bigt.count(t)) {
                if (bigt[t] > m) {
                    m = bigt[t];
                }
            }
        }
    }
    return e;
}

[[maybe_unused]] bool PartitionTable::ViolatedInCleaned(const PartitionTidList& x,
                                                        const PartitionTidList& xa, int rep) {
    std::unordered_map<int, int> bigt;
    bigt.reserve(xa.sets_number);
    bool first = true;
    for (unsigned pi = 0; pi <= xa.tids.size(); pi++) {
        if (pi == xa.tids.size() || xa.tids[pi] == PartitionTidList::SEP) {
            first = true;
        } else if (first) {
            bigt[xa.tids[pi]] = 1;
            first = false;
        }
    }

    bool inCleaned = false;
    for (unsigned cix = 0; cix <= x.tids.size(); cix++) {
        if (cix == x.tids.size() || x.tids[cix] == PartitionTidList::SEP) {
            inCleaned = false;
        } else {
            int t = x.tids[cix];
            if (bigt.count(t)) {
                if (inCleaned && t < rep) {
                    return true;
                } else if (t < rep) {
                    inCleaned = true;
                }
            }
        }
    }
    return false;
}

[[maybe_unused]] SimpleTidList PartitionTable::Violations(const PartitionTidList& x,
                                                          const PartitionTidList& xa) {
    SimpleTidList res;

    std::set<int> bigt;
    for (unsigned pi = 0; pi <= xa.tids.size(); pi++) {
        if (pi == xa.tids.size() || xa.tids[pi] == PartitionTidList::SEP) {
            bigt.insert(xa.tids[pi - 1]);
        }
    }

    int refs = 0;
    SimpleTidList part;
    for (unsigned cix = 0; cix <= x.tids.size(); cix++) {
        if (cix == x.tids.size() || x.tids[cix] == PartitionTidList::SEP) {
            if (refs != 1) {
                res.insert(res.end(), part.begin(), part.end());
            }
            part.clear();
            refs = 0;
        } else {
            int t = x.tids[cix];
            part.push_back(t);
            if (bigt.count(t)) {
                refs++;
            }
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}

[[maybe_unused]] SimpleTidList PartitionTable::Violations(const PartitionTidList& x,
                                                          const PartitionTidList& xa, int& e) {
    SimpleTidList res;
    e = 0;

    std::map<int, int> bigt;
    int count = 0;
    for (unsigned pi = 0; pi <= xa.tids.size(); pi++) {
        if (pi == xa.tids.size() || xa.tids[pi] == PartitionTidList::SEP) {
            bigt[xa.tids[pi - 1]] = count;
            count = 0;
        } else {
            count++;
        }
    }

    count = 0;
    int m = 0;
    int refs = 0;
    SimpleTidList part;
    for (unsigned cix = 0; cix <= x.tids.size(); cix++) {
        if (cix == x.tids.size() || x.tids[cix] == PartitionTidList::SEP) {
            if (refs > 1) {
                res.insert(res.end(), part.begin(), part.end());
            }
            part.clear();
            refs = 0;
            e += count - m;
            m = 0;
            count = 0;
        } else {
            count++;
            int t = x.tids[cix];
            part.push_back(t);
            if (bigt.count(t)) {
                if (bigt[t] > m) {
                    m = bigt[t];
                }
                refs++;
            }
        }
    }
    std::sort(res.begin(), res.end());
    return res;
}

// Считает partition error?
[[maybe_unused]] int PartitionTable::PartitionError(const std::vector<std::vector<int> >& x,
                                                    const std::vector<std::vector<int> >& xa) {
    int e = 0;
    std::map<int, int> bigt;
    for (unsigned pi = 0; pi < xa.size(); pi++) {
        bigt[xa[pi][0]] = xa[pi].size();
    }

    for (const auto& c : x) {
        int m = 1;
        for (int t : c) {
            if (bigt.count(t)) {
                if (bigt[t] > m) {
                    m = bigt[t];
                }
            }
        }
        e += c.size() - m;
    }
    return e;
}
} //namespace algos
