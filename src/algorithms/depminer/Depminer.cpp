#include "Depminer.h"

#include <chrono>
#include <iostream>
#include <iomanip>
#include <list>
#include <memory>

#include "ColumnCombination.h"
#include "ColumnData.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "AgreeSetFactory.h"
#include "logging/easylogging++.h"

using boost::dynamic_bitset, std::make_shared, std::shared_ptr, std::setw, std::vector, std::list, std::dynamic_pointer_cast;

unsigned long long Depminer::executeInternal() {

    const auto startTime = std::chrono::system_clock::now();
    schema = relation_->getSchema();
    
    progressStep = kTotalProgressPercent / schema->getNumColumns();

    //Agree sets
    const util::AgreeSetFactory agreeSetFactory =util:: AgreeSetFactory(relation_.get(), util::AgreeSetFactory::Configuration(), this);
    const auto agreeSets = agreeSetFactory.genAgreeSets();
    toNextProgressPhase();

    //maximal sets
    const std::vector<CMAXSet> cmaxSets = generateCMAXSets(agreeSets);
    toNextProgressPhase();
    
    //LHS
    const auto lhsTime = std::chrono::system_clock::now();
    // 1
    for (auto const& column : schema->getColumns()) {
        lhsForColumn(column, cmaxSets);
        addProgress(progressStep);
    }

    const auto lhs_elapsed_milliseconds 
        = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - lhsTime);
    LOG(INFO) << "> LHS FIND TIME: " << lhs_elapsed_milliseconds.count();
    LOG(INFO) << "> FD COUNT: " << this->fdCollection_.size();
    const auto elapsed_milliseconds 
        = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    return elapsed_milliseconds.count();
}

std::vector<CMAXSet> Depminer::generateCMAXSets(std::unordered_set<Vertical> const& agreeSets) {
    const auto startTime = std::chrono::system_clock::now();

    std::vector<CMAXSet> cmaxSets;

    for (auto const& column : this->schema->getColumns()) {
        CMAXSet result(*column);

        // finding all sets, which doesn't contain column
        for (auto const& ag : agreeSets) {
            if (!ag.contains(*column)) {
                result.addCombination(ag);
            }
        }

        // finding max sets
        std::unordered_set<Vertical> superSets;
        std::unordered_set<Vertical> setsDelete;
        bool toAdd = true;

        for (auto const& set : result.getCombinations()) {
            for (auto const& superSet : superSets) {
                if (set.contains(superSet)) {
                    setsDelete.insert(superSet);
                }
                if (toAdd) {
                    toAdd = !superSet.contains(set);
                }
            }
            for (auto const& toDelete : setsDelete) {
                superSets.erase(toDelete);
            }
            if (toAdd) {
                superSets.insert(set);
            } else {
                toAdd = true;
            }
            setsDelete.clear();
        }
        
        // Inverting MaxSet
        std::unordered_set<Vertical> resultSuperSets;
        for (auto const& combination : superSets) {
            resultSuperSets.insert(combination.invert());
        }
        result.makeNewCombinations(std::move(resultSuperSets));
        cmaxSets.push_back(result);
        addProgress(progressStep);
    }

    const auto elapsed_milliseconds 
        = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    LOG(INFO) << "> CMAX GENERATION TIME: " << elapsed_milliseconds.count();
    LOG(INFO) << "> CMAX SETS COUNT: " << cmaxSets.size();

    return cmaxSets;
}

void Depminer::lhsForColumn(std::unique_ptr<Column> const& column, std::vector<CMAXSet> const& cmaxSets) {
        std::unordered_set<Vertical> level;
        // 3
        CMAXSet correct = genFirstLevel(cmaxSets, *column, level);

        const auto pli = relation_->getColumnData(column->getIndex()).getPositionListIndex();
        bool column_contains_only_equal_values =
            pli->getNumNonSingletonCluster() == 1 && pli->getSize() == relation_->getNumRows();
        if (column_contains_only_equal_values) {
            registerFD(Vertical(), *column);
            return;
        }
        
        //4
        while (!level.empty()) {
            std::unordered_set<Vertical> levelCopy = level;
            //5
            for (auto const& l : level) {
                bool isFD = true;
                for (auto const& combination : correct.getCombinations()) {
                    if (!l.intersects(combination)) {
                        isFD = false;
                        break;
                    }
                }
                //6
                if (isFD) {
                    if (!l.contains(*column)) {
                        this->registerFD(l, *column);
                    }
                    levelCopy.erase(l);
                }
                if (levelCopy.empty()) {
                    break;
                }
            }
            //7
            level = genNextLevel(levelCopy);
        }
}

CMAXSet Depminer::genFirstLevel(std::vector<CMAXSet> const& cmaxSets, Column const& attribute, std::unordered_set<Vertical>& level) {
    CMAXSet correctSet(attribute);
    for (auto const& set : cmaxSets) {
        if (!(set.getColumn() == attribute)) {
            continue;
        }
        correctSet = set;
        for (auto const& combination : correctSet.getCombinations()) {
            for (auto const& column : combination.getColumns()) {
                if (level.count(Vertical(*column)) == 0)
                    level.insert(Vertical(*column));
            }
        }
        break;
    }
    return correctSet;
}

//Apriori-gen function
std::unordered_set<Vertical> Depminer::genNextLevel(std::unordered_set<Vertical> const& prevLevel) {
    std::unordered_set<Vertical> candidates;
    for (auto const& p : prevLevel) {
        for (auto const& q : prevLevel) {
            if (!checkJoin(p, q)) {
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
        for (auto const& column : candidate.getColumns()) {
            candidate = candidate.invert(Vertical(*column));
            if (prevLevel.count(candidate) == 0) {
                prune = true;
                break;
            }
            candidate = candidate.invert(Vertical(*column));
        }
        if (!prune) {
            result.insert(candidate);
        }
    }
    return result;
}

bool Depminer::checkJoin(Vertical const& _p, Vertical const& _q) {
    dynamic_bitset<> p = _p.getColumnIndices();
    dynamic_bitset<> q = _q.getColumnIndices();

    size_t pLast = -1, qLast = -1;

    for (size_t i = 0; i < p.size(); i++) {
        pLast = p[i] ? i : pLast;
        qLast = q[i] ? i : qLast;
    }
    if (pLast >= qLast) return false;
    dynamic_bitset<> intersection = p;
    intersection.intersects(q);
    return p.count() == intersection.count()
            && q.count() == intersection.count();
}
