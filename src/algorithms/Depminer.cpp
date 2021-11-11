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

using boost::dynamic_bitset, std::make_shared, std::shared_ptr, std::cout, std::endl, std::setw, std::vector, std::list, std::dynamic_pointer_cast;

unsigned long long Depminer::execute(){

    relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    schema = relation->getSchema();
    if (relation->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }

    auto startTime = std::chrono::system_clock::now();
    
    progressStep = kTotalProgressPercent / schema->getNumColumns();

    //Agree sets (Написано Михаилом)
    AgreeSetFactory agreeSetFactory = AgreeSetFactory(relation.get(), AgreeSetFactory::Configuration(), this);
    std::unordered_set<Vertical> agreeSets = agreeSetFactory.genAgreeSets();
    toNextProgressPhase();

    //maximal sets
    std::vector<CMAXSet> cmaxSets = generateCMAXSets(agreeSets);
    toNextProgressPhase();
    
    //LHS
    auto lhsTime = std::chrono::system_clock::now();
    // 1
    for(const std::unique_ptr<Column>& column : schema->getColumns()){
        lhsForColumn(column, cmaxSets);
        addProgress(progressStep);
    }

    std::chrono::milliseconds lhs_elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - lhsTime);
    cout << "> FD COUNT: " << this->fdCollection_.size() << "\n";
    std::chrono::milliseconds elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    return elapsed_milliseconds.count();
}

std::vector<CMAXSet> Depminer::generateCMAXSets(std::unordered_set<Vertical>& agreeSets){
    auto startTime = std::chrono::system_clock::now();

    std::vector<CMAXSet> cmaxSets;

    for(auto& column : this->schema->getColumns()){
        MAXSet result(*column);

        // finding all sets, which doesn't contain column
        for(Vertical ag : agreeSets){
            if(!ag.contains(*column)){
                result.addCombination(ag);
            }
        }

        // finding max sets
        std::unordered_set<Vertical> superSets;
        std::unordered_set<Vertical> setsDelete;
        bool toAdd = true;

        for(Vertical set : result.getCombinations()){
            for(Vertical superSet : superSets){
                if(set.contains(superSet)){
                    setsDelete.insert(superSet);
                }
                if(toAdd){
                    toAdd = !superSet.contains(set);
                }
            }
            for(auto toDelete : setsDelete){
                superSets.erase(toDelete);
            }
            if(toAdd){
                superSets.insert(set);
            }else{
                toAdd = true;
            }
            setsDelete.clear();
        }
        
        // Inverting MaxSet
        std::unordered_set<Vertical> resultSuperSets;
        for(Vertical const & combination : superSets){
            resultSuperSets.insert(combination.invert());
        }
        result.makeNewCombinations(resultSuperSets);
        cmaxSets.push_back(result);
        addProgress(progressStep);
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    std::cout << "> CMAX SETS COUNT: " << cmaxSets.size() << "\n";

    return std::move(cmaxSets);
}

void Depminer::lhsForColumn(const std::unique_ptr<Column>& column, std::vector<CMAXSet> const & cmaxSets){
        std::unordered_set<Vertical> li;
        // 3
        CMAXSet correct = genFirstLevel(cmaxSets, *column, li);

        auto pli = relation->getColumnData(column->getIndex()).getPositionListIndex();
        bool column_contains_only_equal_values =
            pli->getNumNonSingletonCluster() == 1 && pli->getSize() == relation->getNumRows();
        if (column_contains_only_equal_values) {
            registerFD(Vertical(), *column);
            return;
        }
        
        //4
        while(!li.empty()){
            std::unordered_set<Vertical> liCopy = li;
            //5
            for(Vertical l : li){
                bool isFD = true;
                for(auto combination : correct.getCombinations()){
                    if(!l.intersects(combination)){
                        isFD = false;
                        break;
                    }
                }
                //6
                if(isFD){
                    if(!l.contains(*column)){
                        this->registerFD(l, *column);
                    }
                    liCopy.erase(l);
                }
                if(liCopy.size() == 0){
                    break;
                }
            }
            //7
            li = genNextLevel(liCopy);
        }
}

CMAXSet Depminer::genFirstLevel(std::vector<CMAXSet> const & cmaxSets, Column attribute, std::unordered_set<Vertical> & li){
    CMAXSet correctSet(attribute);
    for(CMAXSet set : cmaxSets){
        if(!(set.getColumn() == attribute)){
            continue;
        }
        correctSet = set;
        for(Vertical combination : correctSet.getCombinations()){
            for(const Column* column : combination.getColumns()){
                if(li.count(Vertical(*column)) == 0)
                    li.insert(Vertical(*column));
            }
        }
        break;
    }
    return correctSet;
}

//Apriori-gen function
std::unordered_set<Vertical> Depminer::genNextLevel(std::unordered_set<Vertical> const& li){
    std::unordered_set<Vertical> ck;
    for(Vertical p : li){
        for(Vertical q : li){
            if(!checkJoin(p, q)){
                continue;
            }
            Vertical candidate(p);
            candidate = candidate.Union(q);
            ck.insert(candidate);
        }
    }
    std::unordered_set<Vertical> result;
    for(Vertical candidate : ck){
        bool prune = false;
        for(auto column : candidate.getColumns()){
            candidate = candidate.invert(Vertical(*column));
            if(li.count(candidate) == 0){
                prune = true;
                break;
            }
            candidate = candidate.invert(Vertical(*column));
        }
        if(!prune){
            result.insert(candidate);
        }
    }
    return std::move(result);
}

bool Depminer::checkJoin(Vertical const& _p, Vertical const& _q){
    dynamic_bitset<> p = _p.getColumnIndices();
    dynamic_bitset<> q = _q.getColumnIndices();

    int pLast = -1, qLast = -1;

    for(int i = 0; i < p.size(); i++){
        pLast = p[i] ? i : pLast;
        qLast = q[i] ? i : qLast;
    }
    if(pLast >= qLast) return false;
    dynamic_bitset<> intersection = p;
    intersection.intersects(q);
    return p.count() == intersection.count()
            && q.count() == intersection.count();
    
    return true;
}
