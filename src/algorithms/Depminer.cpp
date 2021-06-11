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
#include "algorithms/depminer/util/CMAXGen.h"
#include "AgreeSetFactory.h"

using boost::dynamic_bitset, std::make_shared, std::shared_ptr, std::cout, std::endl, std::setw, std::vector, std::list, std::dynamic_pointer_cast;

std::set<Vertical> unorderedToOrdered(std::unordered_set<Vertical> uset){
    std::set<Vertical> result;
    for(auto el : uset){
        result.insert(el);
    }
    return result;
}

unsigned long long Depminer::execute(){

    std::unique_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    const RelationalSchema* schema = relation->getSchema();
    if (relation->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }
    cout << schema->getName() << " has " << relation->getNumColumns() << " columns, "
         << relation->getNumRows() << " rows, and a maximum NIP of " << setw(2)
         << relation->getMaximumNip() << "." << endl;

    auto startTime = std::chrono::system_clock::now();
    //Agree sets

    AgreeSetFactory agreeSetFactory = AgreeSetFactory(relation.get());
    std::set<Vertical> agreeSets = unorderedToOrdered(agreeSetFactory.genAgreeSets());

    for(auto ag : agreeSets){
        std::cerr << ag.toString() << "\n";
    }

    //maximal sets
    CMAXGen cmaxSet = CMAXGen(schema);
    cmaxSet.execute(agreeSets);

    for(auto set : cmaxSet.getCmaxSets()){
        std::cerr << set.getColumn().toString();
        for(auto comb : set.getCombinations()){
            std::cerr << "\t\t" << comb.toString() << "\n";
        }
        std::cerr << "\n";
    }

    std::cerr << "-----------MAXSETS::::::" << std::endl;

    for(auto set : cmaxSet.getMaxSets()){
        std::cerr << set.getColumn().toString();
        for(auto comb : set.getCombinations()){
            std::cerr << "\t\t" << comb.toString() << "\n";
        }
        std::cerr << "\n";
    }

    //LHS

    for(auto& attribute : schema->getColumns()){
        std::set<Vertical> li;
        CMAXSet correct = genFirstLi(cmaxSet.getCmaxSets(), *attribute, li);
        if(schema->getNumColumns() != 1
            && correct.getCombinations().begin().operator*().getColumns().size()
             == schema->getColumns().size()){
            this->registerFD(Vertical(), *attribute);
            continue;
        }
        while(!li.empty()){
            for(Vertical l : li){
                bool isFD = true;
                for(auto combination : correct.getCombinations()){
                    if(!l.intersects(combination)){
                        isFD = false;
                        break;
                    }
                }
                if(isFD){
                    if(!l.contains(*attribute)){
                        this->registerFD(l, *attribute);
                    }
                    li.erase(l);
                }
                if(li.size() == 0){
                    break;
                }
            }
            li = genNextLi(li);
        }
    }

    cout << "TOTAL FD COUNT: " << this->fdCollection_.size() << "\n";

    // for(auto const& fd : this->fdCollection_){
    //     std::cout << "Registered FD: " << fd.getLhs().toString() << "->" << fd.getRhs().toString() << "\n";
    // }
    std::chrono::milliseconds elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    return elapsed_milliseconds.count();
}

CMAXSet Depminer::genFirstLi(std::set<CMAXSet> cmaxSets, Column attribute, std::set<Vertical> & li){
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

bool checkJoin(Vertical const& _p, Vertical const& _q);

//Apriori-gen function
std::set<Vertical> Depminer::genNextLi(std::set<Vertical> const& li){
    std::set<Vertical> ck;
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
    std::set<Vertical> result;
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
    return result;
}

bool checkJoin(Vertical const& _p, Vertical const& _q){
    dynamic_bitset<> p = _p.getColumnIndices();
    dynamic_bitset<> q = _q.getColumnIndices();

    int pPrev = -1, qPrev = -1;

    for(int i = p.size(); i >= 0; i--){
        pPrev = std::max(pPrev, p[i] ? i : -1);
        qPrev = std::max(qPrev, q[i] ? i : -1);
    }
    if(pPrev >= qPrev) return false;
    for(int i = 0; i < pPrev; i++){
        if(p[i] != q[i]) return false;
    }
    
    return true;
}