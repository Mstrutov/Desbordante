#include <chrono>
#include <iostream>
#include <iomanip>
#include <list>
#include <memory>

#include "ColumnCombination.h"
#include "ColumnData.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "CMAXGen.h"

using std::cerr;

void CMAXGen::execute(std::set<Vertical> agreeSets){
    // auto startTime = std::chrono::system_clock::now();
    this->MaxSetsGenerate(agreeSets);
    // std::chrono::milliseconds elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    // std::cout << "TIME TO GENERATE MAXSETS: " << elapsed_milliseconds.count() << std::endl;
    // auto newStartTime = std::chrono::system_clock::now();
    this->CMaxSetsGenerate();
    // elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - newStartTime);
    // std::cout << "TIME TO GENERATE CMAXSETS: " << elapsed_milliseconds.count() << std::endl;
    std::cout << "TOTAL CMAX SETS: " << this->cmaxSets.size() << "\n";
}

void CMAXGen::MaxSetsGenerate(std::set<Vertical> agreeSets){
    for(auto& column : this->schema->getColumns()){
        MAXSet result(*column);
        for(Vertical ag : agreeSets){
            if(ag.contains(*column)){
                continue;
            }
            result.addCombination(ag);
        }
        //in MAXSet need to add only maximal sets

        std::set<Vertical> superSets;
        std::set<Vertical> setsDelete;
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
        result.makeNewCombinations(superSets);
        this->maxSets.insert(result);
    }
}

void CMAXGen::CMaxSetsGenerate(){
    for(MAXSet maxSet : this->maxSets){
        CMAXSet result = CMAXSet(maxSet.getColumn());
        for(auto combination : maxSet.getCombinations()){
            result.addCombination(combination.invert());
        }
        this->cmaxSets.insert(result);
    }
}

