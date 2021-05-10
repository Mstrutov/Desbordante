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
    this->MaxSetsGenerate(agreeSets);
    this->CMaxSetsGenerate();
}

void CMAXGen::MaxSetsGenerate(std::set<Vertical> agreeSets){
    for(std::shared_ptr<Column> column : this->schema->getColumns()){
        MAXSet result(column);
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
                if(superSet.getColumnIndices().is_subset_of(set.getColumnIndices())){
                    setsDelete.insert(superSet);
                }
                if(toAdd){
                    toAdd = !superSet.getColumnIndices().is_subset_of(set.getColumnIndices());
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
        this->maxSets.insert(result);
    }
}

void CMAXGen::CMaxSetsGenerate(){
    for(MAXSet maxSet : this->maxSets){
        CMAXSet result = CMAXSet(maxSet.getColumn());
        for(auto combination : maxSet.getCombinations()){
            result.addCombination(combination.invert());
        }
        cmaxSets.insert(result);
    }
}

