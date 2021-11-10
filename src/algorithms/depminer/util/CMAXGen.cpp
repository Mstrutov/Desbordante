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

void CMAXGen::execute(std::unordered_set<Vertical>& agreeSets){
    auto startTime = std::chrono::system_clock::now();

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
        this->cmaxSets.insert(result);
    }

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    std::cout << "TIME TO GENERATE CMAXSETS: " << elapsed_milliseconds.count() << std::endl;
    std::cout << "TOTAL CMAX SETS: " << this->cmaxSets.size() << "\n";
}

