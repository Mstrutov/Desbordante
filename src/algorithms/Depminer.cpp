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

using boost::dynamic_bitset, std::make_shared, std::shared_ptr, std::cout, std::endl, std::setw, std::vector, std::list, std::dynamic_pointer_cast;

unsigned long long Depminer::execute(){
    auto startTime = std::chrono::system_clock::now();
    shared_ptr<ColumnLayoutRelationData> relation = ColumnLayoutRelationData::createFrom(inputGenerator_, true);
    shared_ptr<RelationalSchema> schema = relation->getSchema();
    if (relation->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }
    cout << schema->getName() << " has " << relation->getNumColumns() << " columns, "
         << relation->getNumRows() << " rows, and a maximum NIP of " << setw(2)
         << relation->getMaximumNip() << "." << endl;

    //Agree sets
    std::set<Vertical> agreeSets;

    //maximal sets
    CMAXGen cmaxSet = CMAXGen(schema);
    cmaxSet.execute(agreeSets);

    //LHS



    return 0;
}