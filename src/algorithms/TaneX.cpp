#include "TaneX.h"

#include <chrono>
#include <iostream>
#include <iomanip>
#include <list>
#include <memory>

#include "ColumnCombination.h"
#include "ColumnData.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "LatticeLevel.h"
#include "LatticeVertex.h"

double Tane::calculateZeroAryFdError(ColumnData const* rhs, ColumnLayoutRelationData const* relationData) {
    return 1 - rhs->getPositionListIndex()->getNepAsLong() / static_cast<double>(relationData->getNumTuplePairs());
}

double Tane::calculateFdError(PositionListIndex const* lhsPli, PositionListIndex const* jointPli,
                              ColumnLayoutRelationData const* relationData) {
    return (double) (lhsPli->getNepAsLong() - jointPli->getNepAsLong()) / static_cast<double>(relationData->getNumTuplePairs());
}


double Tane::calculateUccError(PositionListIndex const* pli, ColumnLayoutRelationData const* relationData) {
    return pli->getNepAsLong() / static_cast<double>(relationData->getNumTuplePairs());
}

void Tane::registerFD(Vertical const& lhs, Column const* rhs, double error, RelationalSchema const* schema) {
    dynamic_bitset<> lhs_bitset = lhs.getColumnIndices();
    /*std::cout << "Discovered FD: ";
    for (size_t i = lhs_bitset.find_first(); i != dynamic_bitset<>::npos; i = lhs_bitset.find_next(i)) {
        std::cout << schema->getColumn(i)->getName() << " ";
    }
    std::cout << "-> " << rhs->getName() << " - error equals " << error << std::endl;*/
    FDAlgorithm::registerFD(lhs, *rhs);
    countOfFD++;
}

/*void Tane::registerFD(Vertical const* lhs, Column const* rhs, double error, RelationalSchema const* schema) {
    registerFD(*lhs, rhs, error, schema);
}*/

void Tane::registerUCC(Vertical const& key, double error, RelationalSchema const* schema)  {
    /*dynamic_bitset<> key_bitset = key.getColumnIndices();
    std::cout << "Discovered UCC: ";
    for (int i = key_bitset.find_first(); i != -1; i = key_bitset.find_next(i)) {
        std::cout << schema->getColumn(i)->getName() << " ";
    }
    std::cout << "- error equals " << error << std::endl;*/
    countOfUCC++;
}



unsigned long long Tane::executeInternal() {
    RelationalSchema const* schema = relation_->getSchema();

    std::cout << schema->getName() << " has " << relation_->getNumColumns() << " columns, "
         << relation_->getNumRows() << " rows, and a maximum NIP of " << std::setw(2)
         << relation_->getMaximumNip() << "." << std::endl;

    for (auto& column : schema->getColumns()) {
        double avgPartners = relation_->getColumnData(column->getIndex()).
                getPositionListIndex()->getNepAsLong() * 2.0 / relation_->getNumRows();
        std::cout << "* " << column->toString() << ": every tuple has " << std::setw(2)
             << avgPartners << " partners on average." << std::endl;
    }
    auto startTime = std::chrono::system_clock::now();
    double progressStep = 100.0 / (schema->getNumColumns() + 1);

    //Initialize level 0
    std::vector<std::unique_ptr<LatticeLevel>> levels;
    std::unique_ptr<LatticeLevel> level0 = std::make_unique<LatticeLevel>(0);
    // TODO: через указатели кажется надо переделать
    level0->add(std::make_unique<LatticeVertex>(*(schema->emptyVertical)));
    LatticeVertex const* emptyVertex = level0->getVertices().begin()->second.get();
    levels.push_back(std::move(level0));
    addProgress(progressStep);

    //Initialize level1
    dynamic_bitset<> zeroaryFdRhs(schema->getNumColumns());
    std::unique_ptr<LatticeLevel> level1 = std::make_unique<LatticeLevel>(1);
    for (auto& column : schema->getColumns()) {
        //for each attribute set vertex
        ColumnData const& columnData = relation_->getColumnData(column->getIndex());
        auto vertex = std::make_unique<LatticeVertex>(static_cast<Vertical>(*column));

        vertex->addRhsCandidates(schema->getColumns());
        vertex->getParents().push_back(emptyVertex);
        vertex->setKeyCandidate(true);
        vertex->setPositionListIndex(columnData.getPositionListIndex());

        //check FDs: 0->A
        double fdError = calculateZeroAryFdError(&columnData, relation_.get());
        if (fdError <= maxFdError) {  //TODO: max_error
            zeroaryFdRhs.set(column->getIndex());
            registerFD(*schema->emptyVertical, column.get(), fdError, schema);

            vertex->getRhsCandidates().set(column->getIndex(), false);
            if (fdError == 0) {
                vertex->getRhsCandidates().reset();
            }
        }

        level1->add(std::move(vertex));
    }

    for (auto& [key_map, vertex] : level1->getVertices()) {
        Vertical column = vertex->getVertical();
        vertex->getRhsCandidates() &= ~zeroaryFdRhs;  //~ returns flipped copy <- removed already discovered zeroary FDs

        // вот тут костыль, чтобы вытянуть индекс колонки из вершины, в которой только один индекс
        ColumnData const& columnData = relation_->getColumnData(column.getColumnIndices().find_first());
        double uccError = calculateUccError(columnData.getPositionListIndex(), relation_.get());
        if (uccError <= maxUccError) {
            registerUCC(column, uccError, schema);
            vertex->setKeyCandidate(false);
            if (uccError == 0) {
                for (unsigned long rhsIndex = vertex->getRhsCandidates().find_first();
                     rhsIndex < vertex->getRhsCandidates().size();
                     rhsIndex = vertex->getRhsCandidates().find_next(rhsIndex)){
                    if (rhsIndex != column.getColumnIndices().find_first()){
                        registerFD(column, schema->getColumn(rhsIndex), 0, schema);
                    }
                }
                vertex->getRhsCandidates() &= column.getColumnIndices();
                //set vertex invalid if we seek for exact dependencies
                if (maxFdError == 0 && maxUccError == 0) {
                    vertex->setInvalid(true);
                }
            }
        }
    }
    levels.push_back(std::move(level1));
    addProgress(progressStep);

    for (unsigned int arity = 2; arity <= maxArity; arity++) {
        //auto startTime = std::chrono::system_clock::now();
        LatticeLevel::clearLevelsBelow(levels, arity - 1);
        LatticeLevel::generateNextLevel(levels);
        //std::chrono::duration<double> elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
        //aprioriMillis += elapsed_milliseconds.count();

        LatticeLevel* level = levels[arity].get();
        std::cout << "Checking " << level->getVertices().size() << " " << arity << "-ary lattice vertices." << std::endl;
        if (level->getVertices().empty()) {
            break;
        }

        for (auto& [key_map, xaVertex] : level->getVertices()) {
            if (xaVertex->getIsInvalid()){
                continue;
            }

            Vertical xa = xaVertex->getVertical();
            //Calculate XA PLI
            if (xaVertex->getPositionListIndex() == nullptr) {
                auto parentPLI1 = xaVertex->getParents()[0]->getPositionListIndex();
                auto parentPLI2 = xaVertex->getParents()[1]->getPositionListIndex();
                xaVertex->acquirePositionListIndex(parentPLI1->intersect(parentPLI2));
            }

            dynamic_bitset<> xaIndices = xa.getColumnIndices();
            dynamic_bitset<> aCandidates = xaVertex->getRhsCandidates();

            for (const auto& xVertex : xaVertex->getParents()) {
                Vertical const& lhs = xVertex->getVertical();

                // Find index of A in XA. If a is not a candidate, continue. TODO: possible to do it easier??
                //like "aIndex = xaIndices - xIndices;"
                int aIndex = xaIndices.find_first();
                dynamic_bitset<> xIndices = lhs.getColumnIndices();
                while (aIndex >= 0 && xIndices[aIndex]) {
                    aIndex = xaIndices.find_next(aIndex);
                }
                if (!aCandidates[aIndex]) {
                    continue;
                }

                // Check X -> A
                double error = calculateFdError(
                        xVertex->getPositionListIndex(),
                        xaVertex->getPositionListIndex(),
                        relation_.get());
                if (error <= maxFdError) {
                    Column const* rhs = schema->getColumns()[aIndex].get();

                    //TODO: register FD to a file or something
                    registerFD(lhs, rhs, error, schema);
                    xaVertex->getRhsCandidates().set(rhs->getIndex(), false);
                    if (error == 0) {
                        xaVertex->getRhsCandidates() &= lhs.getColumnIndices();
                    }
                }
            }
        }

        //Prune
        //cout << "Pruning level: " << level->getArity() << ". " << level->getVertices().size() << " vertices" << endl;
        std::list<LatticeVertex *> keyVertices;
        for (auto& [map_key, vertex] : level->getVertices()) {
            Vertical columns = vertex->getVertical();            //Originally it's a ColumnCombination

            if (vertex->getIsKeyCandidate()) {
                double uccError = calculateUccError(vertex->getPositionListIndex(), relation_.get());
                if (uccError <= maxUccError) {       //If a key candidate is an approx UCC
                    //TODO: do smth with UCC

                    registerUCC(columns, uccError, schema);
                    vertex->setKeyCandidate(false);
                    if (uccError == 0) {
                        for (size_t rhsIndex = vertex->getRhsCandidates().find_first();
                             rhsIndex != boost::dynamic_bitset<>::npos;
                             rhsIndex = vertex->getRhsCandidates().find_next(rhsIndex)) {
                            Vertical rhs = static_cast<Vertical>(*schema->getColumn((int)rhsIndex));
                            if (!columns.contains(rhs)) {
                                bool isRhsCandidate = true;
                                for (const auto& column : columns.getColumns()) {
                                    Vertical sibling = columns.without(static_cast<Vertical>(*column)).Union(rhs);
                                    auto siblingVertex = level->getLatticeVertex(sibling.getColumnIndices());
                                    if (siblingVertex == nullptr ||
                                         !siblingVertex->getConstRhsCandidates()[rhs.getColumnIndices().find_first()]) {
                                        isRhsCandidate = false;
                                        break;
                                    }
                                    // for each outer rhs: if there is a sibling s.t. it doesn't have this rhs, there is no FD: vertex->rhs
                                }
                                //Found fd: vertex->rhs => register it
                                if (isRhsCandidate){
                                    registerFD(columns, schema->getColumn(rhsIndex), 0, schema);
                                }
                            }
                        }
                        keyVertices.push_back(vertex.get());
                        //cout << "--------------------------" << endl << "KeyVert: " << *vertex;
                    }
                }
            }
            //if we seek for exact FDs then setInvalid
            if (maxFdError == 0 && maxUccError == 0){
                for (auto keyVertex : keyVertices){
                    keyVertex->getRhsCandidates() &= keyVertex->getVertical().getColumnIndices();
                    keyVertex->setInvalid(true);
                }
            }

        }

        //cout << "Pruned level: " << level->getArity() << ". " << level->getVertices().size() << " vertices" << endl;
        //TODO: printProfilingData
        addProgress(progressStep);
    }
    
    setProgress(100);
    std::chrono::milliseconds elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);
    aprioriMillis += elapsed_milliseconds.count();

    std::cout << "Time: " << aprioriMillis << " milliseconds" << std::endl;
    std::cout << "Intersection time: " << PositionListIndex::micros / 1000 << "ms" << std::endl;
    std::cout << "Total intersections: " << PositionListIndex::intersectionCount << std::endl;
    std::cout << "Total FD count: " << countOfFD << std::endl;
    std::cout << "Total UCC count: " << countOfUCC << std::endl;
    // std::cout << "===== FD JSON ========" << getJsonFDs() << std::endl;
    std::cout << "HASH: " << fletcher16() << std::endl;

    return aprioriMillis;
}
