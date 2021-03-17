#include <algorithm>
#include <chrono>
#include <iostream>

#include "ColumnCombination.h"
#include "LatticeLevel.h"

using std::move, std::min, std::shared_ptr, std::vector, std::cout, std::endl, std::sort, std::make_shared;

void LatticeLevel::add(std::unique_ptr<LatticeVertex> vertex) {
    vertices.emplace(vertex->getVertical().getColumnIndices(), std::move(vertex));
}

LatticeVertex const* LatticeLevel::getLatticeVertex(const boost::dynamic_bitset<>& columnIndices) const {
    auto it = vertices.find(columnIndices);
    if (it != vertices.end()){
        return it->second.get();
    }
    else{
        return nullptr;
    }
}

void LatticeLevel::generateNextLevel(std::vector<std::unique_ptr<LatticeLevel>>& levels) {
    unsigned int arity = levels.size() - 1;
    assert(arity >= 1);
    std::cout << "-------------Creating level " << arity + 1 << "...-----------------\n";

    LatticeLevel* currentLevel = levels[arity].get();

    std::vector<LatticeVertex *> currentLevelVertices;
    for (const auto& [map_key, vertex] : currentLevel->getVertices()) {
        currentLevelVertices.push_back(vertex.get());
    }

    std::sort(currentLevelVertices.begin(), currentLevelVertices.end(), LatticeVertex::comparator);
    LatticeLevel nextLevel(arity + 1);

    for (unsigned int vertexIndex1 = 0; vertexIndex1 < currentLevelVertices.size(); vertexIndex1++){
        LatticeVertex* vertex1 = currentLevelVertices[vertexIndex1];

        if (vertex1->getRhsCandidates().none() && !vertex1->getIsKeyCandidate()) {
            continue;
        }

        for (unsigned int vertexIndex2 = vertexIndex1 + 1; vertexIndex2 < currentLevelVertices.size(); vertexIndex2++){
            LatticeVertex* vertex2 = currentLevelVertices[vertexIndex2];

            if (!vertex1->comesBeforeAndSharePrefixWith(*vertex2)){
                break;
            }

            if (!vertex1->getRhsCandidates().intersects(vertex1->getRhsCandidates()) && !vertex2->getIsKeyCandidate()){
                continue;
            }

            Vertical childColumns = vertex1->getVertical().Union(vertex2->getVertical());
            std::unique_ptr<LatticeVertex> childVertex = std::make_unique<LatticeVertex>(childColumns);

            dynamic_bitset<> parentIndices(vertex1->getVertical().getSchema()->getNumColumns());
            parentIndices |= vertex1->getVertical().getColumnIndices();
            parentIndices |= vertex2->getVertical().getColumnIndices();

            childVertex->getRhsCandidates() |= vertex1->getRhsCandidates();
            childVertex->getRhsCandidates() &= vertex2->getRhsCandidates();
            childVertex->setKeyCandidate(vertex1->getIsKeyCandidate() && vertex2->getIsKeyCandidate());
            childVertex->setInvalid(vertex1->getIsInvalid() || vertex2->getIsInvalid());

            for (unsigned int i = 0, skipIndex = parentIndices.find_first(); i < arity - 1; i++, skipIndex = parentIndices.find_next(skipIndex)) {
                parentIndices[skipIndex] = false;
                LatticeVertex const* parentVertex = currentLevel->getLatticeVertex(parentIndices);

                if (parentVertex == nullptr){
                    goto continueMidOuter;
                }
                childVertex->getRhsCandidates() &= parentVertex->getConstRhsCandidates();
                if (childVertex->getRhsCandidates().none()){
                    goto continueMidOuter;
                }
                childVertex->getParents().push_back(parentVertex);
                parentIndices[skipIndex] = true;

                childVertex->setKeyCandidate(childVertex->getIsKeyCandidate() && parentVertex->getIsKeyCandidate());
                childVertex->setInvalid(childVertex->getIsInvalid() || parentVertex->getIsInvalid());

                if (!childVertex->getIsKeyCandidate() && childVertex->getRhsCandidates().none()){
                    goto continueMidOuter;
                }
            }

            childVertex->getParents().push_back(vertex1);
            childVertex->getParents().push_back(vertex2);

            nextLevel.add(std::move(childVertex));

            continueMidOuter:
            continue;
        }
    }

    levels.push_back(std::make_unique<LatticeLevel>(nextLevel));
}

void LatticeLevel::clearLevelsBelow(std::vector<std::unique_ptr<LatticeLevel>>& levels, unsigned int arity) {
    // Clear the levels from the level list
    auto it = levels.begin();

    for (unsigned int i = 0; i < std::min((unsigned int)levels.size(), arity); i++) {
        (*(it++))->getVertices().clear();
    }

    //Clear child references
    if (arity < levels.size()) {
        for (auto& [map_key, retainedVertex] : levels[arity]->getVertices()){
            retainedVertex->getParents().clear();
        }
    }
}
