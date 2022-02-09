#include <algorithm>
#include <chrono>
#include <iostream>

#include "ColumnCombination.h"
#include "LatticeLevel.h"

namespace util {

using std::move, std::min, std::shared_ptr, std::vector, std::cout, std::endl, std::sort, std::make_shared;

void LatticeLevel::Add(std::unique_ptr<LatticeVertex> vertex) {
    vertices_.emplace(vertex->GetVertical().GetColumnIndices(), std::move(vertex));
}

LatticeVertex const* LatticeLevel::GetLatticeVertex(const boost::dynamic_bitset<>& column_indices) const {
    auto it = vertices_.find(column_indices);
    if (it != vertices_.end()) {
        return it->second.get();
    } else {
        return nullptr;
    }
}

void LatticeLevel::GenerateNextLevel(std::vector<std::unique_ptr<LatticeLevel>>& levels) {
    unsigned int arity = levels.size() - 1;
    assert(arity >= 1);
    std::cout << "-------------Creating level " << arity + 1 << "...-----------------\n";

    LatticeLevel* current_level = levels[arity].get();

    std::vector<LatticeVertex*> current_level_vertices;
    for (const auto& [map_key, vertex] : current_level->GetVertices()) {
        current_level_vertices.push_back(vertex.get());
    }

    std::sort(current_level_vertices.begin(), current_level_vertices.end(),
              LatticeVertex::Comparator);
    auto next_level = std::make_unique<LatticeLevel>(arity + 1);

    for (unsigned int vertex_index_1 = 0; vertex_index_1 < current_level_vertices.size();
         vertex_index_1++) {
        LatticeVertex* vertex1 = current_level_vertices[vertex_index_1];

        if (vertex1->GetRhsCandidates().none() && !vertex1->GetIsKeyCandidate()) {
            continue;
        }

        for (unsigned int vertex_index_2 = vertex_index_1 + 1;
             vertex_index_2 < current_level_vertices.size(); vertex_index_2++) {
            LatticeVertex* vertex2 = current_level_vertices[vertex_index_2];

            if (!vertex1->ComesBeforeAndSharePrefixWith(*vertex2)) {
                break;
            }

            if (!vertex1->GetRhsCandidates().intersects(vertex1->GetRhsCandidates()) &&
                !vertex2->GetIsKeyCandidate()) {
                continue;
            }

            Vertical child_columns = vertex1->GetVertical().Union(vertex2->GetVertical());
            std::unique_ptr<LatticeVertex> child_vertex =
                std::make_unique<LatticeVertex>(child_columns);

            dynamic_bitset<> parent_indices(vertex1->GetVertical().GetSchema()->GetNumColumns());
            parent_indices |= vertex1->GetVertical().GetColumnIndices();
            parent_indices |= vertex2->GetVertical().GetColumnIndices();

            child_vertex->GetRhsCandidates() |= vertex1->GetRhsCandidates();
            child_vertex->GetRhsCandidates() &= vertex2->GetRhsCandidates();
            child_vertex->SetKeyCandidate(vertex1->GetIsKeyCandidate() &&
                                          vertex2->GetIsKeyCandidate());
            child_vertex->SetInvalid(vertex1->GetIsInvalid() || vertex2->GetIsInvalid());

            for (unsigned int i = 0, skip_index = parent_indices.find_first(); i < arity - 1;
                 i++, skip_index = parent_indices.find_next(skip_index)) {
                parent_indices[skip_index] = false;
                LatticeVertex const* parent_vertex =
                    current_level->GetLatticeVertex(parent_indices);

                if (parent_vertex == nullptr) {
                    goto continueMidOuter;
                }
                child_vertex->GetRhsCandidates() &= parent_vertex->GetConstRhsCandidates();
                if (child_vertex->GetRhsCandidates().none()) {
                    goto continueMidOuter;
                }
                child_vertex->GetParents().push_back(parent_vertex);
                parent_indices[skip_index] = true;

                child_vertex->SetKeyCandidate(child_vertex->GetIsKeyCandidate() &&
                                              parent_vertex->GetIsKeyCandidate());
                child_vertex->SetInvalid(child_vertex->GetIsInvalid() ||
                                         parent_vertex->GetIsInvalid());

                if (!child_vertex->GetIsKeyCandidate() && child_vertex->GetRhsCandidates().none()) {
                    goto continueMidOuter;
                }
            }

            child_vertex->GetParents().push_back(vertex1);
            child_vertex->GetParents().push_back(vertex2);

            next_level->Add(std::move(child_vertex));

            continueMidOuter:
            continue;
            }
    }

    levels.push_back(std::move(next_level));
}

void LatticeLevel::ClearLevelsBelow(std::vector<std::unique_ptr<LatticeLevel>>& levels,
                                    unsigned int arity) {
    // Clear the levels from the level list
    auto it = levels.begin();

    for (unsigned int i = 0; i < std::min((unsigned int)levels.size(), arity); i++) {
        (*(it++))->GetVertices().clear();
    }

    //Clear child references
    if (arity < levels.size()) {
        for (auto& [map_key, retained_vertex] : levels[arity]->GetVertices()) {
            retained_vertex->GetParents().clear();
        }
    }
}

} // namespace util

