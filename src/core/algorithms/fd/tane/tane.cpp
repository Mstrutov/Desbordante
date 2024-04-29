#include "tane.h"

#include <chrono>
#include <iomanip>
#include <list>
#include <memory>

#include <easylogging++.h>

#include "config/error/option.h"
#include "config/max_lhs/option.h"
#include "lattice_level.h"
#include "lattice_vertex.h"
#include "model/table/column_data.h"
#include "model/table/column_layout_relation_data.h"
#include "model/table/relational_schema.h"

namespace algos {

using boost::dynamic_bitset;

Tane::Tane(std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : PliBasedFDAlgorithm({kDefaultPhaseName}, relation_manager) {
    RegisterOptions();
}

void Tane::RegisterOptions() {
    RegisterOption(config::kErrorOpt(&max_ucc_error_));
}

void Tane::MakeExecuteOptsAvailableFDInternal() {
    MakeOptionsAvailable({config::kErrorOpt.GetName()});
}

void Tane::ResetStateFd() {
    count_of_fd_ = 0;
    count_of_ucc_ = 0;
    apriori_millis_ = 0;
}

double Tane::CalculateZeroAryFdError(ColumnData const* rhs,
                                     ColumnLayoutRelationData const* relation_data) {
    return 1 - rhs->GetPositionListIndex()->GetNepAsLong() /
                       static_cast<double>(relation_data->GetNumTuplePairs());
}

double Tane::CalculateFdError(model::PositionListIndex const* lhs_pli,
                              model::PositionListIndex const* joint_pli,
                              ColumnLayoutRelationData const* relation_data) {
    return (double)(lhs_pli->GetNepAsLong() - joint_pli->GetNepAsLong()) /
           static_cast<double>(relation_data->GetNumTuplePairs());
}

double Tane::CalculateUccError(model::PositionListIndex const* pli,
                               ColumnLayoutRelationData const* relation_data) {
    return pli->GetNepAsLong() / static_cast<double>(relation_data->GetNumTuplePairs());
}

void Tane::RegisterAndCountFd(Vertical const& lhs, Column const& rhs, [[maybe_unused]] double error,
                              [[maybe_unused]] RelationalSchema const* schema) {
    dynamic_bitset<> lhs_bitset = lhs.GetColumnIndices();
    PliBasedFDAlgorithm::RegisterFd(lhs, rhs);
    count_of_fd_++;
}

void Tane::RegisterUcc([[maybe_unused]] Vertical const& key, [[maybe_unused]] double error,
                       [[maybe_unused]] RelationalSchema const* schema) {
    /*dynamic_bitset<> key_bitset = key.getColumnIndices();
    LOG(INFO) << "Discovered UCC: ";
    for (int i = key_bitset.find_first(); i != -1; i = key_bitset.find_next(i)) {
        LOG(INFO) << schema->GetColumn(i)->GetName() << " ";
    }
    LOG(INFO) << "- error equals " << error << std::endl;*/
    count_of_ucc_++;
}

unsigned long long Tane::ExecuteInternal() {
    max_fd_error_ = max_ucc_error_;
    RelationalSchema const* schema = relation_->GetSchema();

    LOG(INFO) << schema->GetName() << " has " << relation_->GetNumColumns() << " columns, "
              << relation_->GetNumRows() << " rows, and a maximum NIP of " << std::setw(2)
              << relation_->GetMaximumNip() << ".";

    for (auto const& column : schema->GetColumns()) {
        double avg_partners =
                relation_->GetColumnData(column.GetIndex()).GetPositionListIndex()->GetNepAsLong() *
                2.0 / relation_->GetNumRows();
        LOG(INFO) << "* " << column.ToString() << ": every tuple has " << std::setw(2)
                  << avg_partners << " partners on average.";
    }
    auto start_time = std::chrono::system_clock::now();
    double progress_step = 100.0 / (schema->GetNumColumns() + 1);

    // Initialize level 0
    std::vector<std::unique_ptr<model::LatticeLevel>> levels;
    auto level0 = std::make_unique<model::LatticeLevel>(0);
    // TODO: через указатели кажется надо переделать
    level0->Add(std::make_unique<model::LatticeVertex>(*(schema->empty_vertical_)));
    model::LatticeVertex const* empty_vertex = level0->GetVertices().begin()->second.get();
    levels.push_back(std::move(level0));
    AddProgress(progress_step);

    // Initialize level1
    dynamic_bitset<> zeroary_fd_rhs(schema->GetNumColumns());
    auto level1 = std::make_unique<model::LatticeLevel>(1);
    for (auto const& column : schema->GetColumns()) {
        // for each attribute set vertex
        ColumnData const& column_data = relation_->GetColumnData(column.GetIndex());
        auto vertex = std::make_unique<model::LatticeVertex>(static_cast<Vertical>(column));

        vertex->AddRhsCandidates(schema->GetColumns());
        vertex->GetParents().push_back(empty_vertex);
        vertex->SetKeyCandidate(true);
        vertex->SetPositionListIndex(column_data.GetPositionListIndex());

        // check FDs: 0->A
        double fd_error = CalculateZeroAryFdError(&column_data, relation_.get());
        if (fd_error <= max_fd_error_) {  // TODO: max_error
            zeroary_fd_rhs.set(column.GetIndex());
            RegisterAndCountFd(*schema->empty_vertical_, column, fd_error, schema);

            vertex->GetRhsCandidates().set(column.GetIndex(), false);
            if (fd_error == 0) {
                vertex->GetRhsCandidates().reset();
            }
        }

        level1->Add(std::move(vertex));
    }

    for (auto& [key_map, vertex] : level1->GetVertices()) {
        Vertical column = vertex->GetVertical();
        vertex->GetRhsCandidates() &=
                ~zeroary_fd_rhs;  //~ returns flipped copy <- removed already discovered zeroary FDs

        // вот тут костыль, чтобы вытянуть индекс колонки из вершины, в которой только один индекс
        ColumnData const& column_data =
                relation_->GetColumnData(column.GetColumnIndices().find_first());
        double ucc_error = CalculateUccError(column_data.GetPositionListIndex(), relation_.get());
        if (ucc_error <= max_ucc_error_) {
            RegisterUcc(column, ucc_error, schema);
            vertex->SetKeyCandidate(false);
            if (ucc_error == 0 && max_lhs_ != 0) {
                for (unsigned long rhs_index = vertex->GetRhsCandidates().find_first();
                     rhs_index < vertex->GetRhsCandidates().size();
                     rhs_index = vertex->GetRhsCandidates().find_next(rhs_index)) {
                    if (rhs_index != column.GetColumnIndices().find_first()) {
                        RegisterAndCountFd(column, schema->GetColumn(rhs_index), 0, schema);
                    }
                }
                vertex->GetRhsCandidates() &= column.GetColumnIndices();
                // set vertex invalid if we seek for exact dependencies
                if (max_fd_error_ == 0 && max_ucc_error_ == 0) {
                    vertex->SetInvalid(true);
                }
            }
        }
    }
    levels.push_back(std::move(level1));
    AddProgress(progress_step);

    unsigned int max_arity =
            max_lhs_ == std::numeric_limits<unsigned int>::max() ? max_lhs_ : max_lhs_ + 1;
    for (unsigned int arity = 2; arity <= max_arity; arity++) {
        // auto start_time = std::chrono::system_clock::now();
        model::LatticeLevel::ClearLevelsBelow(levels, arity - 1);
        model::LatticeLevel::GenerateNextLevel(levels);
        // std::chrono::duration<double> elapsed_milliseconds =
        // std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() -
        // start_time); apriori_millis_ += elapsed_milliseconds.count();

        model::LatticeLevel* level = levels[arity].get();
        LOG(TRACE) << "Checking " << level->GetVertices().size() << " " << arity
                   << "-ary lattice vertices.";
        if (level->GetVertices().empty()) {
            break;
        }

        for (auto& [key_map, xa_vertex] : level->GetVertices()) {
            if (xa_vertex->GetIsInvalid()) {
                continue;
            }

            Vertical xa = xa_vertex->GetVertical();
            // Calculate XA PLI
            if (xa_vertex->GetPositionListIndex() == nullptr) {
                auto parent_pli_1 = xa_vertex->GetParents()[0]->GetPositionListIndex();
                auto parent_pli_2 = xa_vertex->GetParents()[1]->GetPositionListIndex();
                xa_vertex->AcquirePositionListIndex(parent_pli_1->Intersect(parent_pli_2));
            }

            dynamic_bitset<> xa_indices = xa.GetColumnIndices();
            dynamic_bitset<> a_candidates = xa_vertex->GetRhsCandidates();

            for (auto const& x_vertex : xa_vertex->GetParents()) {
                Vertical const& lhs = x_vertex->GetVertical();

                // Find index of A in XA. If a is not a candidate, continue. TODO: possible to do it
                // easier??
                // like "a_index = xa_indices - x_indices;"
                int a_index = xa_indices.find_first();
                dynamic_bitset<> x_indices = lhs.GetColumnIndices();
                while (a_index >= 0 && x_indices[a_index]) {
                    a_index = xa_indices.find_next(a_index);
                }
                if (!a_candidates[a_index]) {
                    continue;
                }

                // Check X -> A
                double error = CalculateFdError(x_vertex->GetPositionListIndex(),
                                                xa_vertex->GetPositionListIndex(), relation_.get());
                if (error <= max_fd_error_) {
                    Column const& rhs = schema->GetColumns()[a_index];

                    // TODO: register FD to a file or something
                    RegisterAndCountFd(lhs, rhs, error, schema);
                    xa_vertex->GetRhsCandidates().set(rhs.GetIndex(), false);
                    if (error == 0) {
                        xa_vertex->GetRhsCandidates() &= lhs.GetColumnIndices();
                    }
                }
            }
        }

        if (arity == max_arity) {
            break;
        }

        // Prune
        // cout << "Pruning level: " << level->GetArity() << ". " << level->GetVertices().size() <<
        // " vertices_" << endl;
        std::list<model::LatticeVertex*> key_vertices;
        for (auto& [map_key, vertex] : level->GetVertices()) {
            Vertical columns = vertex->GetVertical();  // Originally it's a ColumnCombination

            if (vertex->GetIsKeyCandidate()) {
                double ucc_error =
                        CalculateUccError(vertex->GetPositionListIndex(), relation_.get());
                if (ucc_error <= max_ucc_error_) {  // If a key candidate is an approx UCC
                    // TODO: do smth with UCC

                    RegisterUcc(columns, ucc_error, schema);
                    vertex->SetKeyCandidate(false);
                    if (ucc_error == 0) {
                        for (size_t rhs_index = vertex->GetRhsCandidates().find_first();
                             rhs_index != boost::dynamic_bitset<>::npos;
                             rhs_index = vertex->GetRhsCandidates().find_next(rhs_index)) {
                            auto rhs = static_cast<Vertical>(schema->GetColumn(rhs_index));
                            if (!columns.Contains(rhs)) {
                                bool is_rhs_candidate = true;
                                for (auto const& column : columns.GetColumns()) {
                                    Vertical sibling =
                                            columns.Without(static_cast<Vertical>(*column))
                                                    .Union(rhs);
                                    auto sibling_vertex =
                                            level->GetLatticeVertex(sibling.GetColumnIndices());
                                    if (sibling_vertex == nullptr ||
                                        !sibling_vertex->GetConstRhsCandidates()
                                                 [rhs.GetColumnIndices().find_first()]) {
                                        is_rhs_candidate = false;
                                        break;
                                    }
                                    // for each outer rhs: if there is a sibling s.t. it doesn't
                                    // have this rhs, there is no FD: vertex->rhs
                                }
                                // Found fd: vertex->rhs => register it
                                if (is_rhs_candidate) {
                                    RegisterAndCountFd(columns, schema->GetColumn(rhs_index), 0,
                                                       schema);
                                }
                            }
                        }
                        key_vertices.push_back(vertex.get());
                        // cout << "--------------------------" << endl << "KeyVert: " << *vertex;
                    }
                }
            }
            // if we seek for exact FDs then SetInvalid
            if (max_fd_error_ == 0 && max_ucc_error_ == 0) {
                for (auto key_vertex : key_vertices) {
                    key_vertex->GetRhsCandidates() &= key_vertex->GetVertical().GetColumnIndices();
                    key_vertex->SetInvalid(true);
                }
            }
        }

        // TODO: printProfilingData
        AddProgress(progress_step);
    }

    SetProgress(100);
    std::chrono::milliseconds elapsed_milliseconds =
            std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() -
                                                                  start_time);
    apriori_millis_ += elapsed_milliseconds.count();

    LOG(INFO) << "Time: " << apriori_millis_ << " milliseconds";
    LOG(INFO) << "Intersection time: " << model::PositionListIndex::micros_ / 1000 << "ms";
    LOG(INFO) << "Total intersections: " << model::PositionListIndex::intersection_count_
              << std::endl;
    LOG(INFO) << "Total FD count: " << count_of_fd_;
    LOG(INFO) << "Total UCC count: " << count_of_ucc_;
    LOG(INFO) << "HASH: " << Fletcher16();

    return apriori_millis_;
}

}  // namespace algos
