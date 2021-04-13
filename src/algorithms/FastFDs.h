#pragma once

#include <functional>

#include "Vertical.h"
#include "FDAlgorithm.h"
#include "ColumnLayoutRelationData.h"

class FastFDs : public FDAlgorithm {
public:
    using FDAlgorithm::FDAlgorithm;
    unsigned long long execute() override;
private:
    using OrderingComparator = std::function<bool (Column const&, Column const&)>;

    /* Computes all difference sets of `relation_`
     * using simplified algorithm proposed by Wyss et al.
     * TODO: consider using an algorithm with identifier
     *       sets to compute the agree sets. 
     */
    void genDiffSets();

    /* Computes minimal difference sets
     * of `relation_` modulo `col`
     * TODO: Parameter `col` must be Column const&
     * but const forbids construct Vertical from Column
     */
    std::vector<std::shared_ptr<Vertical>> getDiffSetsMod(Column const& col) const;
    std::shared_ptr<Vertical> getAgreeSet(std::vector<int> const& tuple1,
                                          std::vector<int> const& tuple2) const;
    // getOrdering for initial cas,
    // TODO: rename to getInitOrdering
    std::set<Column, OrderingComparator> getOrdering(std::vector<std::shared_ptr<Vertical>> const& diff_sets,
                                                     Column const& attribute) const;
    std::set<Column, OrderingComparator> getOrdering(std::vector<std::shared_ptr<Vertical>> const& diff_sets,
                                                     Column const& attribute,
                                                     std::set<Column, OrderingComparator> const& cur_ordering) const;
    std::set<std::vector<int>> getPLIMaxRepresentation() const;

    void findCovers(Column const& attribute, std::vector<std::shared_ptr<Vertical>> const& diff_sets_mod,
                    std::vector<std::shared_ptr<Vertical>> const& cur_diff_sets, Vertical const& path,
                    std::set<Column, OrderingComparator> const& ordering);
    /* Returns true if `cover` is minimal cover of `diff_sets_mod`
     * false otherwise.
     */
    bool coverMinimal(Vertical const& cover, std::vector<std::shared_ptr<Vertical>> const& diff_sets_mod) const;
    bool isCover(Vertical const& candidate, std::vector<std::shared_ptr<Vertical>> const& sets) const;
    // static non member function?
    void calculateSupersets(std::set<std::vector<int>>& max_representation,
                            std::deque<std::vector<int>> partition) const;
    // need a better name and a description
    void verifyFDs();

    std::shared_ptr<ColumnLayoutRelationData> relation_;
    shared_ptr<RelationalSchema> schema_;
    std::vector<std::shared_ptr<Vertical>> diff_sets_;
};

