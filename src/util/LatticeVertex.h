#pragma once

#include <list>
#include <utility>
#include <variant>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "PositionListIndex.h"
#include "RelationalSchema.h"
#include "Vertical.h"

namespace util {

class LatticeVertex{
private:
    Vertical vertical_;
    // holds either an owned PLI (unique_ptr) or a non-owned one (const*)
    std::variant<std::unique_ptr<PositionListIndex>, PositionListIndex const*> position_list_index_;
    boost::dynamic_bitset<> rhs_candidates_;
    bool is_key_candidate_ = false;
    std::vector<LatticeVertex const*> parents_;
    bool is_invalid_ = false;

public:
    explicit LatticeVertex(Vertical vertical) : vertical_(std::move(vertical)),
                                                rhs_candidates_(vertical_.GetSchema()->GetNumColumns()) {}

    std::vector<LatticeVertex const*>& GetParents() { return parents_; }

    Vertical const& GetVertical() const { return vertical_; }
    boost::dynamic_bitset<>& GetRhsCandidates() { return rhs_candidates_; }
    boost::dynamic_bitset<> const& GetConstRhsCandidates() const { return rhs_candidates_; }

    void AddRhsCandidates(std::vector<std::unique_ptr<Column>> const& candidates);

    bool ComesBeforeAndSharePrefixWith(LatticeVertex const& that) const;
    bool GetIsKeyCandidate() const { return is_key_candidate_; }
    void SetKeyCandidate(bool m_is_key_candidate) { is_key_candidate_ = m_is_key_candidate; }
    bool GetIsInvalid() const { return is_invalid_; }
    void SetInvalid(bool m_is_invalid) { is_invalid_ = m_is_invalid; }

    PositionListIndex const* GetPositionListIndex() const;
    void SetPositionListIndex(PositionListIndex const* position_list_index)
    { position_list_index_ = position_list_index; }
    void AcquirePositionListIndex(std::unique_ptr<PositionListIndex> position_list_index)
    { position_list_index_ = std::move(position_list_index); }

    bool operator> (LatticeVertex const& that) const;

    std::string ToString();

    static bool Comparator(LatticeVertex * v1, LatticeVertex * v2) { return *v2 > *v1; }
    friend std::ostream& operator<<(std::ostream& os, LatticeVertex& lv);
};

} // namespace util

