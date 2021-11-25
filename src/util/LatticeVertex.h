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
    Vertical vertical;
    // holds either an owned PLI (unique_ptr) or a non-owned one (const*)
    std::variant<std::unique_ptr<PositionListIndex>, PositionListIndex const*> positionListIndex_;
    boost::dynamic_bitset<> rhsCandidates;
    bool isKeyCandidate = false;
    std::vector<LatticeVertex const*> parents;
    bool isInvalid = false;

public:
    explicit LatticeVertex(Vertical _vertical) : vertical(std::move(_vertical)),
                                                 rhsCandidates(vertical.getSchema()->getNumColumns()) {}

    std::vector<LatticeVertex const*>& getParents() { return parents; }

    Vertical const& getVertical() const { return vertical; }
    boost::dynamic_bitset<>& getRhsCandidates() { return rhsCandidates; }
    boost::dynamic_bitset<> const& getConstRhsCandidates() const { return rhsCandidates; }

    void addRhsCandidates(std::vector<std::unique_ptr<Column>> const& candidates);

    bool comesBeforeAndSharePrefixWith(LatticeVertex const& that) const;
    bool getIsKeyCandidate() const { return isKeyCandidate; }
    void setKeyCandidate(bool m_isKeyCandidate) { isKeyCandidate = m_isKeyCandidate; }
    bool getIsInvalid() const { return isInvalid; }
    void setInvalid(bool m_isInvalid) { isInvalid = m_isInvalid; }

    PositionListIndex const* getPositionListIndex() const;
    void setPositionListIndex(PositionListIndex const* positionListIndex)
    { positionListIndex_ = positionListIndex; }
    void acquirePositionListIndex(std::unique_ptr<PositionListIndex> positionListIndex)
    { positionListIndex_ = std::move(positionListIndex); }

    bool operator> (LatticeVertex const& that) const;

    std::string toString();

    static bool comparator(LatticeVertex * v1, LatticeVertex * v2) { return *v2 > *v1; }
    friend std::ostream& operator<<(std::ostream& os, LatticeVertex& lv);
};

} // namespace util

