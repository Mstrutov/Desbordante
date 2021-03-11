#include "LatticeVertex.h"

using boost::dynamic_bitset, std::vector, std::shared_ptr, std::make_shared, std::string;

void LatticeVertex::addRhsCandidates(vector<std::unique_ptr<Column>> const& candidates) {
    for (auto& candPtr : candidates){
        rhsCandidates.set(candPtr->getIndex());
    }
}

bool LatticeVertex::comesBeforeAndSharePrefixWith(LatticeVertex const& that) const {
    dynamic_bitset<> thisIndices = vertical.getColumnIndices();
    dynamic_bitset<> thatIndices = that.vertical.getColumnIndices();

    int thisIndex = thisIndices.find_first();
    int thatIndex = thatIndices.find_first();

    int arity = thisIndices.count();
    for (int i = 0; i < arity - 1; i++){
        if (thisIndex != thatIndex) return false;
        thisIndex = thisIndices.find_next(thisIndex);
        thatIndex = thatIndices.find_next(thatIndex);
    }

    return thisIndex < thatIndex;
}


bool LatticeVertex::operator> (LatticeVertex const& that) const {
    if (vertical.getArity() != that.vertical.getArity())
        return vertical.getArity() > that.vertical.getArity();

    dynamic_bitset thisIndices = vertical.getColumnIndices();
    int thisIndex = thisIndices.find_first();
    dynamic_bitset thatIndices = that.vertical.getColumnIndices();
    int thatIndex = thatIndices.find_first();

    int result;
    while (true){
        result = thisIndex - thatIndex;
        if (result)
            return (result > 0);
        thisIndex = thisIndices.find_next(thisIndex);
        thatIndex = thatIndices.find_next(thatIndex);
    }
}

string LatticeVertex::toString() {
    return "Vtx" + vertical.toString();
}

std::ostream& operator<<(std::ostream& os, LatticeVertex& lv) {
    using std::endl;
    os << "Vertex: " << lv.vertical.toString() << endl;

    string rhs;
    for (size_t index = lv.rhsCandidates.find_first();
         index != dynamic_bitset<>::npos;
         index = lv.rhsCandidates.find_next(index)) {
        rhs += std::to_string(index) + " ";
    }
    os << "Rhs Candidates: " << rhs << endl;
    os << "IsKeyCandidate, IsInvalid: " << lv.isKeyCandidate << ", " << lv.isInvalid << endl;

    os << endl;
    return os;
}

PositionListIndex const *LatticeVertex::getPositionListIndex() const {
    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(positionListIndex_)) {
        return std::get<std::unique_ptr<PositionListIndex>>(positionListIndex_).get();
    } else {
        return std::get<PositionListIndex const *>(positionListIndex_);
    }
}
