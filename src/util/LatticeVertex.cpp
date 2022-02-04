#include "LatticeVertex.h"

namespace util {

using boost::dynamic_bitset, std::vector, std::shared_ptr, std::make_shared, std::string;

void LatticeVertex::AddRhsCandidates(vector<std::unique_ptr<Column>> const& candidates) {
    for (auto& cand_ptr : candidates){
        rhs_candidates_.set(cand_ptr->GetIndex());
    }
}

bool LatticeVertex::ComesBeforeAndSharePrefixWith(LatticeVertex const& that) const {
    dynamic_bitset<> this_indices = vertical_.GetColumnIndices();
    dynamic_bitset<> that_indices = that.vertical_.GetColumnIndices();

    int this_index = this_indices.find_first();
    int that_index = that_indices.find_first();

    int arity = this_indices.count();
    for (int i = 0; i < arity - 1; i++){
        if (this_index != that_index) return false;
        this_index = this_indices.find_next(this_index);
        that_index = that_indices.find_next(that_index);
    }

    return this_index < that_index;
}


bool LatticeVertex::operator> (LatticeVertex const& that) const {
    if (vertical_.GetArity() != that.vertical_.GetArity())
        return vertical_.GetArity() > that.vertical_.GetArity();

    dynamic_bitset this_indices = vertical_.GetColumnIndices();
    int this_index = this_indices.find_first();
    dynamic_bitset that_indices = that.vertical_.GetColumnIndices();
    int that_index = that_indices.find_first();

    int result;
    while (true){
        result = this_index - that_index;
        if (result)
            return (result > 0);
        this_index = this_indices.find_next(this_index);
        that_index = that_indices.find_next(that_index);
    }
}

string LatticeVertex::ToString() {
    return "Vtx" + vertical_.ToString();
}

std::ostream& operator<<(std::ostream& os, LatticeVertex& lv) {
    using std::endl;
    os << "Vertex: " << lv.vertical_.ToString() << endl;

    string rhs;
    for (size_t index = lv.rhs_candidates_.find_first();
         index != dynamic_bitset<>::npos;
         index = lv.rhs_candidates_.find_next(index)) {
        rhs += std::to_string(index) + " ";
    }
    os << "Rhs Candidates: " << rhs << endl;
    os << "IsKeyCandidate, IsInvalid: " << lv.is_key_candidate_ << ", " << lv.is_invalid_ << endl;

    os << endl;
    return os;
}

PositionListIndex const *LatticeVertex::GetPositionListIndex() const {
    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(position_list_index_)) {
        return std::get<std::unique_ptr<PositionListIndex>>(position_list_index_).get();
    } else {
        return std::get<PositionListIndex const *>(position_list_index_);
    }
}

} // namespace util

