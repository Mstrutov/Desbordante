#include "CMAXSet.h"

class MAXSet : public CMAXSet{
public:
    MAXSet(Column column) : CMAXSet(column){};
    MAXSet() = default;
    bool operator<(MAXSet const& rhs) const{
        return this->getColumn() < rhs.getColumn();
    }
    void makeNewCombinations(std::set<Vertical> comb){
        this->columnCombinations = comb;
    }
};