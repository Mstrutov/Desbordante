#include "CMAXSet.h"

class MAXSet : public CMAXSet{
public:
    MAXSet(std::shared_ptr<Column> column) : CMAXSet(column){};
    MAXSet() = default;
    bool operator<(MAXSet const& rhs) const{
        return *(this->getColumn()) < *(rhs.getColumn());
    }
};