#include "FDTreeElement.h"
#include "boost/dynamic_bitset.hpp"

FDTreeElement::FDTreeElement(size_t maxAttributeNumber): maxAttributeNumber_(maxAttributeNumber){
    children_.resize(maxAttributeNumber);
}

bool FDTreeElement::checkFd(size_t i) const{
    return this->isFd_[i];
}

FDTreeElement* FDTreeElement::getChild(size_t i) const{
    return this->children_[i].get();
}

void FDTreeElement::addRhsAttribute(size_t i){
    this->rhsAttributes_.set(i);
}

const std::bitset<FDTreeElement::kMaxAttrNum>& FDTreeElement::getRhsAttributes() const{
    return this->rhsAttributes_;
}

void FDTreeElement::markAsLast(size_t i){
    this->isFd_.set(i);
}

bool FDTreeElement::isFinalNode(size_t attr_num) const{
    if (!this->rhsAttributes_[attr_num]){
        return false;
    }
    for (size_t attr = 0; attr < this->maxAttributeNumber_; ++attr){
        if (children_[attr] && children_[attr]->getRhsAttributes()[attr_num]){
            return false;
        }
    }
    return true;
}

bool FDTreeElement::containsGeneralization
(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num, size_t currentAttr) const
{
    if (this->isFd_[attr_num - 1]){
        return true;
    }
    
    size_t nextSetAttr = lhs._Find_next(currentAttr);
    if (nextSetAttr == kMaxAttrNum){
        return false;
    }
    bool found = false;
    if (this->children_[nextSetAttr - 1] && this->children_[nextSetAttr - 1]->getRhsAttributes()[attr_num]){
        found = this->children_[nextSetAttr - 1]->containsGeneralization(lhs, attr_num, nextSetAttr);
    }

    if (found){
        return true;
    }
    return this->containsGeneralization(lhs, attr_num, nextSetAttr);
}

bool FDTreeElement::getGeneralizationAndDelete
(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num, size_t currentAttr, std::bitset<kMaxAttrNum>& specLhs)
{
    if (this->isFd_[attr_num - 1]){
        this->isFd_.reset(attr_num - 1);
        this->rhsAttributes_.reset(attr_num);
        return true;
    }

    size_t nextSetAttr = lhs._Find_next(currentAttr);
    if (nextSetAttr == kMaxAttrNum){
        return false;
    }

    bool found = false;
    if (this->children_[nextSetAttr - 1] && this->children_[nextSetAttr - 1]->getRhsAttributes()[attr_num]){
        found = this->children_[nextSetAttr - 1]->getGeneralizationAndDelete(lhs, attr_num, nextSetAttr, specLhs);
        if (found){
            if (this->isFinalNode(attr_num)){
                this->rhsAttributes_.reset(attr_num);
            }

            specLhs.set(nextSetAttr);
        }
    }
    if (!found){
        found = this->getGeneralizationAndDelete(lhs, attr_num, nextSetAttr, specLhs);
    }
    return found;
}


bool FDTreeElement::getSpecialization
(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num, size_t currentAttr, std::bitset<kMaxAttrNum>& specLhsOut) const
{
    if (!this->rhsAttributes_[attr_num]){
        return false;
    }

    bool found = false;
    size_t attr = (currentAttr > 1 ? currentAttr : 1);
    size_t nextSetAttr = lhs._Find_next(currentAttr);

    if (nextSetAttr == kMaxAttrNum){
        while (!found && attr <= this->maxAttributeNumber_){
            if (this->children_[attr - 1] && this->children_[attr - 1]->getRhsAttributes()[attr_num]){
                found = this->children_[attr - 1]->getSpecialization(lhs, attr_num, currentAttr, specLhsOut);
            }
            ++attr;
        }
        if (found){
            specLhsOut.set(attr - 1);
        }
        return true;
    }

    while (!found && attr < nextSetAttr){
        if (this->children_[attr - 1] && this->children_[attr - 1]->getRhsAttributes()[attr_num]){
            found = this->children_[attr - 1]->getSpecialization(lhs, attr_num, currentAttr, specLhsOut);
        }
        ++attr;
    }
    if (!found && this->children_[nextSetAttr - 1] && this->children_[nextSetAttr - 1]->getRhsAttributes()[attr_num]){
        found = this->children_[nextSetAttr - 1]->getSpecialization(lhs, attr_num, nextSetAttr, specLhsOut);
    }

    specLhsOut.set(attr - 1, found);

    return found;
}

void FDTreeElement::addMostGeneralDependencies(){
    for (size_t i = 1; i <= this->maxAttributeNumber_; ++i){
        this->rhsAttributes_.set(i);
    }

    for (size_t i = 0; i < this->maxAttributeNumber_; ++i){
        this->isFd_[i] = true;
    }
}

void FDTreeElement::addFunctionalDependency(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num){
    FDTreeElement* currentNode = this;
    this->addRhsAttribute(attr_num);

    for (size_t i = lhs._Find_first(); i != kMaxAttrNum; i = lhs._Find_next(i)){
        if (currentNode->children_[i - 1] == nullptr){
            currentNode->children_[i - 1] = std::make_unique<FDTreeElement>(this->maxAttributeNumber_);
        }

        currentNode = currentNode->getChild(i - 1);
        currentNode->addRhsAttribute(attr_num);
    }

    currentNode->markAsLast(attr_num - 1);
}

void FDTreeElement::filterSpecializations(){
    std::bitset<kMaxAttrNum> activePath;
    auto filteredTree = std::make_unique<FDTreeElement>(this->maxAttributeNumber_);

    this->filterSpecializationsHelper(*filteredTree, activePath);

    this->children_ = std::move(filteredTree->children_);
    this->isFd_ = filteredTree->isFd_;
}

void FDTreeElement::filterSpecializationsHelper(FDTreeElement& filteredTree, std::bitset<kMaxAttrNum>& activePath){
    for (size_t attr = 1; attr <= this->maxAttributeNumber_; ++attr){
        if (this->children_[attr - 1]){
            activePath.set(attr);
            this->children_[attr - 1]->filterSpecializationsHelper(filteredTree, activePath);
            activePath.reset(attr);
        }
    }

    for (size_t attr = 1; attr <= this->maxAttributeNumber_; ++attr){
        std::bitset<kMaxAttrNum> specLhsOut; 
        if (this->isFd_[attr - 1] && !filteredTree.getSpecialization(activePath, attr, 0, specLhsOut)){
            filteredTree.addFunctionalDependency(activePath, attr);
        }
    }
}

 void FDTreeElement::printDep(const std::string& fileName, std::vector<std::string>& columnNames) const{
     std::ofstream file;
     file.open(fileName);
     std::bitset<kMaxAttrNum> activePath;
     printDependencies(activePath, file, columnNames);
     file.close();
 }

 void FDTreeElement::printDependencies(std::bitset<kMaxAttrNum>& activePath, std::ofstream& file,
 std::vector<std::string>& columnNames) const {
     std::string columnId;
     if (std::isdigit(columnNames[0][0])){
         columnId = "column";
     }
     std::string out;
     for (size_t attr = 1; attr <= this->maxAttributeNumber_; ++attr){
         if (this->isFd_[attr - 1]){
             out = "{";

             for (size_t i = activePath._Find_first(); i != kMaxAttrNum; i = activePath._Find_next(i)){
                 if (!columnId.empty())
                     out += columnId + std::to_string(std::stoi(columnNames[i - 1]) + 1) + ",";
                 else
                     out += columnNames[i - 1] + ",";
             }

             if (out.size() > 1){
                 out = out.substr(0, out.size() - 1);
             }
             if (!columnId.empty())
                 out += "} -> " + columnId + std::to_string(std::stoi(columnNames[attr - 1]) + 1);
             else
                 out += "} -> " + columnId + columnNames[attr - 1];
             file << out << std::endl;
         }
     }

     for (size_t attr = 1; attr <= this->maxAttributeNumber_; ++attr){
         if (this->children_[attr - 1]){
             activePath.set(attr);
             this->children_[attr - 1]->printDependencies(activePath, file, columnNames);
             activePath.reset(attr);
         }
     }

 }

 void FDTreeElement::fillFdCollection(const RelationalSchema &scheme, std::list<FD> &fdCollection) const {
    std::bitset<kMaxAttrNum> activePath;
    this->transformTreeFdCollection(activePath, fdCollection, scheme);
}

 void FDTreeElement::transformTreeFdCollection(std::bitset<kMaxAttrNum> &activePath, std::list<FD> &fdCollection,
                                               const RelationalSchema &scheme) const {
     for (size_t attr = 1; attr <= this->maxAttributeNumber_; ++attr){
         if (this->isFd_[attr - 1]){
             boost::dynamic_bitset<> lhs_bitset(this->maxAttributeNumber_);
             for (size_t i = activePath._Find_first(); i != kMaxAttrNum; i = activePath._Find_next(i)){
                 lhs_bitset.set(i - 1);
             }
             Vertical lhs(&scheme, lhs_bitset);
             Column rhs(&scheme, scheme.getColumn(attr - 1)->getName(), attr - 1);
             fdCollection.emplace_back(FD{lhs, rhs});
         }
     }

     for (size_t attr = 1; attr <= this->maxAttributeNumber_; ++attr){
         if (this->children_[attr - 1]){
             activePath.set(attr);
             this->children_[attr - 1]->transformTreeFdCollection(activePath, fdCollection, scheme);
             activePath.reset(attr);
         }
     }
}
