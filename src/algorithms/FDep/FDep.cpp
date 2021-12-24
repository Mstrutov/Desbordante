#include "FDep.h"
#include "ColumnLayoutRelationData.h"

#include <chrono>

//#ifndef PRINT_FDS
//#define PRINT_FDS
//#endif

FDep::FDep(const std::filesystem::path &path, char separator, bool hasHeader):
        FDAlgorithm(path, separator, hasHeader){}

unsigned long long FDep::executeInternal(){

    initialize();

    auto startTime = std::chrono::system_clock::now();

    buildNegativeCover();

    this->tuples_.shrink_to_fit();

    this->posCoverTree_ = std::make_unique<FDTreeElement>(this->numberAttributes_);
    this->posCoverTree_->addMostGeneralDependencies();

    std::bitset<FDTreeElement::kMaxAttrNum> activePath;
    calculatePositiveCover(*this->negCoverTree_, activePath);

    posCoverTree_->fillFdCollection(*this->schema_, fdCollection_);

    auto elapsed_milliseconds =
    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - startTime);

#ifdef PRINT_FDS
    posCoverTree_->printDep("recent_call_result.txt", this->columnNames_);
#endif

    return elapsed_milliseconds.count(); 
}

void FDep::initialize(){
    loadData();
}

void FDep::buildNegativeCover() {
    this->negCoverTree_ = std::make_unique<FDTreeElement>(this->numberAttributes_);
    for (auto i = this->tuples_.begin(); i != this->tuples_.end(); ++i){
        for (auto j = i + 1; j != this->tuples_.end(); ++j)
            addViolatedFDs(*i, *j);
    }

    this->negCoverTree_->filterSpecializations();
}

void FDep::addViolatedFDs(std::vector<size_t> const & t1, std::vector<size_t> const & t2){
    std::bitset<FDTreeElement::kMaxAttrNum> equalAttr((2 << this->numberAttributes_) - 1);
    equalAttr.reset(0);
    std::bitset<FDTreeElement::kMaxAttrNum> diffAttr;

    for (size_t attr = 0; attr < this->numberAttributes_; ++attr){
        diffAttr[attr + 1] = (t1[attr] != t2[attr]);
    }

    equalAttr &= (~diffAttr);
    for (size_t attr = diffAttr._Find_first(); attr != FDTreeElement::kMaxAttrNum; attr = diffAttr._Find_next(attr)){
        this->negCoverTree_->addFunctionalDependency(equalAttr, attr);
    }
}

void FDep::calculatePositiveCover(FDTreeElement const & negCoverSubtree, std::bitset<FDTreeElement::kMaxAttrNum> & activePath){
    for (size_t attr = 1; attr <= this->numberAttributes_; ++attr){
        if (negCoverSubtree.checkFd(attr - 1)){
            this->specializePositiveCover(activePath, attr);
        }
    }

    for (size_t attr = 1; attr <= this->numberAttributes_; ++attr){
        if (negCoverSubtree.getChild(attr - 1)){
            activePath.set(attr);
            this->calculatePositiveCover(*negCoverSubtree.getChild(attr - 1), activePath);
            activePath.reset(attr);
        }
    }

}

void FDep::specializePositiveCover(std::bitset<FDTreeElement::kMaxAttrNum> const & lhs, size_t const & a){
    std::bitset<FDTreeElement::kMaxAttrNum> specLhs;

    while (this->posCoverTree_->getGeneralizationAndDelete(lhs, a, 0, specLhs))
    {
        for (size_t attr = this->numberAttributes_; attr > 0; --attr){
            if (!lhs.test(attr) && (attr != a)){
                specLhs.set(attr);
                if (!this->posCoverTree_->containsGeneralization(specLhs, a, 0)){
                    this->posCoverTree_->addFunctionalDependency(specLhs, a);
                }
                specLhs.reset(attr);
            }
        }

        specLhs.reset();
    }
}


void FDep::loadData(){
    this->numberAttributes_ = inputGenerator_.getNumberOfColumns();
    if (this->numberAttributes_ == 0){
        throw std::runtime_error("Unable to work on empty dataset. Check data file");
    }
    this->columnNames_.resize(this->numberAttributes_);

    this->schema_ = std::make_unique<RelationalSchema>(inputGenerator_.getRelationName(), true);

    for (size_t i = 0; i < this->numberAttributes_; ++i){
        this->columnNames_[i] = inputGenerator_.getColumnName(static_cast<int>(i));
        this->schema_->appendColumn(this->columnNames_[i]);
    }

    std::vector<std::string> nextLine; 
    while (inputGenerator_.getHasNext()){
        nextLine = inputGenerator_.parseNext();
        if (nextLine.empty()) break;
        this->tuples_.emplace_back(std::vector<size_t>(this->numberAttributes_));
        for (size_t i = 0; i < this->numberAttributes_; ++i){
            this->tuples_.back()[i] = std::hash<std::string>{}(nextLine[i]);
        }
    } 

}
