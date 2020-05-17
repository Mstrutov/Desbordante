#include "core/ProfilingContext.h"
#include "util/ListAgreeSetSample.h"
#include "util/PLICache.h"
#include "util/VerticalMap.h"

double ProfilingContext::getMaximumEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData) {
    auto columns = relationData->getColumnData();
    auto maxColumn = std::max_element(columns.begin(), columns.end(),
            [](auto ptr1, auto ptr2) {
        return ptr1->getPositionListIndex()->getEntropy() < ptr2->getPositionListIndex()->getEntropy();
    });
    return (*maxColumn)->getPositionListIndex()->getEntropy();
}

double ProfilingContext::getMinEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData) {
    auto columns = relationData->getColumnData();
    auto minColumn = std::min_element(columns.begin(), columns.end(),
                                      [](auto ptr1, auto ptr2) {
                                          return ptr1->getPositionListIndex()->getEntropy() < ptr2->getPositionListIndex()->getEntropy();
                                      });
    return (*minColumn)->getPositionListIndex()->getEntropy();
}

double ProfilingContext::getMedianEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData) {
    std::vector<double> vals;

    for (auto& column : relationData->getColumnData()) {
        if (column->getPositionListIndex()->getEntropy() >= 0.001) {
            vals.push_back(column->getPositionListIndex()->getEntropy());
        }
    }

    std::sort(vals.begin(), vals.end());
    return (vals.size() % 2 == 0) ?
           ((vals[vals.size() / 2] + vals[vals.size() / 2 - 1]) / 2) :
           (vals[vals.size() / 2]);
}

double ProfilingContext::getMedianInvertedEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData) {
    std::vector<double> vals;

    for (auto& column : relationData->getColumnData()) {
        if (column->getPositionListIndex()->getInvertedEntropy() >= 0.001) {
            vals.push_back(column->getPositionListIndex()->getInvertedEntropy());
        }
    }

    std::sort(vals.begin(), vals.end());
    return (vals.size() % 2 == 0) ?
           ((vals[vals.size() / 2] + vals[vals.size() / 2 - 1]) / 2) :
           (vals[vals.size() / 2]);
}

double ProfilingContext::getMeanEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData) {
    double e = 0;

    for (auto& column : relationData->getColumnData()) {
        e += column->getPositionListIndex()->getEntropy();
    }
    return e / relationData->getColumnData().size();
}

double ProfilingContext::getMedianGini(std::shared_ptr<ColumnLayoutRelationData> relationData) {
    std::vector<double> vals;

    for (auto& column : relationData->getColumnData()) {
        if (column->getPositionListIndex()->getEntropy() >= 0.001) {    // getGini?
            vals.push_back(column->getPositionListIndex()->getGiniImpurity());
        }
    }

    std::sort(vals.begin(), vals.end());
    return (vals.size() % 2 == 0) ?
           ((vals[vals.size() / 2] + vals[vals.size() / 2 - 1]) / 2) :
           (vals[vals.size() / 2]);
}

double ProfilingContext::setMaximumEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData,
                                           CachingMethod const &cachingMethod) {
    switch (cachingMethod) {
        case CachingMethod::ENTROPY:
        case CachingMethod::COIN:
        case CachingMethod::NOCACHING:
            return relationData->getMaximumEntropy();
        case CachingMethod::TRUEUNIQUENESSENTROPY:
            return getMaximumEntropy(relationData);
        case CachingMethod::MEANENTROPYTHRESHOLD:
            return getMeanEntropy(relationData);
        case CachingMethod::HEURISTICQ2:
            return getMaximumEntropy(relationData);
        case CachingMethod::GINI:
            return getMedianGini(relationData);
        case CachingMethod::INVERTEDENTROPY:
            return getMedianInvertedEntropy(relationData);
        default:
            return 0;
    }
}

std::shared_ptr<AgreeSetSample>
ProfilingContext::createFocusedSample(std::shared_ptr<Vertical> focus, double boostFactor) {
    std::shared_ptr<ListAgreeSetSample> sample = ListAgreeSetSample::createFocusedFor(
            relationData_,
            focus,
            pliCache_->getOrCreateFor(*focus, *this),
            configuration_.sampleSize * boostFactor
            );
    agreeSetSamples_->put(*focus, sample);
    return sample;
}

// TODO: implement using std::max_element??
std::shared_ptr<AgreeSetSample> ProfilingContext::getAgreeSetSample(std::shared_ptr<Vertical> focus) {
    std::shared_ptr<AgreeSetSample> sample = nullptr;
    for (auto& [key, nextSample] : agreeSetSamples_->getSubsetEntries(*focus)) {
        if (sample == nullptr || nextSample->getSamplingRatio() > sample->getSamplingRatio()) {
            sample = nextSample;
        }
    }
    return std::shared_ptr<AgreeSetSample>();
}

