#include <boost/optional.hpp>

#include "logging/easylogging++.h"

#include "PLICache.h"
#include "VerticalMap.h"


PositionListIndex* PLICache::get(Vertical const &vertical) {
    return index_->get(vertical);
}

PLICache::PLICache(ColumnLayoutRelationData* relationData, CachingMethod cachingMethod,
                   CacheEvictionMethod evictionMethod, double cachingMethodValue, double minEntropy, double meanEntropy,
                   double medianEntropy, double maximumEntropy, double medianGini, double medianInvertedEntropy) :
        relationData_(relationData),
        index_(std::make_unique<CacheMap>(relationData->getSchema())),
        cachingMethod_(cachingMethod),
        evictionMethod_(evictionMethod),
        cachingMethodValue_(cachingMethodValue),
        maximumEntropy_(maximumEntropy),
        meanEntropy_(meanEntropy),
        minEntropy_(minEntropy),
        medianEntropy_(medianEntropy),
        medianGini_(medianGini),
        medianInvertedEntropy_(medianInvertedEntropy) {
    for (auto& column_ptr : relationData->getSchema()->getColumns()) {
        index_->put(
                static_cast<Vertical>(*column_ptr),
                relationData->getColumnData(column_ptr->getIndex()).moveOutPositionListIndex());
    }
}

PLICache::~PLICache() {
    for (auto& column_ptr : relationData_->getSchema()->getColumns()) {
        auto PLI = index_->remove(static_cast<Vertical>(*column_ptr));
        relationData_->getColumnData(column_ptr->getIndex()).moveInPositionListIndex(std::move(PLI));
    }
}

// obtains or calculates a PositionListIndex using cache
PositionListIndex* PLICache::getOrCreateFor(Vertical const &vertical, ProfilingContext* profilingContext) {
    LOG(DEBUG) << boost::format{"PLI for %1% requested: "} % vertical.toString();

    // is PLI already cached?
    PositionListIndex* pli = get(vertical);
    if (pli != nullptr) {
        pli->incFreq();
        LOG(DEBUG) << boost::format{"Served from PLI cache."};
        //addToUsageCounter
        return pli;
    }
    // look for cached PLIs to construct the requested one
    auto subsetEntries = index_->getSubsetEntries(vertical);
    boost::optional<PositionListIndexRank> smallestPliRank;
    std::vector<PositionListIndexRank> ranks;
    ranks.reserve(subsetEntries.size());
    for (auto& [subVertical, subPLI_ptr] : subsetEntries) {
        // TODO: избавиться от таких const_cast, которые сбрасывают константность
        PositionListIndexRank pliRank(&subVertical, const_cast<PositionListIndex*>(subPLI_ptr), subVertical.getArity());
        ranks.push_back(pliRank);
        if (!smallestPliRank
            || smallestPliRank->pli_->getSize() > pliRank.pli_->getSize()
            || (smallestPliRank->pli_->getSize() == pliRank.pli_->getSize() && smallestPliRank->addedArity_ < pliRank.addedArity_)) {
            smallestPliRank = pliRank;
        }
    }
    assert(smallestPliRank);            // check if smallestPliRank is initialized

    std::vector<PositionListIndexRank> operands;
    boost::dynamic_bitset<> cover(relationData_->getNumColumns());
    boost::dynamic_bitset<> coverTester(relationData_->getNumColumns());
    if (smallestPliRank) {
        smallestPliRank->pli_->incFreq();
        operands.push_back(*smallestPliRank);
        cover |= smallestPliRank->vertical_->getColumnIndices();

        while (cover.count() < vertical.getArity() && !ranks.empty()) {
            boost::optional<PositionListIndexRank> bestRank;
            //erase ranks with low addedArity_
            ranks.erase(std::remove_if(ranks.begin(), ranks.end(),
                    [&coverTester, &cover] (auto& rank) {
                        coverTester.reset();
                        coverTester |= rank.vertical_->getColumnIndices();
                        coverTester -= cover;
                        rank.addedArity_ = coverTester.count();
                        return rank.addedArity_ < 2;
                    }), ranks.end());

            for (auto& rank : ranks) {
                if (!bestRank
                    || bestRank->addedArity_ < rank.addedArity_
                    || ( bestRank->addedArity_ == rank.addedArity_ && bestRank->pli_->getSize() > rank.pli_->getSize())) {
                    bestRank = rank;
                }
            }

            if (bestRank) {
                bestRank->pli_->incFreq();
                operands.push_back(*bestRank);
                cover |= bestRank->vertical_->getColumnIndices();
            }
        }
    }

    // TODO: конкретные костыли, надо делать Column : Vertical
    std::vector<std::unique_ptr<Vertical>> verticalColumns;

    for (auto& column : vertical.getColumns()) {
        if (!cover[column->getIndex()]) {
            verticalColumns.push_back(std::make_unique<Vertical>(static_cast<Vertical>(*column)));
            auto columnPLI = index_->get(**verticalColumns.rbegin());
            operands.emplace_back(verticalColumns.rbegin()->get(), columnPLI, 1);
            columnPLI->incFreq();
        }
    }
    // sort operands by ascending order
    std::sort(operands.begin(), operands.end(),
              [](auto& el1, auto& el2) { return el1.pli_->getSize() < el2.pli_->getSize(); });
    // TODO: Profiling context stuff

    LOG(DEBUG) << boost::format {"Intersecting %1%."} % "[UNIMPLEMENTED]";
    // Intersect and cache
    std::unique_ptr<PositionListIndex> intersectionPLI;
    if (operands.size() >= profilingContext->getConfiiguration().naryIntersectionSize) {
        PositionListIndexRank basePliRank = operands[0];
        intersectionPLI = basePliRank.pli_->probeAll(vertical.without(*basePliRank.vertical_), *relationData_);
        cachingProcess(vertical, std::move(intersectionPLI), profilingContext);
    } else {
        Vertical const* currentVertical = nullptr;
        for (auto& operand : operands) {
            if (pli == nullptr) {
                currentVertical = operand.vertical_;
                pli = operand.pli_;
            } else {
                cachingProcess(
                        currentVertical->Union(*operand.vertical_),
                        pli->intersect(operand.pli_),
                        profilingContext);
            }
        }
    }

    LOG(DEBUG) << boost::format {"Calculated from %1% sub-PLIs (saved %2% intersections)."}
        % operands.size() % (vertical.getArity() - operands.size());

    return pli;
}

size_t PLICache::size() const {
    return index_->getSize();
}

void PLICache::cachingProcess(Vertical const &vertical, std::unique_ptr<PositionListIndex> pli, ProfilingContext* profilingContext) {
    switch (cachingMethod_) {
        case CachingMethod::COIN:
            if (profilingContext->nextDouble() < profilingContext->getConfiiguration().cachingProbability) {
                index_->put(vertical, std::move(pli));
            }
            break;
        case CachingMethod::NOCACHING:
            //index_->put(vertical, pli);
            // newUsageInfo - parallel
            break;
        case CachingMethod::ALLCACHING:
            index_->put(vertical, std::move(pli));
            break;
        default:
            //index_->put(vertical, pli);
            // doubts on necessity of statistics => no implementation yet
            throw std::runtime_error("Only NOCACHING and ALLCACHING strategies are currently available");
            //break;
    }
}

