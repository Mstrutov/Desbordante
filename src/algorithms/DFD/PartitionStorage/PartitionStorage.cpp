//
// Created by alex on 22.04.2021.
//

#include <boost/optional.hpp>

#include "logging/easylogging++.h"

#include "PartitionStorage.h"
#include "VerticalMap.h"


std::shared_ptr<PositionListIndex> PartitionStorage::get(Vertical const &vertical) {
    //VerticalMap<std::shared_ptr<PositionListIndex>> obj(nullptr);
    return index_->get(vertical);
}

PartitionStorage::PartitionStorage(std::shared_ptr<ColumnLayoutRelationData> relationData, CachingMethod cachingMethod,
                   CacheEvictionMethod evictionMethod) :
        relationData_(relationData),
        index_(std::make_shared<CacheMap>(relationData->getSchema())),
        cachingMethod_(cachingMethod),
        evictionMethod_(evictionMethod) {
    for (auto& column_ptr : relationData->getSchema()->getColumns()) {
        index_->put(static_cast<Vertical>(*column_ptr), relationData->getColumnData(column_ptr->getIndex())->getPositionListIndex());
    }
    // usageCounter, decayCounter, statistics
}

// obtains or calculates a PositionListIndex using cache
std::shared_ptr<PositionListIndex>
PartitionStorage::getOrCreateFor(Vertical const &vertical) {
    LOG(DEBUG) << boost::format{"PLI for %1% requested: "} % vertical.toString();

    // is PLI already cached?
    std::shared_ptr<PositionListIndex> pli = get(vertical);
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
        PositionListIndexRank pliRank(subVertical, subPLI_ptr, subVertical->getArity());
        ranks.push_back(pliRank);
        if (!smallestPliRank
            || smallestPliRank->pli_->getSize() > pliRank.pli_->getSize()
            || (smallestPliRank->pli_->getSize() == pliRank.pli_->getSize() && smallestPliRank->addedArity_ < pliRank.addedArity_)) {
            smallestPliRank = pliRank;
        }
    }
    assert(smallestPliRank);            // check if smallestPliRank is initialized

    std::vector<PositionListIndexRank> operands;
    boost::dynamic_bitset<> cover(relationData_.lock()->getNumColumns());
    boost::dynamic_bitset<> coverTester(relationData_.lock()->getNumColumns());
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
    for (auto& column : vertical.getColumns()) {
        if (!cover[column->getIndex()]) {
            auto columnData = relationData_.lock()->getColumnData(column->getIndex());
            operands.emplace_back(std::make_shared<Vertical>(static_cast<Vertical>(*column)),
                                  columnData->getPositionListIndex(), 1);
            columnData->getPositionListIndex()->incFreq();
        }
    }
    // sort operands by ascending order
    std::sort(operands.begin(), operands.end(), [](auto& el1, auto& el2) { return el1.pli_->getSize() < el2.pli_->getSize(); });
    // TODO: Profiling context stuff

    LOG(DEBUG) << boost::format {"Intersecting %1%."} % "[UNIMPLEMENTED]";
    // Intersect and cache
    if (operands.size() >= 4) {
        PositionListIndexRank basePliRank = operands[0];
        pli = basePliRank.pli_->probeAll(*vertical.without(*basePliRank.vertical_), *relationData_.lock());
        //cachingProcess(vertical, pli, profilingContext);
        index_->put(vertical, pli);
    } else {
        std::shared_ptr<Vertical> currentVertical = nullptr;
        for (auto& operand : operands) {
            if (pli == nullptr) {
                currentVertical = operand.vertical_;
                pli = operand.pli_;
            } else {
                currentVertical = currentVertical->Union(*operand.vertical_);
                pli = pli->intersect(operand.pli_);
                //cachingProcess(*currentVertical, pli, profilingContext);
                index_->put(*currentVertical, pli);
            }
        }
    }

    LOG(DEBUG) << boost::format {"Calculated from %1% sub-PLIs (saved %2% intersections)."}
                  % operands.size() % (vertical.getArity() - operands.size());

    return pli;
}

size_t PartitionStorage::size() const {
    return index_->getSize();
}

/*void PartitionStorage::cachingProcess(Vertical const &vertical, std::shared_ptr<PositionListIndex> pli, ProfilingContext* profilingContext) {
    switch (cachingMethod_) {
        case CachingMethod::COIN:
            if (profilingContext->customRandom_.nextDouble() < profilingContext->configuration_.cachingProbability) {
                index_->put(vertical, pli);
            }
            break;
        case CachingMethod::NOCACHING:
            //index_->put(vertical, pli);
            // newUsageInfo - parallel
            break;
        case CachingMethod::ALLCACHING:
            index_->put(vertical, pli);
            break;
        default:
            //index_->put(vertical, pli);
            // doubts on necessity of statistics => no implementation yet
            throw std::runtime_error("Only NOCACHING and ALLCACHING strategies are currently available");
            //break;
    }
}*/

