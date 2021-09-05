#include <boost/optional.hpp>

#include "logging/easylogging++.h"

#include "PartitionStorage.h"
#include "VerticalMap.h"


PositionListIndex* PartitionStorage::get(Vertical const &vertical) {
    return index_->get(vertical).get();
}

PartitionStorage::PartitionStorage(ColumnLayoutRelationData* relationData, CachingMethod cachingMethod,
                   CacheEvictionMethod evictionMethod) :
        relationData_(relationData),
        index_(std::make_unique<BlockingVerticalMap<PositionListIndex>>(relationData->getSchema())),
        cachingMethod_(cachingMethod),
        evictionMethod_(evictionMethod) {
    for (auto& column_ptr : relationData->getSchema()->getColumns()) {
        index_->put(
                static_cast<Vertical>(*column_ptr),
                relationData->getColumnData(column_ptr->getIndex()).getPLIOwnership());
    }
}

PartitionStorage::~PartitionStorage() {}

// obtains or calculates a PositionListIndex using cache
std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> PartitionStorage::getOrCreateFor(
        Vertical const &vertical) {
    std::scoped_lock lock(gettingPLIMutex);
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
        PositionListIndexRank pliRank(&subVertical, std::const_pointer_cast<PositionListIndex>(subPLI_ptr), subVertical.getArity());
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

    LOG(DEBUG) << boost::format {"Intersecting %1%."} % "[UNIMPLEMENTED]";

    if (operands.empty()) {
        throw std::logic_error("Current implementation assumes operands.size() > 0");
    }

    // Intersect and cache
    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> variantIntersectionPLI;
    if (operands.size() >= 4) {
        PositionListIndexRank basePliRank = operands[0];
        auto intersectionPLI = basePliRank.pli_->probeAll(vertical.without(*basePliRank.vertical_), *relationData_);
        variantIntersectionPLI = cachingProcess(vertical, std::move(intersectionPLI));
    } else {
        Vertical currentVertical = *operands.begin()->vertical_;
        variantIntersectionPLI = operands.begin()->pli_.get();

        for (size_t i = 1; i < operands.size(); i++) {
            currentVertical = currentVertical.Union(*operands[i].vertical_);
            variantIntersectionPLI =
                    std::holds_alternative<PositionListIndex*>(variantIntersectionPLI)
                    ? std::get<PositionListIndex*>(variantIntersectionPLI)->intersect(operands[i].pli_.get())
                    : std::get<std::unique_ptr<PositionListIndex>>(variantIntersectionPLI)->intersect(operands[i].pli_.get());
            variantIntersectionPLI = cachingProcess(
                    currentVertical,
                    std::move(std::get<std::unique_ptr<PositionListIndex>>(variantIntersectionPLI)));
        }
    }

    LOG(DEBUG) << boost::format {"Calculated from %1% sub-PLIs (saved %2% intersections)."}
                  % operands.size() % (vertical.getArity() - operands.size());

    return variantIntersectionPLI;
}

size_t PartitionStorage::size() const {
    return index_->getSize();
}

std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> PartitionStorage::cachingProcess(
        Vertical const &vertical, std::unique_ptr<PositionListIndex> pli) {
    auto pliPointer = pli.get();
    index_->put(vertical, std::move(pli));
    return pliPointer;
}
