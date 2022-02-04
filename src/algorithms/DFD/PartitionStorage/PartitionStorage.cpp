#include <boost/optional.hpp>

#include "logging/easylogging++.h"

#include "PartitionStorage.h"
#include "VerticalMap.h"


util::PositionListIndex* PartitionStorage::Get(Vertical const &vertical) {
    return index_->Get(vertical).get();
}

PartitionStorage::PartitionStorage(ColumnLayoutRelationData* relation_data,
                                   CachingMethod caching_method,
                                   CacheEvictionMethod eviction_method) :
    relation_data_(relation_data),
    index_(std::make_unique<util::BlockingVerticalMap<util::PositionListIndex>>(relation_data->GetSchema())),
    caching_method_(caching_method),
    eviction_method_(eviction_method) {
    for (auto& column_ptr : relation_data->GetSchema()->GetColumns()) {
      index_->Put(
          static_cast<Vertical>(*column_ptr),
          relation_data->GetColumnData(column_ptr->GetIndex()).GetPliOwnership());
    }
}

PartitionStorage::~PartitionStorage() {}

// obtains or calculates a PositionListIndex using cache
std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>>
PartitionStorage::GetOrCreateFor(Vertical const &vertical) {
    std::scoped_lock lock(getting_pli_mutex_);
    LOG(DEBUG) << boost::format{"PLI for %1% requested: "} % vertical.ToString();

    // is PLI already cached?
    util::PositionListIndex* pli = Get(vertical);
    if (pli != nullptr) {
      pli->IncFreq();
        LOG(DEBUG) << boost::format{"Served from PLI cache."};
        //addToUsageCounter
        return pli;
    }
    // look for cached PLIs to construct the requested one
    auto subset_entries = index_->GetSubsetEntries(vertical);
    boost::optional<PositionListIndexRank> smallest_pli_rank;
    std::vector<PositionListIndexRank> ranks;
    ranks.reserve(subset_entries.size());
    for (auto& [sub_vertical, sub_pli_ptr] : subset_entries) {
        PositionListIndexRank pli_rank(&sub_vertical,
                                       std::const_pointer_cast<util::PositionListIndex>(sub_pli_ptr),
                                       sub_vertical.GetArity());
        ranks.push_back(pli_rank);
        if (!smallest_pli_rank
            || smallest_pli_rank->pli_->GetSize() > pli_rank.pli_->GetSize()
            || (smallest_pli_rank->pli_->GetSize() == pli_rank.pli_->GetSize() && smallest_pli_rank->added_arity_ < pli_rank.added_arity_)) {
            smallest_pli_rank = pli_rank;
        }
    }
    assert(smallest_pli_rank);            // check if smallest_pli_rank is initialized

    std::vector<PositionListIndexRank> operands;
    boost::dynamic_bitset<> cover(relation_data_->GetNumColumns());
    boost::dynamic_bitset<> cover_tester(relation_data_->GetNumColumns());
    if (smallest_pli_rank) {
      smallest_pli_rank->pli_->IncFreq();
        operands.push_back(*smallest_pli_rank);
        cover |= smallest_pli_rank->vertical_->GetColumnIndices();

        while (cover.count() < vertical.GetArity() && !ranks.empty()) {
            boost::optional<PositionListIndexRank> best_rank;
            //erase ranks with low added_arity_
            ranks.erase(std::remove_if(ranks.begin(), ranks.end(),
                                       [&cover_tester, &cover] (auto& rank) {
                                           cover_tester.reset();
                                           cover_tester |= rank.vertical_->GetColumnIndices();
                                           cover_tester -= cover;
                                           rank.added_arity_ = cover_tester.count();
                                           return rank.added_arity_ < 2;
                                       }), ranks.end());

            for (auto& rank : ranks) {
                if (!best_rank
                    || best_rank->added_arity_ < rank.added_arity_
                    || ( best_rank->added_arity_ == rank.added_arity_ && best_rank->pli_->GetSize() > rank.pli_->GetSize())) {
                    best_rank = rank;
                }
            }

            if (best_rank) {
              best_rank->pli_->IncFreq();
                operands.push_back(*best_rank);
                cover |= best_rank->vertical_->GetColumnIndices();
            }
        }
    }

    std::vector<std::unique_ptr<Vertical>> vertical_columns;

    for (auto& column : vertical.GetColumns()) {
        if (!cover[column->GetIndex()]) {
            vertical_columns.push_back(std::make_unique<Vertical>(static_cast<Vertical>(*column)));
            auto column_pli = index_->Get(**vertical_columns.rbegin());
            operands.emplace_back(vertical_columns.rbegin()->get(), column_pli, 1);
          column_pli->IncFreq();
        }
    }
    // sort operands by ascending order
    std::sort(operands.begin(), operands.end(),
              [](auto& el1, auto& el2) { return el1.pli_->GetSize() < el2.pli_->GetSize(); });

    LOG(DEBUG) << boost::format {"Intersecting %1%."} % "[UNIMPLEMENTED]";

    if (operands.empty()) {
        throw std::logic_error("Current implementation assumes operands.size() > 0");
    }

    // Intersect and cache
    std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>> variant_intersection_pli;
    if (operands.size() >= 4) {
        PositionListIndexRank base_pli_rank = operands[0];
        auto intersection_pli = base_pli_rank.pli_->ProbeAll(vertical.Without(*base_pli_rank.vertical_), *relation_data_);
        variant_intersection_pli = CachingProcess(vertical, std::move(intersection_pli));
    } else {
        Vertical current_vertical = *operands.begin()->vertical_;
        variant_intersection_pli = operands.begin()->pli_.get();

        for (size_t i = 1; i < operands.size(); i++) {
            current_vertical = current_vertical.Union(*operands[i].vertical_);
            variant_intersection_pli =
                    std::holds_alternative<util::PositionListIndex*>(variant_intersection_pli)
                    ? std::get<util::PositionListIndex *>(variant_intersection_pli)->Intersect(operands[i].pli_.get())
                    : std::get<std::unique_ptr<util::PositionListIndex>>(variant_intersection_pli)->Intersect(operands[i].pli_.get());
            variant_intersection_pli = CachingProcess(
                current_vertical,
                std::move(std::get<std::unique_ptr<util::PositionListIndex>>(variant_intersection_pli)));
        }
    }

    LOG(DEBUG) << boost::format {"Calculated from %1% sub-PLIs (saved %2% intersections)."}
                  % operands.size() % (vertical.GetArity() - operands.size());

    return variant_intersection_pli;
}

size_t PartitionStorage::Size() const {
    return index_->GetSize();
}

std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>> PartitionStorage::CachingProcess(
        Vertical const &vertical, std::unique_ptr<util::PositionListIndex> pli) {
    auto pli_pointer = pli.get();
  index_->Put(vertical, std::move(pli));
    return pli_pointer;
}
