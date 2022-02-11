#include "PLICache.h"

#include <boost/optional.hpp>
#include <easylogging++.h>

#include "VerticalMap.h"

namespace util {

PositionListIndex* PLICache::Get(Vertical const& vertical) {
    return index_->Get(vertical).get();
}

PLICache::PLICache(ColumnLayoutRelationData* relation_data, CachingMethod caching_method,
                   CacheEvictionMethod eviction_method, double caching_method_value,
                   double min_entropy, double mean_entropy, double median_entropy,
                   double maximum_entropy, double median_gini, double median_inverted_entropy)
    : relation_data_(relation_data),
      // TODO: сделать
      // index_(std::make_unique<VerticalMap<PositionListIndex>>(relation_data->GetSchema())) при
      // одном потоке
      index_(std::make_unique<BlockingVerticalMap<PositionListIndex>>(relation_data->GetSchema())),
      caching_method_(caching_method),
      eviction_method_(eviction_method),
      caching_method_value_(caching_method_value),
      maximum_entropy_(maximum_entropy),
      mean_entropy_(mean_entropy),
      min_entropy_(min_entropy),
      median_entropy_(median_entropy),
      median_gini_(median_gini),
      median_inverted_entropy_(median_inverted_entropy) {
    for (auto& column_ptr : relation_data->GetSchema()->GetColumns()) {
        index_->Put(static_cast<Vertical>(*column_ptr),
                    relation_data->GetColumnData(column_ptr->GetIndex()).GetPliOwnership());
    }
}

PLICache::~PLICache() {
    for (auto& column_ptr : relation_data_->GetSchema()->GetColumns()) {
        //auto PLI =
        index_->Remove(static_cast<Vertical>(*column_ptr));
        //relation_data_->GetColumnData(column_ptr->getIndex()).getPLI(std::move(PLI));
    }
}

// obtains or calculates a PositionListIndex using cache
std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> PLICache::GetOrCreateFor(
    Vertical const& vertical, ProfilingContext* profiling_context) {
    std::scoped_lock lock(getting_pli_mutex_);
    LOG(DEBUG) << boost::format{"PLI for %1% requested: "} % vertical.ToString();

    // is PLI already cached?
    PositionListIndex* pli = Get(vertical);
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
        // TODO: избавиться от таких const_cast, которые сбрасывают константность
        PositionListIndexRank pli_rank(&sub_vertical,
                                       std::const_pointer_cast<PositionListIndex>(sub_pli_ptr),
                                       sub_vertical.GetArity());
        ranks.push_back(pli_rank);
        if (!smallest_pli_rank || smallest_pli_rank->pli_->GetSize() > pli_rank.pli_->GetSize() ||
            (smallest_pli_rank->pli_->GetSize() == pli_rank.pli_->GetSize() &&
             smallest_pli_rank->added_arity_ < pli_rank.added_arity_)) {
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
                                       [&cover_tester, &cover](auto& rank) {
                                           cover_tester.reset();
                                           cover_tester |= rank.vertical_->GetColumnIndices();
                                           cover_tester -= cover;
                                           rank.added_arity_ = cover_tester.count();
                                           return rank.added_arity_ < 2;
                                       }),
                        ranks.end());

            for (auto& rank : ranks) {
                if (!best_rank || best_rank->added_arity_ < rank.added_arity_ ||
                    (best_rank->added_arity_ == rank.added_arity_ &&
                     best_rank->pli_->GetSize() > rank.pli_->GetSize())) {
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

    // TODO: конкретные костыли, надо делать Column : Vertical
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
    // TODO: Profiling context stuff

    LOG(DEBUG) << boost::format{"Intersecting %1%."} % "[UNIMPLEMENTED]";

    if (operands.empty()) {
        throw std::logic_error("Current implementation assumes operands.size() > 0");
    }

    // TODO: тут не очень понятно: CachingProcess может забрать себе PLI, а может и отдать обратно,
    //  поэтому приходится через variant разбирать. Проверить, насколько много платим за обёртку.
    // Intersect and cache
    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> variant_intersection_pli;
    if (operands.size() >= profiling_context->GetConfiguration().nary_intersection_size) {
        PositionListIndexRank base_pli_rank = operands[0];
        auto intersection_pli = base_pli_rank.pli_->ProbeAll(
            vertical.Without(*base_pli_rank.vertical_), *relation_data_);
        variant_intersection_pli =
            CachingProcess(vertical, std::move(intersection_pli), profiling_context);
    } else {
        Vertical current_vertical = *operands.begin()->vertical_;
        variant_intersection_pli = operands.begin()->pli_.get();

        for (size_t i = 1; i < operands.size(); i++) {
            current_vertical = current_vertical.Union(*operands[i].vertical_);
            variant_intersection_pli =
                std::holds_alternative<PositionListIndex*>(variant_intersection_pli)
                    ? std::get<PositionListIndex*>(variant_intersection_pli)
                          ->Intersect(operands[i].pli_.get())
                    : std::get<std::unique_ptr<PositionListIndex>>(variant_intersection_pli)
                          ->Intersect(operands[i].pli_.get());
            variant_intersection_pli = CachingProcess(
                current_vertical,
                std::move(std::get<std::unique_ptr<PositionListIndex>>(variant_intersection_pli)),
                profiling_context);
        }
    }

    LOG(DEBUG) << boost::format{"Calculated from %1% sub-PLIs (saved %2% intersections)."} %
                      operands.size() % (vertical.GetArity() - operands.size());

    return variant_intersection_pli;
}

size_t PLICache::Size() const {
    return index_->GetSize();
}

std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> PLICache::CachingProcess(
    Vertical const& vertical, std::unique_ptr<PositionListIndex> pli,
    ProfilingContext* profiling_context) {
    auto pli_pointer = pli.get();
    switch (caching_method_) {
    case CachingMethod::kCoin:
        if (profiling_context->NextDouble() <
            profiling_context->GetConfiguration().caching_probability) {
            index_->Put(vertical, std::move(pli));
            return pli_pointer;
        } else {
            return pli;
        }
    case CachingMethod::kNoCaching:
        return pli;
    case CachingMethod::kAllCaching:
        index_->Put(vertical, std::move(pli));
        return pli_pointer;
    default:
        throw std::runtime_error(
            "Only kNoCaching and kAllCaching strategies are currently available");
    }
}

} // namespace util

