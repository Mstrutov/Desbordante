#include "ProfilingContext.h"

#include <utility>

#include <easylogging++.h>

#include "ListAgreeSetSample.h"
#include "PLICache.h"
#include "VerticalMap.h"

using std::shared_ptr;

ProfilingContext::ProfilingContext(Configuration configuration,
                                   ColumnLayoutRelationData* relation_data,
                                   std::function<void(const PartialKey&)> const& ucc_consumer,
                                   std::function<void(const PartialFD&)> const& fd_consumer,
                                   CachingMethod const& caching_method,
                                   CacheEvictionMethod const& eviction_method,
                                   double caching_method_value)
    : configuration_(std::move(configuration)),
      relation_data_(relation_data),
      random_(configuration_.seed == 0 ? std::mt19937() : std::mt19937(configuration_.seed)),
      custom_random_(configuration_.seed == 0 ? CustomRandom()
                                              : CustomRandom(configuration_.seed)) {
    ucc_consumer_ = ucc_consumer;
    fd_consumer_ = fd_consumer;
    // TODO: тут проявляется косяк, что unique_ptr<PLI> приходится отбирать у CLRD.
    //       SetSample и MaxEntropy требуют CLRD. Приходится плясать-переставлять методы
    //       да на самом деле в коде подразумевается, что в CLRD есть какая-то ссылка на PLI, так что
    //       надо переделывать
    if (configuration_.sample_size > 0) {
        auto schema = relation_data_->GetSchema();
        agree_set_samples_ =
            std::make_unique<util::BlockingVerticalMap<util::AgreeSetSample>>(schema);
        // TODO: сделать, чтобы при одном потоке agree_set_samples_ = std::make_unique<VerticalMap<AgreeSetSample>>(schema);
        for (auto& column : schema->GetColumns()) {
            CreateColumnFocusedSample(
                static_cast<Vertical>(*column),
                relation_data->GetColumnData(column->GetIndex()).GetPositionListIndex(),
                1);
        }
    } else {
        agree_set_samples_ = nullptr;
    }
    double max_entropy = GetMaximumEntropy(relation_data_);
    pli_cache_ = std::make_unique<util::PLICache>(
        relation_data_, caching_method, eviction_method, caching_method_value,
        GetMinEntropy(relation_data_), GetMeanEntropy(relation_data_),
        GetMedianEntropy(relation_data_), SetMaximumEntropy(relation_data_, caching_method),
        GetMedianGini(relation_data_), GetMedianInvertedEntropy(relation_data_));
    pli_cache_->SetMaximumEntropy(max_entropy);
    // TODO: partialFDScoring - for FD registration
}

ProfilingContext::~ProfilingContext() = default;

double ProfilingContext::GetMaximumEntropy(ColumnLayoutRelationData const* relation_data) {
    auto& columns = relation_data->GetColumnData();
    auto max_column = std::max_element(columns.begin(), columns.end(), [](auto& cd1, auto& cd2) {
        return cd1.GetPositionListIndex()->GetEntropy() < cd2.GetPositionListIndex()->GetEntropy();
    });
    return max_column->GetPositionListIndex()->GetEntropy();
}

double ProfilingContext::GetMinEntropy(ColumnLayoutRelationData const* relation_data) {
    auto& columns = relation_data->GetColumnData();
    auto min_column = std::min_element(columns.begin(), columns.end(), [](auto& cd1, auto& cd2) {
        return cd1.GetPositionListIndex()->GetEntropy() < cd2.GetPositionListIndex()->GetEntropy();
    });
    return min_column->GetPositionListIndex()->GetEntropy();
}

double ProfilingContext::GetMedianEntropy(ColumnLayoutRelationData const* relation_data) {
    std::vector<double> vals;

    for (auto& column : relation_data->GetColumnData()) {
        if (column.GetPositionListIndex()->GetEntropy() >= 0.001) {
            vals.push_back(column.GetPositionListIndex()->GetEntropy());
        }
    }

    return GetMedianValue(std::move(vals), "MedianEntropy");
}

double ProfilingContext::GetMedianInvertedEntropy(ColumnLayoutRelationData const* relation_data) {
    std::vector<double> vals;

    for (auto& column : relation_data->GetColumnData()) {
        if (column.GetPositionListIndex()->GetInvertedEntropy() >= 0.001) {
            vals.push_back(column.GetPositionListIndex()->GetInvertedEntropy());
        }
    }

    return GetMedianValue(std::move(vals), "MedianInvertedEntropy");
}

double ProfilingContext::GetMeanEntropy(ColumnLayoutRelationData const* relation_data) {
    double e = 0;

    for (auto& column : relation_data->GetColumnData()) {
        e += column.GetPositionListIndex()->GetEntropy();
    }
    return e / relation_data->GetColumnData().size();
}

double ProfilingContext::GetMedianGini(ColumnLayoutRelationData const* relation_data) {
    std::vector<double> vals;

    for (auto& column : relation_data->GetColumnData()) {
        if (column.GetPositionListIndex()->GetEntropy() >= 0.001) {    // getGini?
            vals.push_back(column.GetPositionListIndex()->GetGiniImpurity());
        }
    }

    return GetMedianValue(std::move(vals), "MedianGini");
}

double ProfilingContext::SetMaximumEntropy(ColumnLayoutRelationData const* relation_data,
                                           CachingMethod const& caching_method) {
    switch (caching_method) {
    case CachingMethod::kEntropy:
    case CachingMethod::kCoin:
    case CachingMethod::kNoCaching:
        return relation_data->GetMaximumEntropy();
    case CachingMethod::kTrueUniquenessEntropy:
        return GetMaximumEntropy(relation_data);
    case CachingMethod::kMeanEntropyThreshold:
        return GetMeanEntropy(relation_data);
    case CachingMethod::kHeuristicQ2:
        return GetMaximumEntropy(relation_data);
    case CachingMethod::kGini:
        return GetMedianGini(relation_data);
    case CachingMethod::kInvertedEntropy:
        return GetMedianInvertedEntropy(relation_data);
    default:
        return 0;
    }
}

util::AgreeSetSample const* ProfilingContext::CreateFocusedSample(Vertical const& focus,
                                                                  double boost_factor) {
    auto pli = pli_cache_->GetOrCreateFor(focus, this);
    auto pli_pointer = std::holds_alternative<util::PositionListIndex*>(pli)
                           ? std::get<util::PositionListIndex*>(pli)
                           : std::get<std::unique_ptr<util::PositionListIndex>>(pli).get();
    std::unique_ptr<util::ListAgreeSetSample> sample = util::ListAgreeSetSample::CreateFocusedFor(
        relation_data_,
        focus,
        pli_pointer,
        configuration_.sample_size * boost_factor,
        custom_random_
    );
    LOG(TRACE) << boost::format{"Creating sample focused on: %1%"} % focus.ToString();
    auto sample_ptr = sample.get();
    agree_set_samples_->Put(focus, std::move(sample));
    return sample_ptr;
}

util::AgreeSetSample const* ProfilingContext::CreateColumnFocusedSample(
    const Vertical& focus, util::PositionListIndex const* restriction_pli, double boost_factor) {
    std::unique_ptr<util::ListAgreeSetSample> sample = util::ListAgreeSetSample::CreateFocusedFor(
        relation_data_,
        focus,
        restriction_pli,
        configuration_.sample_size * boost_factor,
        custom_random_
    );
    LOG(TRACE) << boost::format{"Creating sample focused on: %1%"} % focus.ToString();
    auto sample_ptr = sample.get();
    agree_set_samples_->Put(focus, std::move(sample));
    return sample_ptr;
}

shared_ptr<util::AgreeSetSample const>
ProfilingContext::GetAgreeSetSample(Vertical const& focus) const {
    shared_ptr<util::AgreeSetSample const> sample = nullptr;
    for (auto& [key, next_sample] : agree_set_samples_->GetSubsetEntries(focus)) {
        if (sample == nullptr || next_sample->GetSamplingRatio() > sample->GetSamplingRatio()) {
            sample = next_sample;
        }
    }
    return sample;
}

double ProfilingContext::GetMedianValue(std::vector<double>&& values,
                                        std::string const& measure_name) {
    if (values.size() <= 1) {
        LOG(WARNING) << "Got " << measure_name << " == 0\n";
        return 0;
    }

    std::sort(values.begin(), values.end());
    return (values.size() % 2 == 0) ?
           ((values[values.size() / 2] + values[values.size() / 2 - 1]) / 2) :
           (values[values.size() / 2]);
}
