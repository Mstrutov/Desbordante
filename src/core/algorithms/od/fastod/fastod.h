#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/algorithm.h"
#include "algorithms/od/fastod/model/attribute_pair.h"
#include "algorithms/od/fastod/model/attribute_set.h"
#include "algorithms/od/fastod/model/canonical_od.h"
#include "algorithms/od/fastod/storage/partition_cache.h"
#include "algorithms/od/fastod/util/timer.h"
#include "config/tabular_data/input_table_type.h"
#include "config/time_limit/type.h"

namespace algos {

class Fastod : public Algorithm {
private:
    using AscCanonicalOD = fastod::AscCanonicalOD;
    using DescCanonicalOD = fastod::DescCanonicalOD;
    using SimpleCanonicalOD = fastod::SimpleCanonicalOD;
    using AttributePair = fastod::AttributePair;
    using AttributeSet = fastod::AttributeSet;
    using PartitionCache = fastod::PartitionCache;
    using DataFrame = fastod::DataFrame;
    using Timer = fastod::Timer;

    config::TimeLimitSecondsType time_limit_seconds_ = 0u;
    bool is_complete_ = true;
    size_t level_ = 1;
    size_t od_count_ = 0;
    size_t fd_count_ = 0;
    size_t ocd_count_ = 0;

    std::vector<AscCanonicalOD> result_asc_;
    std::vector<DescCanonicalOD> result_desc_;
    std::vector<SimpleCanonicalOD> result_simple_;
    std::unordered_set<AttributeSet> context_in_current_level_;
    std::unordered_map<AttributeSet, AttributeSet> cc_;
    std::unordered_map<AttributeSet, std::unordered_set<AttributePair>> cs_asc_;
    std::unordered_map<AttributeSet, std::unordered_set<AttributePair>> cs_desc_;

    PartitionCache partition_cache_;

    AttributeSet schema_;
    std::shared_ptr<DataFrame> data_;
    config::InputTable input_table_;

    Timer timer_;

    void MakeExecuteOptsAvailable() override;
    void LoadDataInternal() override;
    void ResetState() override;
    unsigned long long ExecuteInternal() final;

    void PrepareOptions();
    void RegisterOptions();
    void MakeLoadOptionsAvailable();

    template <bool Ascending>
    void AddToResult(fastod::CanonicalOD<Ascending>&& od) {
        if constexpr (Ascending)
            result_asc_.emplace_back(std::move(od));
        else
            result_desc_.emplace_back(std::move(od));
    }

    void AddToResult(SimpleCanonicalOD&& od) {
        result_simple_.emplace_back(std::move(od));
    }

    bool IsTimeUp() const;

    void CCPut(AttributeSet const& key, AttributeSet attribute_set);
    AttributeSet const& CCGet(AttributeSet const& key);

    template <bool Ascending>
    void CSPut(AttributeSet const& key, AttributePair const& value) {
        if constexpr (Ascending)
            cs_asc_[key].emplace(value);
        else
            cs_desc_[key].emplace(value);
    }

    template <bool Ascending>
    void CSPut(AttributeSet const& key, AttributePair&& value) {
        if constexpr (Ascending)
            cs_asc_[key].emplace(std::move(value));
        else
            cs_desc_[key].emplace(std::move(value));
    }

    template <bool Ascending>
    std::unordered_set<AttributePair>& CSGet(AttributeSet const& key) {
        if constexpr (Ascending)
            return cs_asc_[key];
        else
            return cs_desc_[key];
    }

    void Initialize();
    void ComputeODs();
    void PruneLevels();
    void CalculateNextLevel();

    template <bool Ascending>
    void AddCandidates(AttributeSet const& context, std::vector<AttributeSet> const& delAttrs) {
        if (level_ == 2) {
            for (model::ColumnIndex i = 0; i < data_->GetColumnCount(); i++) {
                for (model::ColumnIndex j = 0; j < data_->GetColumnCount(); j++) {
                    if (i == j) continue;
                    CSPut<Ascending>(fastod::CreateAttributeSet({i, j}, data_->GetColumnCount()),
                                     AttributePair(i, j));
                }
            }
        } else if (level_ > 2) {
            context.Iterate([this, &delAttrs, &context](model::ColumnIndex attr) {
                auto const& candidates = CSGet<Ascending>(delAttrs[attr]);

                for (AttributePair const& attribute_pair : candidates) {
                    const AttributeSet context_delete_ab = fastod::DeleteAttribute(
                            delAttrs[attribute_pair.left], attribute_pair.right);

                    bool add_context = true;

                    context_delete_ab.Iterate([this, &delAttrs, &attribute_pair,
                                               &add_context](model::ColumnIndex attr) {
                        std::unordered_set<AttributePair> const& cs =
                                CSGet<Ascending>(delAttrs[attr]);
                        if (cs.find(attribute_pair) == cs.end()) {
                            add_context = false;
                            return;
                        }
                    });

                    if (add_context) {
                        CSPut<Ascending>(context, attribute_pair);
                    }
                }
            });
        }
    }

    template <bool Ascending>
    void CalculateODs(AttributeSet const& context, std::vector<AttributeSet> const& delAttrs) {
        auto& cs_for_con = CSGet<Ascending>(context);
        for (auto it = cs_for_con.begin(); it != cs_for_con.end();) {
            model::ColumnIndex a = it->left;
            model::ColumnIndex b = it->right;

            if (ContainsAttribute(CCGet(delAttrs[b]), a) &&
                ContainsAttribute(CCGet(delAttrs[a]), b)) {
                fastod::CanonicalOD<Ascending> od(fastod::DeleteAttribute(delAttrs[a], b), a, b);
                if (od.IsValid(data_, partition_cache_)) {
                    ++ocd_count_;
                    AddToResult(std::move(od));
                    cs_for_con.erase(it++);
                } else {
                    ++it;
                }
            } else {
                cs_for_con.erase(it++);
            }
        }
    }

public:
    Fastod();

    void PrintStatistics() const;
    bool IsComplete() const;

    std::vector<AscCanonicalOD> const& GetAscendingDependencies() const;
    std::vector<DescCanonicalOD> const& GetDescendingDependencies() const;
    std::vector<SimpleCanonicalOD> const& GetSimpleDependencies() const;

    std::tuple<std::vector<AscCanonicalOD> const&, std::vector<DescCanonicalOD> const&,
               std::vector<SimpleCanonicalOD> const&>
    Discover();

    std::vector<std::string> DiscoverAsStrings();
};

}  // namespace algos
