#include "algorithms/md/hymd/validator.h"

#include <cassert>
#include <functional>
#include <vector>

#include "algorithms/md/hymd/indexes/records_info.h"
#include "algorithms/md/hymd/lattice/rhs.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/table_identifiers.h"
#include "algorithms/md/hymd/utility/invalidated_rhss.h"
#include "algorithms/md/hymd/utility/java_hash.h"
#include "model/index.h"
#include "util/bitset_utils.h"
#include "util/py_tuple_hash.h"

namespace {
using model::Index, model::md::DecisionBoundary;
using namespace algos::hymd;
using indexes::CompressedRecords;
using indexes::PliCluster;
using indexes::RecSet;
using indexes::SimilarityMatrix;
using utility::InvalidatedRhss;
using RecommendationVector = std::vector<Recommendation>;
using IndexVector = std::vector<Index>;
using AllRecomVecs = std::vector<RecommendationVector>;
using RecIdVec = std::vector<RecordIdentifier>;
using RecPtr = CompressedRecord const*;
using RecordCluster = std::vector<RecPtr>;

template <typename ElementType>
std::vector<ElementType> GetAllocatedVector(std::size_t size) {
    std::vector<ElementType> vec;
    vec.reserve(size);
    return vec;
}
}  // namespace

namespace std {
template <>
struct hash<vector<ValueIdentifier>> {
    size_t operator()(vector<ValueIdentifier> const& p) const {
        constexpr bool use_java_hash = true;
        if constexpr (use_java_hash) {
            return utility::HashIterable(p);
        } else {
            auto hasher = util::PyTupleHash(p.size());
            for (ValueIdentifier el : p) {
                hasher.AddValue(el);
            }
            return hasher.GetResult();
        }
    }
};
}  // namespace std

namespace algos::hymd {

struct WorkingInfo {
    RecommendationVector& recommendations;
    MdElement const old_rhs;
    DecisionBoundary current_bound;
    std::size_t const col_match_values;
    DecisionBoundary interestingness_boundary;
    CompressedRecords const& right_records;
    SimilarityMatrix const& similarity_matrix;
    Index const left_index;
    Index const right_index;

    bool EnoughRecommendations() const {
        return true;
        // <=> return recommendations.size() >= 1 /* was 20 */;
        // I believe this check is no longer needed, as we are only giving "useful" recommendations,
        // which means those that are very likely to actually remove the need to validate MDs.
    }

    bool ShouldStop() const {
        return current_bound == kLowestBound && EnoughRecommendations();
    }

    WorkingInfo(MdElement old_rhs, RecommendationVector& recommendations,
                std::size_t col_match_values, CompressedRecords const& right_records,
                SimilarityMatrix const& similarity_matrix, Index const left_index,
                Index const right_index)
        : recommendations(recommendations),
          old_rhs(old_rhs),
          current_bound(old_rhs.decision_boundary),
          col_match_values(col_match_values),
          right_records(right_records),
          similarity_matrix(similarity_matrix),
          left_index(left_index),
          right_index(right_index) {}
};

RecSet const* Validator::GetSimilarRecords(ValueIdentifier value_id, DecisionBoundary lhs_bound,
                                           Index column_match_index) const {
    assert(lhs_bound != kLowestBound);
    indexes::SimilarityIndex const& similarity_index =
            (*column_matches_info_)[column_match_index].similarity_info.similarity_index;
    indexes::MatchingRecsMapping const& val_index = similarity_index[value_id];
    auto it = val_index.lower_bound(lhs_bound);
    if (it == val_index.end()) return nullptr;
    return &it->second;
}

template <typename PairProvider>
class Validator::SetPairProcessor {
    Validator const* const validator_;
    std::vector<ColumnMatchInfo> const& column_matches_info_ = *validator_->column_matches_info_;
    CompressedRecords const& left_records_ = validator_->GetLeftCompressor().GetRecords();
    CompressedRecords const& right_records_ = validator_->GetRightCompressor().GetRecords();
    InvalidatedRhss& invalidated_;
    lattice::Rhs& lattice_rhs_;
    MdLhs const& lhs_;
    PairProvider pair_provider_;

    enum class Status { kInvalidated, kCheckedAll };

    [[nodiscard]] Status LowerForColumnMatch(WorkingInfo& working_info, PliCluster const& cluster,
                                             RecSet const& similar_records) const;
    [[nodiscard]] Status LowerForColumnMatch(WorkingInfo& working_info,
                                             RecordCluster const& matched_records,
                                             RecIdVec const& similar_records) const;

    template <typename Collection>
    Status LowerForColumnMatchNoCheck(WorkingInfo& working_info,
                                      RecordCluster const& matched_records,
                                      Collection const& similar_records) const;

    std::pair<std::vector<WorkingInfo>, AllRecomVecs> MakeWorkingAndRecs(
            boost::dynamic_bitset<> const& indices_bitset);

    bool Supported(std::size_t support) {
        return validator_->Supported(support);
    }

    Result MakeAllInvalidatedAndSupportedResult(std::vector<WorkingInfo> const& working,
                                                AllRecomVecs&& recommendations) {
        for (WorkingInfo const& working_info : working) {
            invalidated_.PushBack(working_info.old_rhs, kLowestBound);
        }
        return {std::move(recommendations), std::move(invalidated_), false};
    }

    Result MakeOutOfClustersResult(std::vector<WorkingInfo> const& working,
                                   AllRecomVecs&& recommendations, std::size_t support) {
        for (WorkingInfo const& working_info : working) {
            DecisionBoundary const new_bound = working_info.current_bound;
            MdElement old_rhs = working_info.old_rhs;
            if (new_bound == old_rhs.decision_boundary) continue;
            invalidated_.PushBack(old_rhs, new_bound);
        }
        return {std::move(recommendations), std::move(invalidated_), !Supported(support)};
    }

public:
    SetPairProcessor(Validator const* validator, InvalidatedRhss& invalidated,
                     lattice::Rhs& rhs, MdLhs const& lhs)
        : validator_(validator),
          invalidated_(invalidated),
          lattice_rhs_(rhs),
          lhs_(lhs),
          pair_provider_(validator, lhs) {}

    Result ProcessPairs(boost::dynamic_bitset<> const& indices_bitset) {
        auto [working, recommendations] = MakeWorkingAndRecs(indices_bitset);
        std::size_t support = 0;
        while (pair_provider_.TryGetNextPair()) {
            auto const& cluster = pair_provider_.GetCluster();
            auto const& similar = pair_provider_.GetSimilarRecords();
            support += cluster.size() * similar.size();
            bool all_invalid = true;
            for (WorkingInfo& working_info : working) {
                Status const status = LowerForColumnMatch(working_info, cluster, similar);
                if (status == Status::kCheckedAll) all_invalid = false;
            }
            if (all_invalid && Supported(support))
                return MakeAllInvalidatedAndSupportedResult(working, std::move(recommendations));
        }
        return MakeOutOfClustersResult(working, std::move(recommendations), support);
    }
};

template <typename PairProvider>
std::pair<std::vector<WorkingInfo>, AllRecomVecs>
Validator::SetPairProcessor<PairProvider>::MakeWorkingAndRecs(
        boost::dynamic_bitset<> const& indices_bitset) {
    std::pair<std::vector<WorkingInfo>, AllRecomVecs> working_and_recs;
    auto& [working, recommendations] = working_and_recs;
    std::size_t const working_size = indices_bitset.count();
    working.reserve(working_size);
    recommendations.reserve(working_size);
    IndexVector indices = util::BitsetToIndices<Index>(indices_bitset);
    for (Index index : indices) {
        RecommendationVector& last_recs = recommendations.emplace_back();
        auto const& [sim_info, left_index, right_index] = column_matches_info_[index];
        MdElement rhs{index, lattice_rhs_[index]};
        working.emplace_back(rhs, last_recs, validator_->GetLeftValueNum(index), right_records_,
                             sim_info.similarity_matrix, left_index, right_index);
    }

    auto for_each_working = [&](auto f) { std::for_each(working.begin(), working.end(), f); };
    for_each_working([&](WorkingInfo const& w) { lattice_rhs_[w.old_rhs.index] = kLowestBound; });
    std::vector<DecisionBoundary> const gen_max_rhs =
            validator_->lattice_->GetRhsInterestingnessBounds(lhs_, indices);
    for_each_working([&](WorkingInfo const& w) {
        lattice_rhs_[w.old_rhs.index] = w.old_rhs.decision_boundary;
    });

    auto it = working.begin();
    auto set_advance = [&it](DecisionBoundary bound) { it++->interestingness_boundary = bound; };
    std::for_each(gen_max_rhs.begin(), gen_max_rhs.end(), set_advance);
    return working_and_recs;
}

template <typename PairProvider>
template <typename Collection>
auto Validator::SetPairProcessor<PairProvider>::LowerForColumnMatchNoCheck(
        WorkingInfo& working_info, RecordCluster const& matched_records,
        Collection const& similar_records) const -> Status {
    assert(!similar_records.empty());
    assert(!matched_records.empty());

    std::unordered_map<ValueIdentifier, RecordCluster> grouped(
            std::min(matched_records.size(), working_info.col_match_values));
    for (RecPtr left_record_ptr : matched_records) {
        grouped[(*left_record_ptr)[working_info.left_index]].push_back(left_record_ptr);
    }
    DecisionBoundary& current_rhs_bound = working_info.current_bound;
    SimilarityMatrix const& similarity_matrix = working_info.similarity_matrix;
    for (auto const& [left_value_id, records_left] : grouped) {
        for (RecordIdentifier record_id_right : similar_records) {
            CompressedRecord const& right_record = working_info.right_records[record_id_right];
            auto add_recommendations = [&records_left, &right_record, &working_info]() {
                for (RecPtr left_record_ptr : records_left) {
                    working_info.recommendations.emplace_back(left_record_ptr, &right_record);
                }
            };
            auto const& row = similarity_matrix[left_value_id];
            ValueIdentifier const right_value_id = right_record[working_info.right_index];
            auto it_right = row.find(right_value_id);
            if (it_right == row.end()) {
            rhs_not_valid:
                add_recommendations();
            // rhs_not_valid:
                current_rhs_bound = kLowestBound;
                if (working_info.EnoughRecommendations()) return Status::kInvalidated;
                continue;
            }

            preprocessing::Similarity const pair_similarity = it_right->second;
            // NOTE: I believe the purpose of inference from record pairs is to cut down on the
            // number of costly validations. But a validation still happens if the value was lowered
            // but not to zero.
            // if (pair_similarity < working_info.old_rhs.decision_boundary) add_recommendations();
            if (pair_similarity < current_rhs_bound) current_rhs_bound = pair_similarity;
            if (current_rhs_bound <= working_info.interestingness_boundary) goto rhs_not_valid;
        }
    }
    return Status::kCheckedAll;
}

template <typename PairProvider>
auto Validator::SetPairProcessor<PairProvider>::LowerForColumnMatch(
        WorkingInfo& working_info, RecordCluster const& matched_records,
        RecIdVec const& similar_records) const -> Status {
    if (working_info.ShouldStop()) return Status::kInvalidated;
    return LowerForColumnMatchNoCheck(working_info, matched_records, similar_records);
}

template <typename PairProvider>
auto Validator::SetPairProcessor<PairProvider>::LowerForColumnMatch(
        WorkingInfo& working_info, PliCluster const& cluster, RecSet const& similar_records) const
        -> Status {
    if (working_info.ShouldStop()) return Status::kInvalidated;

    assert(!similar_records.empty());
    RecordCluster cluster_records = GetAllocatedVector<RecPtr>(cluster.size());
    for (RecordIdentifier left_record_id : cluster) {
        cluster_records.push_back(&left_records_[left_record_id]);
    }
    return LowerForColumnMatchNoCheck(working_info, cluster_records, similar_records);
}

class Validator::OneCardPairProvider {
    Validator const* const validator_;
    ValueIdentifier value_id_ = ValueIdentifier(-1);
    Index const non_zero_index_;
    DecisionBoundary const decision_boundary_;
    std::vector<PliCluster> const& clusters_ =
            validator_->GetLeftCompressor()
                    .GetPli(validator_->GetLeftPliIndex(non_zero_index_))
                    .GetClusters();
    std::size_t const clusters_size_ = clusters_.size();
    RecSet const* similar_records_ptr_{};

public:
    OneCardPairProvider(Validator const* validator, MdLhs const& lhs)
        : validator_(validator),
          non_zero_index_(lhs.begin()->child_array_index),
          decision_boundary_(lhs.begin()->decision_boundary) {}

    bool TryGetNextPair() {
        for (++value_id_; value_id_ != clusters_size_; ++value_id_) {
            similar_records_ptr_ =
                    validator_->GetSimilarRecords(value_id_, decision_boundary_, non_zero_index_);
            if (similar_records_ptr_ != nullptr) return true;
        }
        return false;
    }

    PliCluster const& GetCluster() const {
        return clusters_[value_id_];
    }

    RecSet const& GetSimilarRecords() const {
        return *similar_records_ptr_;
    }
};

class Validator::MultiCardPairProvider {
    using MatchingInfo = std::vector<std::tuple<Index, Index, DecisionBoundary>>;

    struct InitInfo {
        Validator const* validator;
        MatchingInfo matching_info;
        IndexVector non_first_indices;
        Index first_pli_index;
        std::size_t plis_involved = 1;

        InitInfo(Validator const* validator, MdLhs const& lhs) : validator(validator) {
            using AbsoluteLhsElements = std::vector<MdElement>;
            std::size_t const cardinality = lhs.Cardinality();
            matching_info.reserve(cardinality);

            std::size_t const left_pli_number = validator->GetLeftCompressor().GetPliNumber();
            non_first_indices.reserve(std::min(cardinality, left_pli_number));
            std::vector<AbsoluteLhsElements> pli_map(left_pli_number);
            Index col_match_index = 0;
            for (auto const& [child_array_index, bound] : lhs) {
                col_match_index += child_array_index;
                pli_map[validator->GetLeftPliIndex(col_match_index)].emplace_back(col_match_index,
                                                                                  bound);
                ++col_match_index;
            }

            Index pli_idx = 0;
            while (pli_map[pli_idx].empty()) ++pli_idx;

            Index value_ids_index = 0;
            auto fill_for_value_ids_idx = [this,
                                           &value_ids_index](AbsoluteLhsElements const& elements) {
                for (auto const& [col_match_idx, bound] : elements) {
                    matching_info.emplace_back(col_match_idx, value_ids_index, bound);
                }
                ++value_ids_index;
            };
            first_pli_index = pli_idx;
            fill_for_value_ids_idx(pli_map[pli_idx]);
            for (++pli_idx; pli_idx != left_pli_number; ++pli_idx) {
                AbsoluteLhsElements const& pli_col_matches = pli_map[pli_idx];
                if (pli_col_matches.empty()) continue;
                ++plis_involved;
                non_first_indices.push_back(pli_idx);
                fill_for_value_ids_idx(pli_col_matches);
            }
        }
    };

    using GroupMap = std::unordered_map<std::vector<ValueIdentifier>, RecordCluster>;

    Validator const* const validator_;
    GroupMap grouped_;
    GroupMap::iterator cur_group_iter_ = grouped_.begin();
    GroupMap::iterator end_group_iter_ = grouped_.end();
    ValueIdentifier first_value_id_ = ValueIdentifier(-1);
    std::vector<ValueIdentifier> value_ids_;
    std::vector<RecSet const*> rec_sets_;
    IndexVector const non_first_indices_;
    IndexVector::const_iterator non_first_start_ = non_first_indices_.begin();
    IndexVector::const_iterator non_first_end_ = non_first_indices_.end();
    std::vector<PliCluster> const& first_pli_;
    std::size_t first_pli_size_ = first_pli_.size();
    CompressedRecords const& left_records_ = validator_->GetLeftCompressor().GetRecords();
    MatchingInfo const matching_info_;
    RecIdVec similar_records_;
    RecordCluster const* cluster_ptr_;

    MultiCardPairProvider(InitInfo init_info)
        : validator_(init_info.validator),
          non_first_indices_(std::move(init_info.non_first_indices)),
          first_pli_(
                  validator_->GetLeftCompressor().GetPli(init_info.first_pli_index).GetClusters()),
          matching_info_(std::move(init_info.matching_info)) {
        value_ids_.reserve(init_info.plis_involved);
        rec_sets_.reserve(matching_info_.size());
    }

    bool TryGetNextGroup() {
        if (++first_value_id_ == first_pli_size_) return false;
        grouped_.clear();
        PliCluster const& cluster = first_pli_[first_value_id_];
        value_ids_.push_back(first_value_id_);
        for (RecordIdentifier record_id : cluster) {
            std::vector<ValueIdentifier> const& record = left_records_[record_id];
            for (auto ind_it = non_first_start_; ind_it != non_first_end_; ++ind_it) {
                value_ids_.push_back(record[*ind_it]);
            }
            grouped_[value_ids_].push_back(&record);
            value_ids_.erase(++value_ids_.begin(), value_ids_.end());
        }
        value_ids_.clear();
        cur_group_iter_ = grouped_.begin();
        end_group_iter_ = grouped_.end();
        return true;
    }

public:
    MultiCardPairProvider(Validator const* validator, MdLhs const& lhs)
        : MultiCardPairProvider(InitInfo{validator, lhs}) {}

    bool TryGetNextPair() {
        similar_records_.clear();
        do {
            while (cur_group_iter_ != end_group_iter_) {
                auto const& [val_ids, cluster] = *cur_group_iter_;
                rec_sets_.clear();
                ++cur_group_iter_;
                for (auto const& [column_match_index, value_ids_index, bound] : matching_info_) {
                    RecSet const* similar_records_ptr = validator_->GetSimilarRecords(
                            val_ids[value_ids_index], bound, column_match_index);
                    if (similar_records_ptr == nullptr) goto no_similar_records;
                    rec_sets_.push_back(similar_records_ptr);
                }
                goto matched_on_all;
            no_similar_records:
                continue;
            matched_on_all:
                auto check_set_begin = rec_sets_.begin();
                auto check_set_end = rec_sets_.end();
                using CRecSet = RecSet const;
                auto size_cmp = [](CRecSet* p1, CRecSet* p2) { return p1->size() < p2->size(); };
                std::sort(check_set_begin, check_set_end, size_cmp);
                CRecSet& first = **check_set_begin;
                ++check_set_begin;
                for (RecordIdentifier rec : first) {
                    auto rec_cont = [rec](CRecSet* set_ptr) { return set_ptr->contains(rec); };
                    if (std::all_of(check_set_begin, check_set_end, rec_cont)) {
                        similar_records_.push_back(rec);
                    }
                }
                if (similar_records_.empty()) continue;
                cluster_ptr_ = &cluster;
                return true;
            }
        } while (TryGetNextGroup());
        return false;
    }

    RecordCluster const& GetCluster() const {
        return *cluster_ptr_;
    }

    RecIdVec const& GetSimilarRecords() const {
        return similar_records_;
    }
};

Validator::Result Validator::Validate(lattice::ValidationInfo& info) const {
    MdLhs const& lhs = info.messenger->GetLhs();
    lattice::Rhs& lattice_rhs = info.messenger->GetRhs();
    // After a call to this method, info.rhs_indices must not be used
    boost::dynamic_bitset<>& indices_bitset = info.rhs_indices;
    std::size_t const cardinality = lhs.Cardinality();
    InvalidatedRhss invalidated;
    if (cardinality == 0) [[unlikely]] {
        util::ForEachIndex(indices_bitset, [&](auto index) {
            DecisionBoundary const old_bound = lattice_rhs[index];
            DecisionBoundary const new_bound =
                    (*column_matches_info_)[index].similarity_info.lowest_similarity;
            if (old_bound == new_bound) [[unlikely]]
                return;
            invalidated.PushBack({index, old_bound}, new_bound);
        });
        return {{}, std::move(invalidated), !Supported(GetTotalPairsNum())};
    }

    if (cardinality == 1) {
        Index const non_zero_index = lhs.begin()->child_array_index;
        // Never happens when disjointedness pruning is on.
        if (indices_bitset.test_set(non_zero_index, false)) {
            invalidated.PushBack({non_zero_index, lattice_rhs[non_zero_index]}, kLowestBound);
        }
        SetPairProcessor<OneCardPairProvider> processor(this, invalidated, lattice_rhs, lhs);
        return processor.ProcessPairs(indices_bitset);
    }

    SetPairProcessor<MultiCardPairProvider> processor(this, invalidated, lattice_rhs, lhs);
    return processor.ProcessPairs(indices_bitset);
}

}  // namespace algos::hymd
