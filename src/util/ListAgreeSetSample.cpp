#include "ListAgreeSetSample.h"
#include "logging/easylogging++.h"

namespace util {

std::unique_ptr<ListAgreeSetSample> ListAgreeSetSample::CreateFocusedFor(
    ColumnLayoutRelationData const* relation, Vertical const& restriction_vertical,
    PositionListIndex const* restriction_p_li, unsigned int sample_size, CustomRandom& random) {
    return AgreeSetSample::CreateFocusedFor<ListAgreeSetSample>(
        relation, restriction_vertical, restriction_p_li, sample_size, random);
}

std::unique_ptr<std::vector<unsigned long long>> ListAgreeSetSample::BitSetToLongLongVector(
    boost::dynamic_bitset<> const& bitset) {
    auto result = std::make_unique<std::vector<unsigned long long>>(
        std::vector<unsigned long long>((bitset.size() + 63) / 64, 0));
    for (size_t i = 0; i < bitset.size(); i++) {
        //idea is: long long ~ 64 bits. shift i-th bit in the bitset i mod 64 times and set the corresponding bit
        (*result)[i / 64] |= bitset[i] << i % 64;
    }
    return result;
}

ListAgreeSetSample::ListAgreeSetSample(
    ColumnLayoutRelationData const* relation, Vertical const& focus, unsigned int sample_size,
    unsigned long long population_size,
    std::unordered_map<boost::dynamic_bitset<>, int> const& agree_set_counters)
    : AgreeSetSample(relation, focus, sample_size, population_size) {
    for (auto& el : agree_set_counters) {
        agree_set_counters_.emplace_back(Entry(BitSetToLongLongVector(el.first), el.second));
    }
}

unsigned long long ListAgreeSetSample::GetNumAgreeSupersets(Vertical const& agreement) const {
    unsigned long long count = 0;
    std::vector<unsigned long long> min_agree_set =
        *BitSetToLongLongVector(agreement.GetColumnIndices());

    for (const auto& agree_set_counter : agree_set_counters_) {
        std::vector<unsigned long long> agree_set = *agree_set_counter.agree_set;
        unsigned int i = 0;
        unsigned int min_fields = std::min(agree_set.size(), min_agree_set.size());
        while (i < min_fields) {
            if ((agree_set[i] & min_agree_set[i]) != min_agree_set[i]) goto Entries;
            i++;
        }
        while (i < min_agree_set.size()) {
            if (min_agree_set[i] != 0) goto Entries;
            i++;
        }
        count += agree_set_counter.count;
    Entries:
        continue;
    }
    return count;
}

unsigned long long ListAgreeSetSample::GetNumAgreeSupersets(Vertical const& agreement,
                                                            Vertical const& disagreement) const {
    unsigned long long count = 0;
    std::vector<unsigned long long> min_agree_set =
        *BitSetToLongLongVector(agreement.GetColumnIndices());
    std::vector<unsigned long long> min_disagree_set =
        *BitSetToLongLongVector(disagreement.GetColumnIndices());
    //std::cout << "-----------------------------------\n";
    for (const auto& agree_set_counter : agree_set_counters_) {
        /*for (auto const& el : *agree_set_counter.agreeSet_)
            std::cout << el << ' ';
        std::cout << agree_set_counter.count_ << "\n";*/
        std::vector<unsigned long long> agree_set = *agree_set_counter.agree_set;
        //check the agreement
        unsigned int i = 0;
        unsigned int min_fields = std::min(agree_set.size(), min_agree_set.size());
        while (i < min_fields) {
            if ((agree_set[i] & min_agree_set[i]) != min_agree_set[i]) goto Entries;
            i++;
        }
        while (i < min_agree_set.size()) {
            if (min_agree_set[i] != 0) goto Entries;
            i++;
        }
        //check the disagreement
        i = 0;
        min_fields = std::min(agree_set.size(), min_disagree_set.size());
        while (i < min_fields) {
            if ((agree_set[i] & min_disagree_set[i]) != 0) goto Entries;
            i++;
        }

        count += agree_set_counter.count;
    Entries:
        continue;
    }
    LOG(DEBUG) << boost::format{"AgreeSetSample for %1% against %2% returned %3% "} %
                      agreement.ToString() % disagreement.ToString() % count;
    //std::cout << '\n';
    //_numQueries
    //_nanoQueries
    return count;
}

std::unique_ptr<std::vector<unsigned long long>> ListAgreeSetSample::GetNumAgreeSupersetsExt(
    Vertical const& agreement, Vertical const& disagreement) const {
    unsigned long long count = 0, count_agreements = 0;
    std::vector<unsigned long long> min_agree_set =
        *BitSetToLongLongVector(agreement.GetColumnIndices());
    std::vector<unsigned long long> min_disagree_set =
        *BitSetToLongLongVector(disagreement.GetColumnIndices());

    for (const auto& agree_set_counter : agree_set_counters_) {
        std::vector<unsigned long long> agree_set = *agree_set_counter.agree_set;
        //check the agreement
        unsigned int i = 0;
        unsigned int min_fields = std::min(agree_set.size(), min_agree_set.size());
        while (i < min_fields) {
            if ((agree_set[i] & min_agree_set[i]) != min_agree_set[i]) goto Entries;
            i++;
        }
        while (i < min_agree_set.size()) {
            if (min_agree_set[i] != 0) goto Entries;
            i++;
        }
        count_agreements += agree_set_counter.count;
        //check the disagreement
        i = 0;
        min_fields = std::min(agree_set.size(), min_disagree_set.size());
        while (i < min_fields) {
            if ((agree_set[i] & min_disagree_set[i]) != 0) goto Entries;
            i++;
        }

        count += agree_set_counter.count;
        Entries:
        continue;
    }
    return std::make_unique<std::vector<unsigned long long>>(
        std::vector<unsigned long long>{count_agreements, count});
}

} // namespace util
