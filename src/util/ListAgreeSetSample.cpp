#include "ListAgreeSetSample.h"

std::shared_ptr<ListAgreeSetSample> ListAgreeSetSample::createFocusedFor(std::shared_ptr<ColumnLayoutRelationData> relation,
                                                        std::shared_ptr<Vertical> restrictionVertical,
                                                        std::shared_ptr<PositionListIndex> restrictionPLi,
                                                        unsigned int sampleSize, CustomRandom& random) {
    return AgreeSetSample::createFocusedFor<ListAgreeSetSample>(relation, restrictionVertical, restrictionPLi, sampleSize, random);
}

std::shared_ptr<std::vector<long long>> ListAgreeSetSample::bitSetToLongLongVector(boost::dynamic_bitset<> const& bitset) {
    auto result = std::make_shared<std::vector<long long>>(std::vector<long long>( (bitset.size() + 63 ) / 64, 0));
    for (size_t i = 0; i < bitset.size(); i++) {
        (*result)[i / 64] |= bitset[i] << i % 64;  //idea is: long long ~ 64 bits. shift i-th bit in the bitset i mod 64 times and set the corresponding bit
    }
    return result;
}

ListAgreeSetSample::ListAgreeSetSample(std::shared_ptr<ColumnLayoutRelationData> relation, std::shared_ptr<Vertical> focus,
        unsigned int sampleSize, unsigned long long populationSize,
        std::map<boost::dynamic_bitset<>, int> const & agreeSetCounters) : AgreeSetSample(relation, focus, sampleSize, populationSize) {
    for(auto el : agreeSetCounters) {
        agreeSetCounters_.emplace_back(Entry(bitSetToLongLongVector(el.first), el.second));
    }
}

//TODO: for now, ignoring profiling data (namely, _nanoQueries and _numQueries)
long long ListAgreeSetSample::getNumAgreeSupersets(std::shared_ptr<Vertical> agreement) {
    long long count = 0;
    std::vector<long long> minAgreeSet = *bitSetToLongLongVector(agreement->getColumnIndices());

    for (const auto& agreeSetCounter : agreeSetCounters_) {
        std::vector<long long> agreeSet = *agreeSetCounter.agreeSet_;
        int i = 0;
        int minFields = std::min(agreeSet.size(), minAgreeSet.size());
        while (i < minFields) {
            if ((agreeSet[i] & minAgreeSet[i]) != minAgreeSet[i]) goto Entries;
            i++;
        }
        while (i < minAgreeSet.size()) {
            if (minAgreeSet[i] != 0) goto Entries;
            i++;
        }
        count += agreeSetCounter.count_;
    Entries:
        continue;
    }
    //_numQueries
    //_nanoQueries
    return count;
}

long long ListAgreeSetSample::getNumAgreeSupersets(std::shared_ptr<Vertical> agreement, std::shared_ptr<Vertical> disagreement) {
    long long count = 0;
    std::vector<long long> minAgreeSet = *bitSetToLongLongVector(agreement->getColumnIndices());
    std::vector<long long> minDisagreeSet = *bitSetToLongLongVector(disagreement->getColumnIndices());
    //std::cout << "-----------------------------------\n";
    for (const auto& agreeSetCounter : agreeSetCounters_) {
        /*for (auto const& el : *agreeSetCounter.agreeSet_)
            std::cout << el << ' ';
        std::cout << agreeSetCounter.count_ << "\n";*/
        std::vector<long long> agreeSet = *agreeSetCounter.agreeSet_;
        //check the agreement
        int i = 0;
        int minFields = std::min(agreeSet.size(), minAgreeSet.size());
        while (i < minFields) {
            if ((agreeSet[i] & minAgreeSet[i]) != minAgreeSet[i]) goto Entries;
            i++;
        }
        while (i < minAgreeSet.size()) {
            if (minAgreeSet[i] != 0) goto Entries;
            i++;
        }
        //check the disagreement
        i = 0;
        minFields = std::min(agreeSet.size(), minDisagreeSet.size());
        while (i < minFields) {
            if ((agreeSet[i] & minDisagreeSet[i]) != 0) goto Entries;
            i++;
        }

        count += agreeSetCounter.count_;
    Entries:
        continue;
    }
    //std::cout << '\n';
    //_numQueries
    //_nanoQueries
    return count;
}

std::shared_ptr<std::vector<long long>> ListAgreeSetSample::getNumAgreeSupersetsExt(std::shared_ptr<Vertical> agreement,
                                                                                    std::shared_ptr<Vertical> disagreement) {
    long long count = 0, countAgreements = 0;
    std::vector<long long> minAgreeSet = *bitSetToLongLongVector(agreement->getColumnIndices());
    std::vector<long long> minDisagreeSet = *bitSetToLongLongVector(disagreement->getColumnIndices());

    for (const auto& agreeSetCounter : agreeSetCounters_) {
        std::vector<long long> agreeSet = *agreeSetCounter.agreeSet_;
        //check the agreement
        int i = 0;
        int minFields = std::min(agreeSet.size(), minAgreeSet.size());
        while (i < minFields) {
            if ((agreeSet[i] & minAgreeSet[i]) != minAgreeSet[i]) goto Entries;
            i++;
        }
        while (i < minAgreeSet.size()) {
            if (minAgreeSet[i] != 0) goto Entries;
            i++;
        }
        countAgreements += agreeSetCounter.count_;
        //check the disagreement
        i = 0;
        minFields = std::min(agreeSet.size(), minDisagreeSet.size());
        while (i < minFields) {
            if ((agreeSet[i] & minDisagreeSet[i]) != 0) goto Entries;
            i++;
        }

        count += agreeSetCounter.count_;
        Entries:
        continue;
    }
    //_numQueries
    //_nanoQueries
    return std::make_shared<std::vector<long long>> (std::vector<long long> {countAgreements, count});
}
