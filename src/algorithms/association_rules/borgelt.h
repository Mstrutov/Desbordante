#pragma once

#include "algorithms/ar_algorithm.h"

namespace algos {

class Borgelt : public ARAlgorithm {
public:
    std::string algo_ = "fpgrowth";
    std::string sep_ = ",";
    std::filesystem::path path_;
    std::filesystem::path const kOutputFile = "tmp";
    static const config::OptionType<decltype(algo_)> AlgoOpt;

private:
    double GetSupport(const std::vector<unsigned int> &frequent_itemset) const override { return frequent_itemset.size(); }
    unsigned long long GenerateAllRules() override { return 0; }
    unsigned long long FindFrequent() override { return 0; }
    std::list<std::set<std::string>> GetFrequentList() const override { return {}; }

    void ConvertAndFillResults();
    void FitInternal(model::IDatasetStream &data_stream) override;
    void MakeExecuteOptsAvailable() override;
    unsigned long long ExecuteInternal() override;
public:
    explicit Borgelt()
        : ARAlgorithm({}) {
        RegisterOption(AlgoOpt.GetOption(&algo_)); }

    ~Borgelt() override;
};

}  // namespace algos
