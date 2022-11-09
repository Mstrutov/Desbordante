#pragma once

#include <list>
#include <set>
#include <stack>
#include <vector>

#include <boost/any.hpp>

#include "algorithms/ar_algorithm_enums.h"
#include "algorithms/options/type.h"
#include "algorithms/primitive.h"
#include "model/ar.h"
#include "model/transactional_data.h"

namespace algos {

class ARAlgorithm : public algos::Primitive {
protected:
    using MinSupType = double;
    double minconf_;
    InputFormat input_format_ = InputFormat::singular;
    unsigned int tid_column_index_;
    unsigned int item_column_index_;
    bool first_column_tid_;
    std::list<model::ArIDs> ar_collection_;

    static const config::OptionType<decltype(input_format_)> InputFormatOpt;
    static const config::OptionType<decltype(tid_column_index_)> TidColumnIndexOpt;
    static const config::OptionType<decltype(item_column_index_)> ItemColumnIndexOpt;
    static const config::OptionType<decltype(first_column_tid_)> FirstColumnTidOpt;
    static const config::OptionType<MinSupType> MinSupportOpt;
    static const config::OptionType<decltype(minconf_)> MinConfidenceOpt;

    struct RuleNode {
        model::ArIDs rule;
        std::list<RuleNode> children;
        RuleNode() = default;

        /* Temporary fix. Now we allocate generated AR twice -- in ar_collection_
         * and also in a rule node by copying it.
         * */
        explicit RuleNode(model::ArIDs const& rule)
            : rule(rule) {}
    };

    RuleNode root_;

    bool GenerateRuleLevel(std::vector<unsigned> const& frequent_itemset,
                           double support, unsigned level_number);
    bool MergeRules(std::vector<unsigned> const& frequent_itemset, double support, RuleNode* node);
    static void UpdatePath(std::stack<RuleNode*>& path, std::list<RuleNode>& vertices);
    void RegisterOptions();

protected:
    std::unique_ptr<model::TransactionalData> transactional_data_;
    MinSupType minsup_;

    void GenerateRulesFrom(std::vector<unsigned> const& frequent_itemset, double support);

    virtual double GetSupport(std::vector<unsigned> const& frequent_itemset) const = 0;
    virtual unsigned long long GenerateAllRules() = 0;
    virtual unsigned long long FindFrequent() = 0;
    void FitInternal(model::IDatasetStream &data_stream) override;
    void MakeExecuteOptsAvailable() override;
    unsigned long long ExecuteInternal() override;

public:
    explicit ARAlgorithm(std::vector<std::string_view> phase_names);

    std::list<model::ArIDs> const& GetArIDsList() const noexcept { return ar_collection_; };
    std::vector<std::string> const& GetItemNamesVector() const noexcept {
        return transactional_data_->GetItemUniverse();
    }

    virtual std::list<std::set<std::string>> GetFrequentList() const = 0;  // for debugging and testing
    std::list<model::ARStrings> GetArStringsList() const;

    virtual ~ARAlgorithm() = default;
};

}  // namespace algos
