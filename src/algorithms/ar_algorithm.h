#pragma once

#include <list>
#include <set>
#include <stack>
#include <vector>

#include <boost/any.hpp>

#include "algorithms/algorithm.h"
#include "algorithms/ar_algorithm_enums.h"
#include "model/ar.h"
#include "model/transactional_data.h"
#include "util/config/tabular_data/input_table_type.h"

namespace algos {

class ARAlgorithm : public Algorithm {
private:
    util::config::InputTable input_table_;

    double minconf_;
    InputFormat input_format_ = InputFormat::singular;
    unsigned int tid_column_index_;
    unsigned int item_column_index_;
    bool first_column_tid_;
    std::list<model::ArIDs> ar_collection_;

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

    void ResetState() final;
    virtual void ResetStateAr() = 0;

protected:
    std::unique_ptr<model::TransactionalData> transactional_data_;
    double minsup_;

    void GenerateRulesFrom(std::vector<unsigned> const& frequent_itemset, double support);

    virtual double GetSupport(std::vector<unsigned> const& frequent_itemset) const = 0;
    virtual unsigned long long GenerateAllRules() = 0;
    virtual unsigned long long FindFrequent() = 0;
    void LoadDataInternal() final;
    unsigned long long ExecuteInternal() final;

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
