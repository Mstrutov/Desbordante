#pragma once

#include <list>
#include <set>
#include <stack>
#include <vector>

#include <boost/any.hpp>

#include "AR.h"
#include "Primitive.h"
#include "TransactionalData.h"

namespace algos {

class ARAlgorithm : public algos::Primitive {
public:
    struct Config {
        std::filesystem::path data{}; /* Path to input file */
        char separator = ',';         /* Separator for csv */
        bool has_header = true;       /* Indicates if input file has header */
        std::shared_ptr<model::InputFormat> input_format;
        double minsup = 0;
        double minconf = 0;
    };

private:
    double minconf_;
    std::list<model::ArIDs> ar_collection_;
    std::shared_ptr<model::InputFormat> input_format_;

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

protected:
    std::unique_ptr<model::TransactionalData> transactional_data_;
    double minsup_;

    void GenerateRulesFrom(std::vector<unsigned> const& frequent_itemset, double support);

    virtual double GetSupport(std::vector<unsigned> const& frequent_itemset) const = 0;
    virtual unsigned long long GenerateAllRules() = 0;
    virtual unsigned long long FindFrequent() = 0;

public:
    explicit ARAlgorithm(Config const& config, std::vector<std::string_view> phase_names)
        : Primitive(config.data, config.separator, config.has_header, std::move(phase_names)),
          minconf_(config.minconf),
          input_format_(config.input_format),
          minsup_(config.minsup) {}

    std::list<model::ArIDs> const& GetArIDsList() const noexcept { return ar_collection_; };
    std::vector<std::string> const& GetItemNamesVector() const noexcept {
        return transactional_data_->GetItemUniverse();
    }

    virtual std::list<std::set<std::string>> GetFrequentList() const = 0;  // for debugging and testing
    std::list<model::ARStrings> GetArStringsList() const;

    unsigned long long Execute() override;
    virtual ~ARAlgorithm() = default;
};

} // namespace algos
