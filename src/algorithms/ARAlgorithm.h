#pragma once

#include <vector>
#include <list>
#include <set>

#include <boost/any.hpp>

#include "TransactionalData.h"
#include "AR.h"
#include "Primitive.h"

class ARAlgorithm : public algos::Primitive {
public:
    struct Config {
        std::filesystem::path data{};   /* Path to input file */
        char separator = ',';           /* Separator for csv */
        bool has_header = true;         /* Indicates if input file has header */
        std::shared_ptr<InputFormat> input_format;
        double minsup = 0;
        double minconf = 0;
    };

private:
    double minconf_;
    std::list<ARStrings> ar_collection_;
    std::shared_ptr<InputFormat> input_format_;

    struct RuleNode {
        ArIDs rule;
        std::list<RuleNode> children;
        RuleNode() = default;
        RuleNode(std::vector<unsigned>&& left, std::vector<unsigned>&& right, double confidence)
            : rule(std::move(left), std::move(right), confidence) {}
    };

    RuleNode root_;

    bool GenerateRuleLevel(std::vector<unsigned> const& frequent_itemset,
                           double support, unsigned level_number);
    bool MergeRules(std::vector<unsigned> const& frequent_itemset, double support, RuleNode* node);
    static void UpdatePath(std::stack<RuleNode*>& path, std::list<RuleNode>& vertices);
protected:
    std::unique_ptr<TransactionalData> transactional_data_;
    double minsup_;

    void GenerateRulesFrom(std::vector<unsigned> const& frequent_itemset, double support);
    void registerARStrings(ArIDs const& rule);

    virtual double GetSupport(std::vector<unsigned> const& frequent_itemset) const = 0;
    virtual unsigned long long GenerateAllRules() = 0;
    virtual unsigned long long FindFrequent() = 0;
public:
    explicit ARAlgorithm(Config const& config, std::vector<std::string_view> phase_names)
        : Primitive(config.data, config.separator, config.has_header, std::move(phase_names)),
          minconf_(config.minconf),
          input_format_(config.input_format),
          minsup_(config.minsup) {}

    std::list<ARStrings> GetArList() const noexcept { return ar_collection_; }
    virtual std::list<std::set<std::string>> GetFrequentList() const = 0;   //for debugging and testing

    unsigned long long Execute() override;
    virtual ~ARAlgorithm() = default;
};
