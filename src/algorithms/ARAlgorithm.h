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
    double minconf;
    std::list<ARStrings> arCollection;
    std::shared_ptr<InputFormat> inputFormat;

    struct RuleNode {
        AR rule;
        std::list<RuleNode> children;
        RuleNode() = default;
        RuleNode(std::vector<unsigned>&& left, std::vector<unsigned>&& right, double confidence)
            : rule(std::move(left), std::move(right), confidence) {}
    };

    RuleNode root;

    bool generateRuleLevel(std::vector<unsigned> const& frequentItemset, double support, unsigned levelNumber);
    bool mergeRules(std::vector<unsigned> const& frequentItemset, double support, RuleNode* node);
    static void updatePath(std::stack<RuleNode*> & path, std::list<RuleNode> & vertices);
protected:
    std::unique_ptr<TransactionalData> transactionalData;
    double minsup;

    void generateRulesFrom(std::vector<unsigned> const& frequentItemset, double support);
    void registerARStrings(AR const& rule);

    virtual double getSupport(std::vector<unsigned> const& frequentItemset) const = 0;
    virtual unsigned long long generateAllRules() = 0;
    virtual unsigned long long findFrequent() = 0;
public:
    explicit ARAlgorithm(Config const& config, std::vector<std::string_view> phase_names)
        : Primitive(config.data, config.separator, config.has_header, std::move(phase_names)),
          minconf(config.minconf), inputFormat(config.input_format), minsup(config.minsup) {}

    std::list<ARStrings> arList() const noexcept { return arCollection; }
    virtual std::list<std::set<std::string>> getAllFrequent() const = 0;   //for debugging and testing

    unsigned long long Execute() override;
    virtual ~ARAlgorithm() = default;
};
