#pragma once

#include <vector>
#include <list>

#include "TransactionalData.h"
#include "AR.h"

class ARAlgorithm {
private:
    double minconf;
    std::list<AR> arCollection;
    CSVParser inputGenerator;
    TransactionalInputFormat inputFormat;
    bool hasTransactionID;

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
    void registerAR(std::vector<unsigned>&& left, std::vector<unsigned>&& right,
                    double confidence);

    virtual double getSupport(std::vector<unsigned> const& frequentItemset) const = 0;
    virtual unsigned long long generateAllRules() = 0;
    virtual unsigned long long findFrequent() = 0;
public:
    ARAlgorithm(double minsup, double minconf,
                std::filesystem::path const& path,
                TransactionalInputFormat inputFormat = TransactionalInputFormat::TwoColumns,
                bool hasTransactionID = false,
                char separator = ',',
                bool hasHeader = true)
            : minconf(minconf), inputGenerator(path, separator, hasHeader),
              inputFormat(inputFormat), hasTransactionID(hasTransactionID), minsup(minsup) {}

    std::list<AR> const& arIDsList() const noexcept { return arCollection; }
    std::list<std::vector<std::string>> arList() const; //TODO возвращает айтемы строками
    virtual std::list<std::set<std::string>> getAllFrequent() const = 0;   //for debugging and tests

    unsigned long long execute();
    virtual ~ARAlgorithm() = default;
};
