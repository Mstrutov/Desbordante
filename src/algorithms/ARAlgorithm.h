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
protected:
    std::unique_ptr<TransactionalData> transactionalData;

    double minsup;

    void generateRules(Itemset const& frequentItemset);  //где применяем minconf?

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
