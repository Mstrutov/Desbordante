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
protected:
    std::unique_ptr<TransactionalData> transactionalData;

    double minsup;

    void generateRules(Itemset const& frequentItemset);  //где применяем minconf?

    virtual unsigned long long generateAllRules() = 0;
    virtual unsigned long long findFrequent() = 0;
public:
    ARAlgorithm(double minsup, double minconf,
                std::filesystem::path const& path,
                char separator = ',',
                bool hasHeader = true)
            : minconf(minconf), inputGenerator(path, separator, hasHeader), minsup(minsup) {}

    std::list<AR> const& arIDsList() const noexcept { return arCollection; }
    std::list<std::vector<std::string>> arList() const; //TODO возвращает айтемы строками

    unsigned long long execute();
    virtual ~ARAlgorithm() = default;
};
