#include "ARAlgorithm.h"

unsigned long long ARAlgorithm::execute() {
    transactionalData = TransactionalData::createFrom(inputGenerator, inputFormat, hasTransactionID);
    if (transactionalData->getNumTransactions() == 0) {
        throw std::runtime_error("Got an empty .csv file: AR mining is meaningless.");
    }

    auto time = findFrequent();
    time += generateAllRules();

    return time;
}
