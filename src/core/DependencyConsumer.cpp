
#include "DependencyConsumer.h"

PartialFD DependencyConsumer::registerFd(std::shared_ptr<Vertical> lhs, std::shared_ptr<Column> rhs, double error, double score) {
    PartialFD partialFd(lhs, rhs, error, score);
    fdConsumer(partialFd);
    return partialFd;
}

PartialKey DependencyConsumer::registerUcc(std::shared_ptr<Vertical> keyVertical, double error, double score) {
    PartialKey partialKey(keyVertical, error, score);
    uccConsumer(partialKey);
    return partialKey;
}


