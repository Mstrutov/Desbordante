#pragma once
#include <functional>
#include "model/PartialFD.h"
#include "model/PartialKey.h"
class DependencyConsumer {
protected:
    std::function<void (PartialFD const&)> fdConsumer;
    std::function<void (PartialKey const&)> uccConsumer;
public:
    PartialFD registerFd(std::shared_ptr<Vertical> lhs, std::shared_ptr<Column> rhs, double error, double score) const;
    PartialKey registerUcc(std::shared_ptr<Vertical> keyVertical, double error, double score) const;

};