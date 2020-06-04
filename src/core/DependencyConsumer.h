#pragma once
#include <functional>
#include "PartialFD.h"
#include "PartialKey.h"
class DependencyConsumer {
protected:
    std::function<void (PartialFD const&)> fdConsumer_;
    std::function<void (PartialKey const&)> uccConsumer_;
public:
    PartialFD registerFd(std::shared_ptr<Vertical> lhs, std::shared_ptr<Column> rhs, double error, double score) const;
    PartialKey registerUcc(std::shared_ptr<Vertical> keyVertical, double error, double score) const;

};