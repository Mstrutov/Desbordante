#pragma once
#include <functional>
#include <list>
#include "PartialFD.h"
#include "PartialKey.h"

class DependencyConsumer {
protected:
    std::function<void (PartialFD const&)> fdConsumer_;
    std::function<void (PartialKey const&)> uccConsumer_;
public:
    std::list<PartialFD> discoveredFDs_;
    std::list<PartialKey> discoveredUCCs_;

    PartialFD registerFd(Vertical const& lhs, Column const& rhs, double error, double score) const;
    PartialKey registerUcc(Vertical const& keyVertical, double error, double score) const;

    std::string fdsToString() const;
    std::string uccsToString() const;

    virtual std::string getJsonFDs();
};