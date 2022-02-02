#pragma once
#include <functional>
#include <list>
#include <mutex>

#include "PartialFD.h"
#include "PartialKey.h"

class DependencyConsumer {
private:
    std::mutex mutable discover_fd_mutex_;
    std::mutex mutable discover_ucc_mutex_;

    std::list<PartialFD> discoveredFDs_;
    std::list<PartialKey> discoveredUCCs_;

protected:
    std::function<void (PartialFD const&)> fdConsumer_;
    std::function<void (PartialKey const&)> uccConsumer_;

    void discoverFD(PartialFD const& fd) {
        std::scoped_lock lock(discover_fd_mutex_);
        discoveredFDs_.push_back(fd);
    }

    void discoverUCC(PartialKey const& key) {
        std::scoped_lock lock(discover_ucc_mutex_);
        discoveredUCCs_.push_back(key);
    }

public:
    PartialFD registerFd(Vertical const& lhs, Column const& rhs, double error,
                         double score) const;
    PartialKey registerUcc(Vertical const& keyVertical, double error,
                           double score) const;

    std::string fdsToString() const;
    std::string uccsToString() const;

    virtual std::string getJsonFDs();
};
