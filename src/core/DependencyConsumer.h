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

    std::list<PartialFD> discovered_fds_;
    std::list<PartialKey> discovered_uccs_;

protected:
    std::function<void (PartialFD const&)> fd_consumer_;
    std::function<void (PartialKey const&)> ucc_consumer_;

    void DiscoverFd(PartialFD const& fd) {
        std::scoped_lock lock(discover_fd_mutex_);
        discovered_fds_.push_back(fd);
    }

    void DiscoverUcc(PartialKey const& key) {
        std::scoped_lock lock(discover_ucc_mutex_);
        discovered_uccs_.push_back(key);
    }

public:
    PartialFD RegisterFd(Vertical const& lhs, Column const& rhs,
                                  double error, double score) const;
    PartialKey RegisterUcc(Vertical const& key_vertical, double error,
                           double score) const;

    std::string FDsToString() const;
    std::string UCCsToString() const;

    virtual std::string GetJsonFDs();
};
