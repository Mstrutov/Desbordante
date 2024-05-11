#pragma once

#include <queue>
#include <memory>
#include <vector>
#include <algorithm>

#include "cluster.h"

namespace algos {

class MLFQ {
private:
    using Queue = std::pair<std::queue<Cluster *>, double>;

    struct LastQueueElement {
        Cluster *cluster {};
        bool operator<(const LastQueueElement &other) const;
    };

    std::vector<Queue> queues_;
    size_t effective_size_ = 0;
    std::priority_queue<LastQueueElement> last_queue_ {};

    // actual_queue is index of not empty queue with the highest prioritet
    int actual_queue_ = -1;

public:
    MLFQ(size_t queues_number);

    void Add(Cluster *cluster, double prioritet, bool add_if_zero = false);
    void AddAtLast(Cluster *cluster);
    [[nodiscard]] Cluster * Get();
    [[nodiscard]] Cluster * GetFromLast();

    size_t GetEffectiveSize() const;
    size_t GetLastQueueSize() const;
    double MaxEffectInLastQueue() const;
    void clear();
};
}
