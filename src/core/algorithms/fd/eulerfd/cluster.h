#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <vector>

namespace algos {

class Cluster {
public:
    using RegisterTuplesFunction = std::function<size_t(size_t, size_t)>;
    using RandomStrategy = std::function<int()>;

private:
    std::vector<size_t> cluster_data_;

    constexpr static size_t kInitialWindow = 3;
    std::array<double, kInitialWindow> hist_effects_;
    size_t window_ = 0;

    size_t sample_number_ = 0;
    size_t new_tuples_pairs_ = 0;
    size_t new_non_fds_ = 0;

    void ShuffleData(RandomStrategy rand);

public:
    Cluster(std::vector<size_t> &&cluster_data, RandomStrategy rand);

    Cluster(Cluster &other) = delete;
    Cluster &operator=(Cluster &other) = delete;

    Cluster(Cluster &&other) = default;
    Cluster &operator=(Cluster &&other) = default;

    double Sample(RegisterTuplesFunction handle_tuples);
    double GetAverage() const;
};
}  // namespace algos
