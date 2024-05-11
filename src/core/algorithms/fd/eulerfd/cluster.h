#include <vector>
#include <array>
#include <functional>
#include <algorithm>
#include <numeric>

namespace algos {

class Cluster {
public:
    using RegisterTuplesFunction = std::function<size_t (size_t, size_t)>;
    using RandomStrategy = std::function<int ()>;

private:
    std::vector<size_t> cluster_data_;

    constexpr static const size_t initial_window_ = 3;
    std::array<double, initial_window_> hist_effects_;
    size_t window_ = 0;

    size_t sample_number_ = 0;
    size_t new_tuples_pairs_ = 0;
    size_t new_non_fds_ = 0;

    void shuffle_data(RandomStrategy rand);

public:
    Cluster(std::vector<size_t> &&cluster_data, RandomStrategy rand);

    Cluster(Cluster &other) = delete;
    Cluster & operator=(Cluster &other) = delete;

    Cluster(Cluster &&other) = default;
    Cluster & operator=(Cluster &&other) = default;

    double Sample(RegisterTuplesFunction handle_tuples);
    double GetAverage() const;
};
}
