#include "cluster.h"

namespace algos {

void Cluster::ShuffleData(RandomStrategy rand) {
    if (cluster_data_.size() == 0) {
        return;
    }
    for (size_t i = cluster_data_.size() - 1; i > 0; i--) {
        std::swap(cluster_data_[i], cluster_data_[rand() % (i + 1)]);
    }
}

Cluster::Cluster(std::vector<size_t> &&cluster_data, RandomStrategy rand)
    : cluster_data_(std::move(cluster_data)) {
    ShuffleData(rand);
    hist_effects_.fill(1);
}

double Cluster::GetAverage() const {
    double sum = std::accumulate(hist_effects_.begin(), hist_effects_.end(), 0.0);
    return sum / kInitialWindow;
}

double Cluster::Sample(RegisterTuplesFunction handle_tuples) {
    new_tuples_pairs_ = new_non_fds_ = 0;
    window_++;

    for (int i = 0; i < (int)cluster_data_.size() - (int)window_; i++) {
        new_non_fds_ += handle_tuples(cluster_data_[i], cluster_data_[i + window_]);
    }
    new_tuples_pairs_ = cluster_data_.size() - window_;

    double cur_eff = new_tuples_pairs_ == 0 ? 0 : (double)new_non_fds_ / new_tuples_pairs_;
    hist_effects_[sample_number_ % kInitialWindow] = cur_eff;

    sample_number_++;

    return cur_eff;
}

}  // namespace algos
