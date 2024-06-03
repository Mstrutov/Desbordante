#pragma once

#include <unordered_map>

#include "types.h"

namespace algos {
class Statistic {
    bool has_value_ = false;
    std::byte const *data_ = nullptr;
    std::unique_ptr<model::Type> type_ = nullptr;

public:
    Statistic() noexcept = default;
    Statistic(std::byte const *, model::Type const *, bool clone_data);
    Statistic(std::byte const *data, model::TypeId id, bool clone_data = false);
    Statistic(Statistic const &);
    Statistic(Statistic &&);
    ~Statistic();

    Statistic &operator=(Statistic const &);
    Statistic &operator=(Statistic &&);

    bool HasValue() const noexcept;
    std::byte const *GetData() const;
    model::Type const *GetType() const;
    std::byte const *ReleaseData();
    std::string ToString() const;
};

struct ColumnStats {
    std::string type;
    size_t count;
    size_t distinct;
    bool is_categorical;
    Statistic avg, STD, skewness, kurtosis, min, max, sum, quantile25, quantile50, quantile75,
            num_zeros, num_negatives, sum_of_squares, geometric_mean, mean_ad, median, median_ad,
            vocab, num_non_letter_chars, num_digit_chars, num_lowercase_chars, num_uppercase_chars,
            num_chars, num_avg_chars, min_num_chars, max_num_chars, min_num_words, max_num_words,
            num_words, num_entirely_uppercase, num_entirely_lowercase;

    std::string ToString() const;
    std::unordered_map<std::string, std::string> ToKeyValueMap() const;
};

}  // namespace algos
