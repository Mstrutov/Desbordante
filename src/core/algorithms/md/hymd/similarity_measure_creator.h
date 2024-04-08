#pragma once

#include <string>

#include "algorithms/md/hymd/preprocessing/similarity_measure/similarity_measure.h"

namespace algos::hymd {

class SimilarityMeasureCreator {
private:
    std::string const similarity_measure_name_;

public:
    SimilarityMeasureCreator(std::string similarity_measure_name) noexcept
        : similarity_measure_name_(std::move(similarity_measure_name)) {}

    std::string const& GetSimilarityMeasureName() const noexcept {
        return similarity_measure_name_;
    }

    virtual std::unique_ptr<preprocessing::similarity_measure::SimilarityMeasure> MakeMeasure()
            const = 0;

    virtual ~SimilarityMeasureCreator() = default;
};

}  // namespace algos::hymd
