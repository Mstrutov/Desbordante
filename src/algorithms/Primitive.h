#pragma once

#include <filesystem>
#include <mutex>
#include <string_view>
#include <vector>

#include "CSVParser.h"

namespace algos {

class Primitive {
private:
    std::mutex mutable progress_mutex_;
    double cur_phase_progress_ = 0;
    uint8_t cur_phase_id_ = 0;

protected:
    CSVParser input_generator_;
    /* Vector of names of algorithm phases, should be initialized in a constructor
     * if algorithm has more than one phase. This vector is used to determine the
     * total number of phases.
     * Use empty vector to intialize this field if your algorithm does not have
     * implemented progress bar.
     */
    std::vector<std::string_view> const phase_names_;

    void AddProgress(double const val) noexcept;
    void SetProgress(double const val) noexcept;
    void ToNextProgressPhase() noexcept;

public:
    constexpr static double kTotalProgressPercent = 100.0;

    Primitive(Primitive const& other) = delete;
    Primitive& operator=(Primitive const& other) = delete;
    Primitive(Primitive&& other) = delete;
    Primitive& operator=(Primitive&& other) = delete;
    virtual ~Primitive() = default;

    explicit Primitive(std::vector<std::string_view> phase_names)
        : phase_names_(std::move(phase_names)) {}
    Primitive(CSVParser input_generator, std::vector<std::string_view> phase_names)
        : input_generator_(std::move(input_generator)), phase_names_(std::move(phase_names)) {}
    Primitive(std::filesystem::path const& path, char const separator, bool const has_header,
              std::vector<std::string_view> phase_names)
        : input_generator_(path, separator, has_header), phase_names_(std::move(phase_names)) {}

    virtual unsigned long long Execute() = 0;

    /* Returns pair with current progress state.
     * Pair has the form <current phase id, current phase progess>
     */
    std::pair<uint8_t, double> GetProgress() const noexcept;

    std::vector<std::string_view> const& GetPhaseNames() const noexcept {
        return phase_names_;
    }
};

}  // namespace algos
