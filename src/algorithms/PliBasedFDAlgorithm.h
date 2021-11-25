#pragma once

#include "ColumnLayoutRelationData.h"
#include "FDAlgorithm.h"

class PliBasedFDAlgorithm : public FDAlgorithm {
private:
    void initialize() override;

protected:
    std::unique_ptr<ColumnLayoutRelationData> relation_;

    ColumnLayoutRelationData const& getRelation() const noexcept {
        // getRelation should be called after input file is parsed i.e. after algorithm execution
        assert(relation_ != nullptr);
        return *relation_;
    }

public:
    explicit PliBasedFDAlgorithm(std::filesystem::path const& path,
                                 char separator = ',', bool hasHeader = true,
                                 bool const is_null_equal_null = true,
                                 std::vector<std::string_view> phase_names = { "FD mining" }) :
                                 FDAlgorithm(path, separator, hasHeader, is_null_equal_null, std::move(phase_names)) {}
};
