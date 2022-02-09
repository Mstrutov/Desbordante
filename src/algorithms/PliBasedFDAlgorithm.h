#pragma once

#include "ColumnLayoutRelationData.h"
#include "FDAlgorithm.h"

class PliBasedFDAlgorithm : public FDAlgorithm {
private:
    void Initialize() override;

protected:
    std::unique_ptr<ColumnLayoutRelationData> relation_;

    ColumnLayoutRelationData const& GetRelation() const noexcept {
        // GetRelation should be called after input file is parsed i.e. after algorithm execution
        assert(relation_ != nullptr);
        return *relation_;
    }

public:
    explicit PliBasedFDAlgorithm(std::filesystem::path const& path, char separator = ',',
                                 bool has_header = true, bool const is_null_equal_null = true,
                                 std::vector<std::string_view> phase_names = {"FD mining"})
        : FDAlgorithm(path, separator, has_header, is_null_equal_null, std::move(phase_names)) {}
    std::vector<Column const*> GetKeys() const override;
};
