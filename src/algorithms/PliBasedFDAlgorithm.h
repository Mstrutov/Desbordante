#pragma once

#include "ColumnLayoutRelationData.h"
#include "FDAlgorithm.h"

class PliBasedFDAlgorithm : public FDAlgorithm {
private:
    using algos::Primitive::input_generator_;
    using FDAlgorithm::is_null_equal_null_;

    void Initialize() override;

protected:
    std::shared_ptr<ColumnLayoutRelationData> relation_;

    ColumnLayoutRelationData const& GetRelation() const noexcept {
        // GetRelation should be called after input file is parsed i.e. after algorithm execution
        assert(relation_ != nullptr);
        return *relation_;
    }

public:
    explicit PliBasedFDAlgorithm(std::filesystem::path const& path, char separator = ',',
                                 bool has_header = true, bool const is_null_equal_null = true,
                                 std::vector<std::string_view> phase_names = {kDefaultPhaseName})
        : FDAlgorithm(path, separator, has_header, is_null_equal_null, std::move(phase_names)) {}

    explicit PliBasedFDAlgorithm(std::shared_ptr<ColumnLayoutRelationData> relation,
                                 std::vector<std::string_view> phase_names = {kDefaultPhaseName})
        : FDAlgorithm(std::move(phase_names)), relation_(std::move(relation)) {}

    std::vector<Column const*> GetKeys() const override;
};
