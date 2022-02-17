#pragma once

#include "ColumnLayoutRelationData.h"
#include "FDAlgorithm.h"

class PliBasedFDAlgorithm : public FDAlgorithm {
private:
    using algos::Primitive::input_generator_;

    void Initialize() override;

protected:
    std::shared_ptr<ColumnLayoutRelationData> relation_;

    ColumnLayoutRelationData const& GetRelation() const noexcept {
        // GetRelation should be called after input file is parsed i.e. after algorithm execution
        assert(relation_ != nullptr);
        return *relation_;
    }

public:
    explicit PliBasedFDAlgorithm(Config const& config, std::vector<std::string_view> phase_names)
        : FDAlgorithm(config, std::move(phase_names)) {}

    explicit PliBasedFDAlgorithm(std::shared_ptr<ColumnLayoutRelationData> relation,
                                 Config const& config, std::vector<std::string_view> phase_names)
        : FDAlgorithm(config, std::move(phase_names)), relation_(std::move(relation)) {}

    std::vector<Column const*> GetKeys() const override;
};
