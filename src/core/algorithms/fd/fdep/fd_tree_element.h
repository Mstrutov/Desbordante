#pragma once

#include <bitset>
#include <list>
#include <memory>
#include <vector>

// For printing Dependencies
#include <fstream>
#include <string>

#include "algorithms/fd/fd.h"
#include "model/table/relational_schema.h"

class FDTreeElement {
public:
    // The maximum number of columns in the dataset. Using in std::bitset template.
    static constexpr int kMaxAttrNum = 256;
    explicit FDTreeElement(size_t max_attribute_number);

    FDTreeElement(FDTreeElement const&) = delete;
    FDTreeElement& operator=(FDTreeElement const&) = delete;

    void AddMostGeneralDependencies();

    // Using in cover-trees as post filtration of functional dependencies with redundant left-hand
    // side.
    void FilterSpecializations();

    [[nodiscard]] bool CheckFd(size_t index) const;

    [[nodiscard]] FDTreeElement* GetChild(size_t index) const;

    void AddFunctionalDependency(std::bitset<kMaxAttrNum> const& lhs, size_t attr_num);

    // Searching for generalization of functional dependency in cover-trees.
    bool GetGeneralizationAndDelete(std::bitset<kMaxAttrNum> const& lhs, size_t attr_num,
                                    size_t current_attr, std::bitset<kMaxAttrNum>& spec_lhs);

    [[nodiscard]] bool ContainsGeneralization(std::bitset<kMaxAttrNum> const& lhs, size_t attr_num,
                                              size_t current_attr) const;

    // Printing found dependencies in output file.
    void PrintDep(std::string const& file, std::vector<std::string>& column_names) const;

    void FillFdCollection(RelationalSchema const& scheme, std::list<FD>& fd_collection,
                          unsigned int max_lhs_) const;

private:
    std::vector<std::unique_ptr<FDTreeElement>> children_;
    std::bitset<kMaxAttrNum> rhs_attributes_;
    size_t max_attribute_number_;
    std::bitset<kMaxAttrNum> is_fd_;

    void AddRhsAttribute(size_t index);

    [[nodiscard]] std::bitset<kMaxAttrNum> const& GetRhsAttributes() const;

    void MarkAsLast(size_t index);

    // Checking whether node is a leaf or not.
    [[nodiscard]] bool IsFinalNode(size_t attr_num) const;

    // Searching for specialization of functional dependency in cover-trees.
    bool GetSpecialization(std::bitset<kMaxAttrNum> const& lhs, size_t attr_num,
                           size_t current_attr, std::bitset<kMaxAttrNum>& spec_lhs_out) const;

    void FilterSpecializationsHelper(FDTreeElement& filtered_tree,
                                     std::bitset<kMaxAttrNum>& active_path);

    // Helper function for PrintDep.
    void PrintDependencies(std::bitset<kMaxAttrNum>& active_path, std::ofstream& file,
                           std::vector<std::string>& column_names) const;

    void TransformTreeFdCollection(std::bitset<kMaxAttrNum>& active_path,
                                   std::list<FD>& fd_collection, RelationalSchema const& scheme,
                                   unsigned int max_lhs_) const;
};
