#pragma once

#include <memory>
#include <bitset>
#include <vector>
#include <list>

// For printing Dependencies
#include <string>
#include <fstream>
#include <iostream>

#include "RelationalSchema.h"
#include "FD.h"

class FDTreeElement {
public:
    // The maximum number of columns in the dataset. Using in std::bitset template.
    static constexpr int kMaxAttrNum = 256;
    explicit FDTreeElement(size_t max_attribute_number);

    FDTreeElement(const FDTreeElement&) = delete;
    FDTreeElement& operator=(const FDTreeElement&) = delete;

    void AddMostGeneralDependencies();

    // Using in cover-trees as post filtration of functional dependencies with redundant left-hand side.
    void FilterSpecializations();

    [[nodiscard]] bool CheckFd(size_t index) const;

    [[nodiscard]] FDTreeElement* GetChild(size_t index) const;

    void AddFunctionalDependency(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num);

    // Searching for generalization of functional dependency in cover-trees.
    bool GetGeneralizationAndDelete(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num,
                                    size_t current_attr, std::bitset<kMaxAttrNum>& spec_lhs);

    [[nodiscard]] bool ContainsGeneralization(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num,
                                              size_t current_attr) const;

    // Printing found dependencies in output file.
    void PrintDep(const std::string& file, std::vector<std::string>& column_names) const;

    void FillFdCollection(const RelationalSchema& scheme, std::list<FD>& fd_collection) const;

private:
    std::vector<std::unique_ptr<FDTreeElement>> children_;
    std::bitset<kMaxAttrNum> rhs_attributes_;
    size_t max_attribute_number_;
    std::bitset<kMaxAttrNum> is_fd_;

    void AddRhsAttribute(size_t index);

    [[nodiscard]] const std::bitset<kMaxAttrNum>& GetRhsAttributes() const;

    void MarkAsLast(size_t index);

    // Checking whether node is a leaf or not.
    [[nodiscard]] bool IsFinalNode(size_t attr_num) const;

    // Searching for specialization of functional dependency in cover-trees.
    bool GetSpecialization(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num,
                           size_t current_attr, std::bitset<kMaxAttrNum>& spec_lhs_out) const;

    void FilterSpecializationsHelper(FDTreeElement& filtered_tree,
                                     std::bitset<kMaxAttrNum>& active_path);

    // Helper function for PrintDep.
    void PrintDependencies(std::bitset<kMaxAttrNum>& active_path, std::ofstream& file,
                           std::vector<std::string>& column_names) const;

    void TransformTreeFdCollection(std::bitset<kMaxAttrNum>& active_path,
                                   std::list<FD>& fd_collection,
                                   const RelationalSchema& scheme) const;
};
