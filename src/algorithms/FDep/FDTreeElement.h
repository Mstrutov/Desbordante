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

class FDTreeElement{
 public:
    // The maximum number of columns in the dataset. Using in std::bitset template.
    static constexpr int kMaxAttrNum = 256;
    explicit FDTreeElement(size_t maxAttributeNumber);

    FDTreeElement (const FDTreeElement&) = delete;
    FDTreeElement& operator=(const FDTreeElement&) = delete;

    void addMostGeneralDependencies();

    // Using in cover-trees as post filtration of functional dependencies with redundant left-hand side.
    void filterSpecializations();

    [[nodiscard]] bool checkFd(size_t index) const;

    [[nodiscard]] FDTreeElement* getChild(size_t index) const;

    void addFunctionalDependency(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num);

    // Searching for generalization of functional dependency in cover-trees.
    bool getGeneralizationAndDelete(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num,
                                    size_t currentAttr, std::bitset<kMaxAttrNum>& specLhs);

    [[nodiscard]] bool containsGeneralization(const std::bitset<kMaxAttrNum>& lhs, size_t attr_num,
                                              size_t currentAttr) const;

    // Printing found dependencies in output file.
    void printDep(const std::string& file, std::vector<std::string>& columnNames) const;

    void fillFdCollection(const RelationalSchema & scheme, std::list<FD> & fdCollection) const;
 private:
    std::vector<std::unique_ptr<FDTreeElement>> children_;
    std::bitset<kMaxAttrNum> rhsAttributes_;
    size_t maxAttributeNumber_;
    std::bitset<kMaxAttrNum> isFd_;

    void addRhsAttribute(size_t index);

    [[nodiscard]] const std::bitset<kMaxAttrNum>& getRhsAttributes() const;

    void markAsLast(size_t index);

    // Checking whether node is a leaf or not.
    [[nodiscard]] bool isFinalNode(size_t attr_num) const;

    // Searching for specialization of functional dependency in cover-trees.
    bool getSpecialization (const std::bitset<kMaxAttrNum>& lhs, size_t attr_num,
                            size_t currentAttr, std::bitset<kMaxAttrNum>& specLhsOut) const;

    void filterSpecializationsHelper(FDTreeElement& filteredTree, std::bitset<kMaxAttrNum>& activePath);

    // Helper function for printDep.
    void printDependencies(std::bitset<kMaxAttrNum>& activePath, std::ofstream& file,
                           std::vector<std::string>& columnNames) const;

    void transformTreeFdCollection(std::bitset<kMaxAttrNum>& activePath, std::list<FD> & fdCollection,
                                   const RelationalSchema & scheme) const;
};
