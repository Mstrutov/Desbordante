#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>

#include "algorithms/gfd/gfd_validation.h"
#include "algorithms/gfd/pattern.h"
#include "algorithms/gfd/parser.h"

using edge = std::pair<int, int>;
using attr = std::map<std::string, std::string>;
using field = std::pair<int, std::string>;
using namespace GFDtools;
auto current_path = std::filesystem::current_path() / "input_data" / "graph_data";

TEST(GFDValidationTest, TestTrivially) {
  auto graph_path = current_path / "quadrangle.txt";
  std::ifstream f;
  f.open(graph_path);
  Graph quadrilaterals = Parser::parse_graph(f);
  f.close();
  auto gfd_path = current_path / "family_tree_gfd3.txt";
  f.open(gfd_path);
  GFD gfd = Parser::parse_gfd(f);
  f.close();
  
  int expected_size = 1;

  auto algorithm = algos::GFDValidation(quadrilaterals, std::vector<GFD>{gfd});
  algorithm.Execute();
  std::vector<GFD> GFDList = algorithm.GFDList();

  EXPECT_EQ(expected_size, GFDList.size());
  EXPECT_EQ(gfd.getPattern(), GFDList.at(0).getPattern());
}

TEST(GFDValidationTest, TestExistingMatches0) {
  auto graph_path = current_path / "directors.txt";
  std::ifstream f;
  f.open(graph_path);
  Graph directors = Parser::parse_graph(f);
  f.close();
  auto gfd_path = current_path / "directors_gfd.txt";
  f.open(gfd_path);
  GFD connection_director_film = Parser::parse_gfd(f);
  f.close();
  
  int expected_size = 0;

  auto algorithm = algos::GFDValidation(
      directors, std::vector<GFD>{connection_director_film});
  algorithm.Execute();
  std::vector<GFD> GFDList = algorithm.GFDList();

  EXPECT_EQ(expected_size, GFDList.size());
}

TEST(GFDValidationTest, TestExistingMatches1) {
  std::ifstream f;
  auto graph_path = current_path / "family_tree.txt";
  f.open(graph_path);
  Graph family = Parser::parse_graph(f);
  f.close();
  auto gfd_path1 = current_path / "family_tree_gfd1.txt";
  f.open(gfd_path1);
  GFD eyes_simple = Parser::parse_gfd(f);
  f.close();
  auto gfd_path2 = current_path / "family_tree_gfd2.txt";
  f.open(gfd_path2);
  GFD eyes_complex = Parser::parse_gfd(f);
  f.close();
  auto gfd_path3 = current_path / "family_tree_gfd3.txt";
  f.open(gfd_path3);
  GFD body_types = Parser::parse_gfd(f);
  f.close();
    
  int expected_size = 2;

  auto algorithm = algos::GFDValidation(
      family, std::vector<GFD>{eyes_simple, eyes_complex, body_types});
  algorithm.Execute();
  std::vector<GFD> GFDList = algorithm.GFDList();

  EXPECT_EQ(expected_size, GFDList.size());

  for (const GFD &current : GFDList) {
    EXPECT_TRUE(current.getPattern() != eyes_simple.getPattern());
  }
  EXPECT_TRUE(GFDList.begin()->getPattern() !=
               (++GFDList.begin())->getPattern());

  for (const GFD &current : GFDList) {
    EXPECT_TRUE(current.getPattern() == eyes_complex.getPattern() ||
                current.getPattern() == body_types.getPattern());
  }
}
