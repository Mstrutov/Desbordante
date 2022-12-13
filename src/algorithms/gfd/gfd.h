#pragma once
#include <vector>

#include "pattern.h"

typedef std::pair<int, std::string> Token;
typedef std::pair<Token, Token> Literal;

class GFD {

private:
  Pattern pattern;
  std::vector<Literal> premises;
  std::vector<Literal> conclusion;
  void print(std::vector<Literal>) const;

public:
  GFD() {};
  GFD(Pattern &pattern_, std::vector<Literal> &premises_,
      std::vector<Literal> &conclusion_)
      : pattern(pattern_), premises(premises_), conclusion(conclusion_) {}

  Pattern getPattern() const { return this->pattern; }
  std::vector<Literal> getPremises() const { return this->premises; }
  std::vector<Literal> getConclusion() const { return this->conclusion; }

  void print() const;
};
