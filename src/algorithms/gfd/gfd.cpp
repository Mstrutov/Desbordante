#include <iostream>

#include "gfd.h"

void GFD::print(std::vector<Literal> literals) const {
  for (const Literal &l : literals) {
    auto fst_token = l.first;
    auto snd_token = l.second;
    if (fst_token.first != -1) {
      std::cout << fst_token.first << ".";
    }
    std::cout << fst_token.second << "=";
    if (snd_token.first != -1) {
      std::cout << snd_token.first << ".";
    }
    std::cout << snd_token.second << "; ";
  }
}

void GFD::print() const {
  this->pattern.print();
  std::cout << std::endl << "Premises:   ";
  print(this->premises);
  std::cout << std::endl << "=>" << std::endl << "Conclusion: ";
  print(this->conclusion);
}
