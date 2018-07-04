#pragma once

#include<stdexcept>
#include<iostream>
#include<unordered_set>
#include"bdmc_pseudo.hpp"

///This class takes care of hashing pseudomomenta in the bold code
///pseudo-momenta are between zero and one.
class hash{
public:
  ///add a boson momentum line
  bool add_momentum(pseudo_mom mom);

  ///check for existence of a boson momentum line
  bool has_momentum(pseudo_mom mom) const;

private:
  // Note that we base all the hashing and comparing on floats. This will take
  // care of round-off errors
  struct hasher {
      std::size_t operator() (pseudo_mom x) const {
          return std::hash<float>()(x());
      }
  };
  struct comparer {
      bool operator() (pseudo_mom x, pseudo_mom y) const {
          return float(x()) == float(y());
      }
  };

  ///hash for the propagators.
  std::unordered_set<pseudo_mom, hasher, comparer> momenta_;
};
