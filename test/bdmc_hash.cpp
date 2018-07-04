#include"bdmc_hash.hpp"

bool hash::add_momentum(pseudo_mom d){
  if (d() >= 1. || d() < 0.)
    throw std::range_error("pseudo momentum is out of range");

  return momenta_.insert(d).second;
}

bool hash::has_momentum(pseudo_mom d) const{
  return momenta_.find(d) != momenta_.end();
}
