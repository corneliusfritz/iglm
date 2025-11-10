#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_NEWARP
#define ARMA_DONT_USE_ARPACK

#ifndef attribute_H
#define attribute_H

#include "helper_functions.h"
#include <string>
// #include <random>
// #include <set>
// #include <iostream>
// #include <unordered_map>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



class Attribute {
public:
  Attribute(int, std::string, double); // Constructor
  Attribute(int, arma::vec vector_tmp, std::string, double); // Alternative Constructor
  // Attributes (public)
  arma::vec attribute;
  double scale; 
  std::string type;
  // Member functions
  // void initialize();
  void print();
  double get_val(int from) const;
  double get_scale() const;
  int get_sum() const;
  bool check() const;
  void set_attr_0(int from);
  void set_attr_1(int from);
  void set_attr_value(int, double);
  void set_attr( arma::vec vector_tmp);
  void set_attr_from_armavec(arma::vec vwector_tmp);
  void change_attr(int from);

private:
  // Attributes (private)
  int n_actors;
};

#endif

