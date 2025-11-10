#pragma once

// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <random>
#include <set>
#include <unordered_map>
#include "xz_class.h"
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_NEWARP
#define ARMA_DONT_USE_ARPACK
#define DARMA_USE_CURRENT

class XYZ_class: public XZ_class {
public:
  // Additional Member
  Attribute y_attribute;
  // arma::vec gradient;
  // Constructor
  XYZ_class (int, bool, std::string, std::string, double, double); 
  // Empty constructor without neighborhood
  XYZ_class (int, bool, arma::mat, arma::mat, std::string, std::string, 
             double, double); 
  // Empty network and attribute but with neighborhood
  XYZ_class (int, bool, std::unordered_map< int, std::unordered_set<int>>, 
             std::unordered_map< int, std::unordered_set<int>>, 
             arma::mat, std::string, std::string, 
             double, double); 
  // Empty network and attribute but with neighborhood (given as set of integers)
  XYZ_class (int, bool, arma::vec, arma::vec, arma::mat, arma::mat, arma::mat,
             std::string, std::string, double, double); // Matrix constructor
  // Additional Member functions (all are overloaded)
  void copy_from(XYZ_class);
  void print();
  void set_info(arma::vec, arma::vec, std::unordered_map< int, std::unordered_set<int>>);
  void set_info_arma(arma::vec,  arma::vec, arma::mat);
  // void initialize_from_zero(arma::mat);
  // void initialize(arma::mat, arma::vec,arma::vec, arma::mat);
};
