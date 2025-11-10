#pragma once

// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <random>
#include <set>
#include <unordered_map>
#include "attribute_class.h"
#include "network_class.h"
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_NEWARP
#define ARMA_DONT_USE_ARPACK
#define DARMA_USE_CURRENT

class XZ_class {
public:
  // Member
  Attribute x_attribute;
  Network z_network;
  std::unordered_map< int, std::unordered_set<int>> neighborhood;
  std::unordered_map< int, std::unordered_set<int>> overlap;
  arma::mat overlap_mat;
  
  int n_actor;
  // Constructor
  XZ_class (int, bool, std::string, double); 
  // Empty constructor without neighborhood
  XZ_class (int, bool, arma::mat, arma::mat, std::string, double); 
  // Empty network and attribute but with neighborhood and overlap structure
  XZ_class (int, bool, std::unordered_map< int, std::unordered_set<int>>, 
            std::unordered_map< int, std::unordered_set<int>>, 
            arma::mat, std::string, double); 
  // Matrix constructor
  XZ_class (int, bool, arma::mat, arma::vec, arma::mat, arma::mat, std::string, double);
  
  // Member functions
  void assign_neighborhood(std::unordered_map< int, std::unordered_set<int>> new_neighborhood);
  void change_neighborhood(int actor, std::unordered_set<int> new_neighborhood);
  // void initialize(arma::mat, arma::vec, arma::mat);
  // void initialize_from_zero(arma::mat);
  bool check_if_full_neighborhood();
  void copy_from(XZ_class);
  void neighborhood_initialize();
  void print();
  void neighborhood_randomization(int seed, int K);
  void set_neighborhood_from_mat(arma::mat mat);
  bool get_val_neighborhood(int, int) const;
  bool get_val_overlap(int, int) const;
  void set_info(arma::vec, std::unordered_map< int, std::unordered_set<int>>);
  void set_info_arma(arma::vec, arma::mat);
};
