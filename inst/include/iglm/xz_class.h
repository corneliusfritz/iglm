#pragma once

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include "attribute_class.h"
#include "network_class.h"
#define DARMA_USE_CURRENT

class IGLM_API XZ_class {
public:
  // Member
  int n_actor;
  Network z_network;
  std::vector<std::vector<int>> overlap;
  Attribute x_attribute;
  // Other members
  std::vector<std::vector<int>> neighborhood;
  std::vector<std::vector<int>> adj_list_nb;
  std::vector<std::vector<int>> adj_list_in_nb;
  std::vector<char> overlap_bool_mat;
  std::vector<char> neighborhood_bool_mat;
  arma::mat overlap_mat;
  std::vector<int> all_actors;

  int N_total_overlap;
  int N_1_overlap;
  inline size_t get_mat_idx(int from, int to) const {
    return (from - 1) * n_actor + (to - 1);
  }

  // Constructors
  XZ_class(int n_actor_, bool directed_, std::string type_, double scale_);
  XZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_);
  XZ_class(int n_actor_, bool directed_, std::vector<std::vector<int>> neighborhood_,
           std::vector<std::vector<int>> overlap_,
           arma::mat overlap_mat_, std::string type_, double scale_);
  XZ_class(int n_actor_, bool directed_, arma::mat z_network_, arma::vec x_attribute_,
           arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_);
  
  // Member functions
  void set_network_from_mat(int n_actor_, bool directed_, arma::mat mat);
  void initialize_overlap_counts();
  void add_edge(int from, int to);
  void delete_edge(int from, int to);

  inline bool get_val_overlap(int from, int to )const {
    return overlap_bool_mat[get_mat_idx(from, to)] || overlap_bool_mat[get_mat_idx(to, from)];
  }

  double count_edges() const;
  double count_nb_edges() const;
  
  std::vector<int> get_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const;
  size_t count_common_partners(unsigned int from, unsigned int to, std::string type = "OSP") const;
  
  
  std::vector<int> get_common_partners_nb(unsigned int from,unsigned int to, std::string type = "OSP")const;
  size_t count_common_partners_nb(unsigned int from, unsigned int to, std::string type = "OSP") const;
  
  inline bool get_val_neighborhood(int from, int to ) const {
    return neighborhood_bool_mat[get_mat_idx(from, to)];
  }
  
  bool check_if_full_neighborhood() const;
  void print();
  void copy_from(XZ_class obj);
  void set_neighborhood_from_mat(arma::mat mat);
  void neighborhood_initialize();
  void assign_neighborhood(const std::unordered_map< int, std::unordered_set<int>>& new_neighborhood);
  void change_neighborhood(int actor, std::unordered_set<int> new_neighborhood);
};
