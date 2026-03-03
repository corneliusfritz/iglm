#pragma once

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include "attribute_class.h"
#include "network_class.h"
#define DARMA_USE_CURRENT

class XZ_class {
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

  inline size_t get_mat_idx(int from, int to) const {
    return (from - 1) * n_actor + (to - 1);
  }

  // Constructors
  // Constructor 1
  XZ_class(int n_actor_, bool directed_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),              
    x_attribute(n_actor_, type_, scale_)         
  {
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
    overlap_mat = arma::zeros<arma::mat>(2, 0);
    
    for (int i = 1; i <= n_actor; ++i) { 
      all_actors.push_back(i);
    } 
  }
  
  // Constructor 2
  XZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),             
    x_attribute(n_actor_, type_, scale_)        
  {
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);

    mat_to_map_vec(neighborhood_, n_actor, directed_, neighborhood, neighborhood, neighborhood_bool_mat);
    mat_to_map_vec(overlap_, n_actor, directed_, overlap, overlap, overlap_bool_mat);
    overlap_mat = overlap_;
    
    for (int i = 1; i <= n_actor; i++){
      all_actors.push_back(i);
    } 
  }
  
  // Constructor 3 
  XZ_class(int n_actor_, bool directed_, std::vector<std::vector<int>> neighborhood_,
           std::vector<std::vector<int>> overlap_,
           arma::mat overlap_mat_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),                
    x_attribute(n_actor_, type_, scale_),         
    overlap_mat(overlap_mat_)                     
  {
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);

    for (int i = 1; i <= n_actor; i++){ 
      for(int neighbor : neighborhood_[i]) {
        neighborhood[i].push_back(neighbor);
        neighborhood_bool_mat[get_mat_idx(i, neighbor)] = 1;
      }
      for(int over : overlap_[i]) {
        overlap[i].push_back(over);
        overlap_bool_mat[get_mat_idx(i, over)] = 1;
      }
      std::sort(neighborhood[i].begin(), neighborhood[i].end());
      std::sort(overlap[i].begin(), overlap[i].end());

      adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
      }
      all_actors.push_back(i);
    }
  }
  
  // Constructor 4 
  XZ_class(int n_actor_, bool directed_, arma::mat z_network_, arma::vec x_attribute_,
           arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
    n_actor(n_actor_),                               
    z_network(n_actor_, directed_, z_network_),      
    x_attribute(n_actor_, x_attribute_, type_, scale_)
  {
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);

    mat_to_map_vec(neighborhood_, n_actor, directed_, neighborhood, neighborhood, neighborhood_bool_mat);
    mat_to_map_vec(overlap_, n_actor, directed_, overlap, overlap, overlap_bool_mat);
    overlap_mat = overlap_;
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
      }
      all_actors.push_back(i);
    }
  } 
  
  // Member functions
  void set_network_from_mat(int n_actor_, bool directed_, arma::mat mat){
    z_network.set_network_from_mat(n_actor_, directed_, mat); 
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
      } 
    }
  }
  
  void add_edge(int from, int to) {
    if(z_network.directed){
      if(!z_network.get_val(from, to)){
        z_network.add_edge(from, to);
        if(overlap_bool_mat[get_mat_idx(from, to)]){
          auto& l_nb = adj_list_nb[from];
          l_nb.insert(std::lower_bound(l_nb.begin(), l_nb.end(), to), to);
          auto& li_nb = adj_list_in_nb[to];
          li_nb.insert(std::lower_bound(li_nb.begin(), li_nb.end(), from), from);
        }
      }
    } else{
      if(!z_network.get_val(from, to)){
        z_network.add_edge(from, to);
        if(overlap_bool_mat[get_mat_idx(from, to)]){
          auto& l_f = adj_list_nb[from];
          l_f.insert(std::lower_bound(l_f.begin(), l_f.end(), to), to);
          auto& l_t = adj_list_nb[to];
          l_t.insert(std::lower_bound(l_t.begin(), l_t.end(), from), from);
        }
      }
    } 
  } 

  void delete_edge(int from, int to) {
    if(z_network.directed){
      if(z_network.get_val(from, to)){
        z_network.delete_edge(from, to);
        if(overlap_bool_mat[get_mat_idx(from, to)]){
          adj_list_nb[from].erase(std::remove(adj_list_nb[from].begin(), adj_list_nb[from].end(), to), adj_list_nb[from].end());
          adj_list_in_nb[to].erase(std::remove(adj_list_in_nb[to].begin(), adj_list_in_nb[to].end(), from), adj_list_in_nb[to].end());
        }
      }
    } else{ 
      if(z_network.get_val(from, to)){
        z_network.delete_edge(from, to);
        if(overlap_bool_mat[get_mat_idx(from, to)]){
          adj_list_nb[from].erase(std::remove(adj_list_nb[from].begin(), adj_list_nb[from].end(), to), adj_list_nb[from].end());
          adj_list_nb[to].erase(std::remove(adj_list_nb[to].begin(), adj_list_nb[to].end(), from), adj_list_nb[to].end());
        }
      }
    } 
  }

  inline bool get_val_overlap(int from, int to )const {
    return overlap_bool_mat[get_mat_idx(from, to)] || overlap_bool_mat[get_mat_idx(to, from)];
  }

  double count_edges() const{
    return z_network.count_edges();
  }

  double count_nb_edges() const{
    double count = 0.0;
    for(int i=1; i<=n_actor; i++){
      count += adj_list_nb[i].size();
    }
    return(count);
  }
  
  std::vector<int> get_common_partners_nb(unsigned int from,unsigned int to, std::string type = "OSP")const {
    if(type == "OTP"){
      return(get_intersection_vec(adj_list_nb[from], adj_list_in_nb[to])); 
    } else  if(type == "ISP"){ 
      return(get_intersection_vec(adj_list_in_nb[from], adj_list_in_nb[to])); 
    }else  if(type == "OSP"){ 
      return(get_intersection_vec(adj_list_nb[from], adj_list_nb[to])); 
    } 
    else  if(type == "ITP"){
      return(get_intersection_vec(adj_list_in_nb[from], adj_list_nb[to])); 
    } 
    std::vector<int> res; 
    return(res);
  } 
  
  size_t count_common_partners_nb(unsigned int from, unsigned int to, std::string type = "OSP") const {
    const std::vector<int>* l1_ptr = nullptr;
    const std::vector<int>* l2_ptr = nullptr;
    if (type == "OTP") {
      l1_ptr = &adj_list_nb[from];
      l2_ptr = &adj_list_in_nb[to];
    } else if (type == "ISP") {
      l1_ptr = &adj_list_in_nb[from];
      l2_ptr = &adj_list_in_nb[to];
    } else if (type == "OSP") {
      l1_ptr = &adj_list_nb[from];
      l2_ptr = &adj_list_nb[to];
    } else if (type == "ITP") {
      l1_ptr = &adj_list_in_nb[from];
      l2_ptr = &adj_list_nb[to];
    } else {
      return 0;
    }

    const std::vector<int>& l1 = *l1_ptr;
    // l2_ptr is retained because it points to the potential matches, but we use the adjacency matrix for lookups instead.
    
    // We iterate over the smaller list and check indices in the larger list.
    // However, since we have the adjacency matrix and overlap matrix,
    // we can just check those regardless of l2's size.
    // Checking adj_mat/overlap_mat is O(1).
    size_t count = 0;
    if (type == "OTP") {
      for (int k : l1) if (z_network.adj_mat[get_mat_idx(k, to)] && overlap_bool_mat[get_mat_idx(k, to)]) count++;
    } else if (type == "ISP") {
      // ISP: k -> from and k -> to
      for (int k : l1) if (z_network.adj_mat[get_mat_idx(k, to)] && overlap_bool_mat[get_mat_idx(k, to)]) count++;
    } else if (type == "OSP") {
      // OSP: from -> k and to -> k
      for (int k : l1) if (z_network.adj_mat[get_mat_idx(to, k)] && overlap_bool_mat[get_mat_idx(to, k)]) count++;
    } else if (type == "ITP") {
      // ITP: k -> from and to -> k
      for (int k : l1) if (z_network.adj_mat[get_mat_idx(to, k)] && overlap_bool_mat[get_mat_idx(to, k)]) count++;
    }
    return count;
  }
  
  inline bool get_val_neighborhood(int from, int to ) const{
    return neighborhood_bool_mat[get_mat_idx(from, to)];
  }
  
  bool check_if_full_neighborhood() {
    int sum_tmp= 0;
    for(int i = 1; i <= n_actor; ++i) {
      if((int) neighborhood[i].size() == n_actor){
        sum_tmp ++;
      }
    }
    if(sum_tmp == n_actor){
      return(true);
    } else {
      return(false);
    }
  }

  void print() {
    Rcout << "Network: Implemented as dense flat structures" << std::endl;
  }
  
  void copy_from(XZ_class obj) {
    z_network = obj.z_network;
    x_attribute = obj.x_attribute;
    neighborhood = obj.neighborhood;
    overlap = obj.overlap;
    adj_list_nb = obj.adj_list_nb;
    adj_list_in_nb = obj.adj_list_in_nb;
    overlap_bool_mat = obj.overlap_bool_mat;
    neighborhood_bool_mat = obj.neighborhood_bool_mat;
    n_actor = obj.n_actor;
  }
  
  void set_neighborhood_from_mat(arma::mat mat) {
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
    for (int i = 1; i <= n_actor; i++){
      neighborhood[i].clear();
      for(int j = 1; j <= n_actor; j++) {
        if(mat(i-1, j-1) == 1) {
          neighborhood[i].push_back(j);
          neighborhood_bool_mat[get_mat_idx(i, j)] = 1;
        }
      }
    }
  }
  
  void neighborhood_initialize() {
    for (int i = 1; i <= n_actor; i++){
      neighborhood[i].clear();
    }
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
  }
  
  void assign_neighborhood(const std::unordered_map< int, std::unordered_set<int>>& new_neighborhood) {
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
    for (int i = 1; i <= n_actor; i++){
      neighborhood[i].clear();
      for(int neighbor : new_neighborhood.at(i)) {
        neighborhood[i].push_back(neighbor);
        neighborhood_bool_mat[get_mat_idx(i, neighbor)] = 1;
      }
      std::sort(neighborhood[i].begin(), neighborhood[i].end());
    }
  }
  
  void change_neighborhood(int actor, std::unordered_set<int> new_neighborhood) {
    for(int old_n : neighborhood[actor]) {
      neighborhood_bool_mat[get_mat_idx(actor, old_n)] = 0;
    }
    neighborhood[actor].clear();
    for(int new_n : new_neighborhood) {
      neighborhood[actor].push_back(new_n);
      neighborhood_bool_mat[get_mat_idx(actor, new_n)] = 1;
    }
    std::sort(neighborhood[actor].begin(), neighborhood[actor].end());
  }
};
