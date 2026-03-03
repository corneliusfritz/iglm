// Defines a header file containing function signatures for functions in src/

// Protect signatures using an inclusion guard.
#ifndef network_class_H
#define network_class_H
#define DARMA_USE_CURRENT
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include "iglm/helper_functions.h"

class Network {
public:
  // Members
  bool directed;
  unsigned int number_edges; 
  // Note that adj_list refers to the outgoing and adj_list_in to the ingoing ties
  // for undirected networks adj_list_in is not initialized and never used 
  std::vector<std::vector<int>> adj_list;
  std::vector<std::vector<int>> adj_list_in;
  std::vector<char> adj_mat; // Flat representation for O(1) existence checks

  inline size_t get_mat_idx(int from, int to) const {
    return (from - 1) * n_actor + (to - 1);
  }

  // Constructors
  Network (int n_actor_, bool directed_) {
    n_actor = n_actor_;
    directed = directed_;
    adj_list.resize(n_actor + 1); // 1-indexed
    if (directed) {
      adj_list_in.resize(n_actor + 1);
    }
    adj_mat.assign(n_actor * n_actor, 0);
    number_edges = 0;
  }
  
  Network (int n_actor_, bool directed_, arma::mat mat) {
    n_actor = n_actor_;
    directed = directed_;
    adj_list.resize(n_actor + 1);
    if (directed) {
      adj_list_in.resize(n_actor + 1);
    }
    adj_mat.assign(n_actor * n_actor, 0);
    mat_to_map_vec(mat, n_actor, directed, adj_list, adj_list_in, adj_mat);
    number_edges = count_edges(); 
  }
  
  void set_network_from_mat(int n_actor_, bool directed_, arma::mat mat){
    n_actor = n_actor_;
    directed = directed_;
    adj_list.clear();
    adj_list.resize(n_actor + 1);
    if (directed) {
      adj_list_in.clear();
      adj_list_in.resize(n_actor + 1);
    }
    adj_mat.assign(n_actor * n_actor, 0);
    mat_to_map_vec(mat, n_actor, directed, adj_list, adj_list_in, adj_mat);
    number_edges = count_edges(); 
  }
  
  // This is a member function to change the state of a particular pair in the network
  // (either from 1->0 or 0->1)
  void change_edge(int from, int to) {
    if(directed){
      if(adj_mat[get_mat_idx(from, to)]) {
        delete_edge(from, to);
      } else {
        add_edge(from, to);
      }
    } else {
      if(adj_mat[get_mat_idx(from, to)]) {
        delete_edge(from, to);
      } else {
        add_edge(from, to);
      }
    }
  }
  
  size_t count_common_partners(unsigned int from, unsigned int to, std::string type = "OSP") const {
    if (type == "OTP") {
      const std::vector<int>& l1 = adj_list[from];
      const std::vector<int>& l2 = adj_list_in[to];
      if (l1.size() < l2.size()) {
        size_t count = 0;
        for (int k : l1) if (adj_mat[get_mat_idx(k, to)]) count++;
        return count;
      } else {
        size_t count = 0;
        for (int k : l2) if (adj_mat[get_mat_idx(from, k)]) count++;
        return count;
      }
    } else if (type == "ISP") {
      const std::vector<int>& l1 = adj_list_in[from];
      const std::vector<int>& l2 = adj_list_in[to];
      if (l1.size() < l2.size()) {
        size_t count = 0;
        for (int k : l1) if (adj_mat[get_mat_idx(k, to)]) count++;
        return count;
      } else {
        size_t count = 0;
        for (int k : l2) if (adj_mat[get_mat_idx(k, from)]) count++;
        return count;
      }
    } else if (type == "OSP") {
      const std::vector<int>& l1 = adj_list[from];
      const std::vector<int>& l2 = adj_list[to];
      if (l1.size() < l2.size()) {
        size_t count = 0;
        for (int k : l1) if (adj_mat[get_mat_idx(to, k)]) count++;
        return count;
      } else {
        size_t count = 0;
        for (int k : l2) if (adj_mat[get_mat_idx(from, k)]) count++;
        return count;
      }
    } else if (type == "ITP") {
      const std::vector<int>& l1 = adj_list_in[from];
      const std::vector<int>& l2 = adj_list[to];
      if (l1.size() < l2.size()) {
        size_t count = 0;
        for (int k : l1) if (adj_mat[get_mat_idx(to, k)]) count++;
        return count;
      } else {
        size_t count = 0;
        for (int k : l2) if (adj_mat[get_mat_idx(k, from)]) count++;
        return count;
      }
    }
    return 0;
  }
  
  std::vector<int> get_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const {
    if(type == "OTP"){
      return(get_intersection_vec(adj_list[from], adj_list_in[to])); 
    } else  if(type == "ISP"){
      return(get_intersection_vec(adj_list_in[from], adj_list_in[to])); 
    }else  if(type == "OSP"){
      return(get_intersection_vec(adj_list[from], adj_list[to])); 
    }
    else  if(type == "ITP"){
      return(get_intersection_vec(adj_list_in[from], adj_list[to])); 
    }
    std::vector<int> res; 
    return(res);
  }
  
  double count_edges() const{
    double count = 0.0;
    for(int i=1; i<=n_actor; i++){
      count += adj_list[i].size();
    }
    return(count);
  }
  
  inline double get_val(int from, int to) const {
    return adj_mat[get_mat_idx(from, to)] ? 1.0 : 0.0;
  }
  
  // This is a member function to add a particular matrix
  void add_edge(int from, int to) {
    if(directed){
      if (!adj_mat[get_mat_idx(from, to)]) {
        auto& o_list = adj_list[from];
        o_list.insert(std::lower_bound(o_list.begin(), o_list.end(), to), to);
        auto& i_list = adj_list_in[to];
        i_list.insert(std::lower_bound(i_list.begin(), i_list.end(), from), from);
        adj_mat[get_mat_idx(from, to)] = 1;
        number_edges ++;
      }
    } else{
      if (!adj_mat[get_mat_idx(from, to)]) {
        auto& list_f = adj_list[from];
        list_f.insert(std::lower_bound(list_f.begin(), list_f.end(), to), to);
        auto& list_t = adj_list[to];
        list_t.insert(std::lower_bound(list_t.begin(), list_t.end(), from), from);
        adj_mat[get_mat_idx(from, to)] = 1;
        adj_mat[get_mat_idx(to, from)] = 1;
        number_edges ++;
      }
    }
  }
  
  // This is a member function to delete a particular matrix
  void delete_edge(int from, int to) {
    if(directed){
      if (adj_mat[get_mat_idx(from, to)]) {
        adj_list[from].erase(std::remove(adj_list[from].begin(), adj_list[from].end(), to), adj_list[from].end());
        adj_list_in[to].erase(std::remove(adj_list_in[to].begin(), adj_list_in[to].end(), from), adj_list_in[to].end());
        adj_mat[get_mat_idx(from, to)] = 0;
        number_edges --;
      }
    } else{
      if (adj_mat[get_mat_idx(from, to)]) {
        adj_list[from].erase(std::remove(adj_list[from].begin(), adj_list[from].end(), to), adj_list[from].end());
        adj_list[to].erase(std::remove(adj_list[to].begin(), adj_list[to].end(), from), adj_list[to].end());
        adj_mat[get_mat_idx(from, to)] = 0;
        adj_mat[get_mat_idx(to, from)] = 0;
        number_edges --;
      }
    }
  }
  
  // This is a member function to add edges to a network through an adjacency matrix
  void add_edges_from_mat(arma::mat mat) {
    mat_to_map_vec(mat, n_actor, directed, adj_list, adj_list_in, adj_mat);
    number_edges = count_edges();
  }
  
private:
  // Private members
  int n_actor;
};
#endif
