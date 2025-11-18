// Defines a header file containing function signatures for functions in src/

// Protect signatures using an inclusion guard.
#ifndef network_class_H
#define network_class_H
#define DARMA_USE_CURRENT
#include <RcppArmadillo.h>
#include <set>
#include <unordered_map>
#include "iglm/helper_functions.h"

class Network {
public:
  // Members
  bool directed;
  unsigned int number_edges; 
  // Note that adj_list refers to the outgoing and adj_list_in to the ingoing ties
  // for undirected networks adj_list_in is not initialized and never used 
  std::unordered_map< int, std::unordered_set<int>> adj_list;
  std::unordered_map< int, std::unordered_set<int>> adj_list_in;
  // Constructors
  Network (int n_actors_, bool directed_) {
    n_actors = n_actors_;
    directed = directed_;
    for (int i = 1; i <= n_actors; i++){
      adj_list[i] = std::unordered_set<int>();
      adj_list_in[i] = std::unordered_set<int>();
    }
    number_edges = 0;
  }
  
  Network (int n_actors_, bool directed_, arma::mat mat) {
    n_actors = n_actors_;
    directed = directed_;
    mat_to_map(mat,n_actors, directed, adj_list,adj_list_in);
    number_edges = count_edges(); 
  }
  // This is a member function to change the state of a particular pair in the network
  // (either from 1->0 or 0->1)
  void change_edge(int from, int to) {
    if(directed){
      if(adj_list.at(from).count(to)){
        adj_list.at(from).erase(to);
        adj_list_in.at(to).erase(from);
        number_edges--;
      } else {
        adj_list.at(from).insert(to);
        adj_list_in.at(to).insert(from);
        number_edges ++;
      }
    } else {
      if(adj_list.at(from).count(to)){
        adj_list.at(from).erase(to);
        adj_list.at(to).erase(from);
        number_edges--;
      } else {
        adj_list.at(from).insert(to);
        adj_list.at(to).insert(from);
        number_edges ++;
      }
    }
    
  }
  
  size_t count_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const {
    if(type == "OTP"){
      return(count_intersection(adj_list.at(from), adj_list_in.at(to))); 
    } else  if(type == "ISP"){
      return(count_intersection(adj_list_in.at(from), adj_list_in.at(to))); 
    }else  if(type == "OSP"){
      return(count_intersection(adj_list.at(from), adj_list.at(to))); 
    } 
    else  if(type == "ITP"){
      return(count_intersection(adj_list_in.at(from), adj_list.at(to))); 
    } 
    size_t res = 0; 
    return(res);
  } 
  
  std::unordered_set<int> get_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const {
    if(type == "OTP"){
      return(get_intersection(adj_list.at(from), adj_list_in.at(to))); 
    } else  if(type == "ISP"){
      return(get_intersection(adj_list_in.at(from), adj_list_in.at(to))); 
    }else  if(type == "OSP"){
      return(get_intersection(adj_list.at(from), adj_list.at(to))); 
    }
    else  if(type == "ITP"){
      return(get_intersection(adj_list_in.at(from), adj_list.at(to))); 
    }
    std::unordered_set<int> res; 
    return(res);
  }
  double count_edges() const{
    double count = 0.0;
    for(int i=1; i<=n_actors; i++){
      count += adj_list.at(i).size();
    }
    return(count);
  }
  
  double get_val(int from, int to)const {
    if(adj_list.at(from).count(to)){
      return(1.0);
    } else {
      return(0.0);
    }
  }
  
  // This is a member function to add a particular matrix
  void add_edge(int from, int to) {
    if(directed){
      adj_list.at(from).insert(to);
      adj_list_in.at(to).insert(from);
    } else{
      adj_list.at(from).insert(to);
      adj_list.at(to).insert(from);
    }
    number_edges ++;
  }
  // This is a member function to delete a particular matrix
  void delete_edge(int from, int to) {
    if(directed){
      adj_list.at(from).erase(to);
      adj_list_in.at(to).erase(from);
    } else{
      adj_list.at(from).erase(to);
      adj_list.at(to).erase(from);
    }
    number_edges --;
  }
  // This is a member function to initialize the matrix
  // void initialize() {
  //   for (int i = 1; i <= n_actors; i++){
  //     adj_list[i] = std::unordered_set<int>();
  //   }
  // }
  // This is a member function to add edges to a network through an adjacency matrix
  void add_edges_from_mat(arma::mat mat) {
    mat_to_map(mat,n_actors, directed, adj_list,adj_list_in);
  }
  
private:
  // Private members
  int n_actors;
};
#endif

