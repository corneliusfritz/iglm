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
  int n_actor;
  Network z_network;
  std::unordered_map< int, std::unordered_set<int>> overlap;
  Attribute x_attribute;
  // Other members
  std::unordered_map< int, std::unordered_set<int>> neighborhood;
  std::unordered_map< int, std::unordered_set<int>> adj_list_nb;
  std::unordered_map< int, std::unordered_set<int>> adj_list_in_nb;
  arma::mat overlap_mat;
  std::unordered_set<int> all_actors;
  // Constructors
  // Constructor 1 (Corrected)
  XZ_class(int n_actor_, bool directed_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),              
    overlap(),                                  
    x_attribute(n_actor_, type_, scale_)         
  {
    for (int i = 1; i <= n_actor; i++){ 
      // Rcout << i  << std::endl;
      neighborhood[i] = std::unordered_set<int>();
      adj_list_nb[i] = std::unordered_set<int>();
      adj_list_in_nb[i] = std::unordered_set<int>();
      overlap[i] = std::unordered_set<int>();
    } 
    overlap_mat = arma::zeros<arma::mat>(2, 0);
    
    for (int i = 1; i <= n_actor; ++i) { 
      all_actors.insert(i);
    } 
  }
  
  // Constructor 2
  XZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),             
    overlap(),                                 
    x_attribute(n_actor_, type_, scale_)        
  {
    mat_to_map_neighborhood(neighborhood_, overlap_, n_actor, directed_,
                            neighborhood,
                            overlap);
    overlap_mat = overlap_;
    
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = std::unordered_set<int>();
      adj_list_in_nb[i] = std::unordered_set<int>();
      all_actors.insert(i);
    } 
  }
  
  // Constructor 3 
  XZ_class(int n_actor_, bool directed_, std::unordered_map< int, std::unordered_set<int>> neighborhood_,
           std::unordered_map< int, std::unordered_set<int>> overlap_,
           arma::mat overlap_mat_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),                
    overlap(overlap_),                            
    x_attribute(n_actor_, type_, scale_),         
    neighborhood(neighborhood_),                  
    overlap_mat(overlap_mat_)                     
  {
    // Rcout << "Hello"  << std::endl;
    // This loop is correct and uses the initialized members.
    for (int i = 1; i <= n_actor; i++){ // Use the member n_actor
      // Rcout << i  << std::endl;
      adj_list_nb[i] = get_intersection(z_network.adj_list.at(i),overlap.at(i));
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection(z_network.adj_list_in.at(i),overlap.at(i));
      }
      all_actors.insert(i);
    }
  }
  
  // Constructor 4 
  XZ_class(int n_actor_, bool directed_, arma::mat z_network_, arma::vec x_attribute_,
           arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
    n_actor(n_actor_),                               
    z_network(n_actor_, directed_, z_network_),      
    overlap(),                                       
    x_attribute(n_actor_, x_attribute_, type_, scale_)
  {
    // Rcout << "Hello 4"  << std::endl;
    mat_to_map_neighborhood(neighborhood_, overlap_, n_actor, directed_,
                            neighborhood,
                            overlap);
    overlap_mat = overlap_;
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = get_intersection(z_network.adj_list.at(i),overlap.at(i));
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection(z_network.adj_list_in.at(i),overlap.at(i));
      }
      all_actors.insert(i);
    }
  } 
  // Member functions
  void set_network_from_mat(int n_actor_, bool directed_, arma::mat mat){
    z_network.set_network_from_mat(n_actor_, directed_, mat); 
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = get_intersection(z_network.adj_list.at(i),overlap.at(i));
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection(z_network.adj_list_in.at(i),overlap.at(i));
      } 
    }
    
  }
  
  void add_edge(int from, int to) {
    if(z_network.directed){
      z_network.adj_list.at(from).insert(to);
      z_network.adj_list_in.at(to).insert(from);
      if(overlap.at(from).count(to)){
        adj_list_nb.at(from).insert(to);
        adj_list_in_nb.at(to).insert(from);
      }
    } else{
      z_network.adj_list.at(from).insert(to);
      z_network.adj_list.at(to).insert(from);
      if(overlap.at(from).count(to)){
        adj_list_nb.at(from).insert(to);
        adj_list_nb.at(to).insert(from);
      }
    } 
    
  } 
  // This is a member function to delete a particular matrix
  void delete_edge(int from, int to) {
    if(z_network.directed){
      z_network.adj_list.at(from).erase(to);
      z_network.adj_list_in.at(to).erase(from);
      if(overlap.at(from).count(to)){
        adj_list_nb.at(from).erase(to);
        adj_list_in_nb.at(to).erase(from);
      }
    } else{ 
      z_network.adj_list.at(from).erase(to);
      z_network.adj_list.at(to).erase(from);
      if(overlap.at(from).count(to)){
        adj_list_nb.at(from).erase(to);
        adj_list_nb.at(to).erase(from);
      }
    } 
  }
  bool get_val_overlap(int from, int to )const {
    if(overlap.at(from).count(to)){
      return(true);
    } else if(overlap.at(to).count(from)) {
      return(true);
    } else {
      return(false);
    }
  }
  double count_edges() const{
    double count = 0.0;
    for(int i=1; i<=n_actor; i++){
      count += z_network.adj_list.at(i).size();
    }
    return(count);
  }
  double count_nb_edges() const{
    double count = 0.0;
    for(int i=1; i<=n_actor; i++){
      count += adj_list_nb.at(i).size();
    }
    return(count);
  }
  
  std::unordered_set<int> get_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const {
    return(z_network.get_common_partners(from,to,type));
  }
  size_t count_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const {
    return(z_network.count_common_partners(from,to,type));
  }
  std::unordered_set<int> get_common_partners_nb(unsigned int from,unsigned int to, std::string type = "OSP")const {
    if(type == "OTP"){
      return(get_intersection(adj_list_nb.at(from),
                              adj_list_in_nb.at(to))); 
    } else  if(type == "ISP"){ 
      return(get_intersection(adj_list_in_nb.at(from),
                              adj_list_in_nb.at(to))); 
    }else  if(type == "OSP"){ 
      return(get_intersection(adj_list_nb.at(from),
                              adj_list_nb.at(to))); 
    } 
    else  if(type == "ITP"){
      return(get_intersection(adj_list_in_nb.at(from), 
                              adj_list_nb.at(to))); 
    } 
    std::unordered_set<int> res; 
    return(res);
  } 
  
  size_t count_common_partners_nb(unsigned int from,unsigned int to, std::string type = "OSP")const {
    if(type == "OTP"){
      return(count_intersection(adj_list_nb.at(from),
                              adj_list_in_nb.at(to))); 
    } else  if(type == "ISP"){ 
      return(count_intersection(adj_list_in_nb.at(from),
                              adj_list_in_nb.at(to))); 
    }else  if(type == "OSP"){  
      return(count_intersection(adj_list_nb.at(from),
                              adj_list_nb.at(to))); 
    }  
    else  if(type == "ITP"){
      return(count_intersection(adj_list_in_nb.at(from), 
                              adj_list_nb.at(to))); 
    } else{
      return(0);
    }
  } 
  
  bool get_val_neighborhood(int from, int to ) const{
    if(neighborhood.at(from).count(to)){
      return(true);
    } else {
      return(false);
    }
  }
  
  // This function checks if the neighborhood is full or not
  bool check_if_full_neighborhood() {
    int sum_tmp= 0;
    for(int i: seq(1,n_actor)){
      if((int) neighborhood.at(i).size() == n_actor){
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
    Rcout << "Network" << std::endl;
    Rcout << map_to_mat(z_network.adj_list, n_actor) << std::endl;
    Rcout << "Attribute" << std::endl;
    x_attribute.print();
    Rcout << "Neighborhood Matrix" << std::endl;
    Rcout << map_to_mat(neighborhood, n_actor) << std::endl;
    // print_set(neighborhood,n_ac);
  }
  
  void copy_from(XZ_class obj) {
    z_network = obj.z_network;
    x_attribute = obj.x_attribute;
    neighborhood = obj.neighborhood;
    n_actor = obj.n_actor;
  }
  
  
  void set_neighborhood_from_mat(arma::mat mat) {
    arma::rowvec tmp_row;
    for (int i = 1; i <= n_actor; i++){
      tmp_row = mat.row(i-1);
      arma::uvec ids = find(tmp_row == 1) + 1; // Find indices
      // Rcout << ids << std::endl;
      neighborhood.at(i)= std::unordered_set<int>(ids.begin(),ids.end());
    }
  }
  
  void neighborhood_initialize() {
    for (int i = 1; i <= n_actor; i++){
      neighborhood[i] = std::unordered_set<int>();
    }
  }
  
  void assign_neighborhood(std::unordered_map< int, std::unordered_set<int>> new_neighborhood) {
    neighborhood = new_neighborhood;
    
  }
  void set_info(arma::vec x_attribute_, std::unordered_map< int, std::unordered_set<int>> z_network_) {
    x_attribute.attribute = x_attribute_;
    z_network.adj_list = z_network_;
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = get_intersection(z_network.adj_list.at(i),overlap.at(i));
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection(z_network.adj_list_in.at(i),overlap.at(i));
      } 
    }
    
  }
  
  void set_info_arma(arma::vec x_attribute_, arma::mat z_network_) {
    x_attribute.attribute = x_attribute_;
    mat_to_map(z_network_,n_actor, z_network.directed, z_network.adj_list,z_network.adj_list_in);
    // Rcout << "B" << std::endl;
    // Rcout << n_actor << std::endl;
    for (int i = 1; i <= n_actor; i++){
      // Rcout << "Updating adj_list_nb for actor " << i << std::endl;
      adj_list_nb[i] = get_intersection(z_network.adj_list.at(i),overlap.at(i));
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection(z_network.adj_list_in.at(i),overlap.at(i));
      } 
    }
    
  }
  
  void change_neighborhood(int actor, std::unordered_set<int> new_neighborhood) {
    neighborhood.at(actor -1 ) = new_neighborhood;
  }
};
