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
  // Constructors
  XYZ_class(int n_actor_, bool directed_, std::string type_x_,std::string type_y_, 
                       double scale_x_, double scale_y_): XZ_class(n_actor_,directed_, type_x_, scale_x_), y_attribute(n_actor_, type_y_, scale_y_){
    
  } 
  XYZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_x_,std::string type_y_, double scale_x_, double scale_y_): 
    XZ_class(n_actor_,directed_, neighborhood_, overlap_, type_x_, scale_x_), y_attribute(n_actor_, type_y_, scale_y_){
    n_actor = n_actor_;
  }
  
  XYZ_class(int n_actor_, bool directed_,std::unordered_map< int, std::unordered_set<int>> neighborhood_,
                       std::unordered_map< int, std::unordered_set<int>> overlap_, 
                       arma::mat overlap_mat_, 
                       std::string type_x_,std::string type_y_,
                       double scale_x_, double scale_y_): XZ_class(n_actor_,directed_, neighborhood_, overlap_, overlap_mat_, type_x_, scale_x_),  y_attribute(n_actor_, type_y_, scale_y_){
    n_actor = n_actor_;
  }
  
  
  XYZ_class(int n_actor_, bool directed_,  arma::vec x_attribute_, arma::vec y_attribute_, arma::mat z_network_,arma::mat neighborhood_, 
                       arma::mat overlap_, std::string type_x_,std::string type_y_, double scale_x_, double scale_y_):
    XZ_class(n_actor_,directed_, z_network_,x_attribute_, neighborhood_,overlap_, type_x_, scale_x_),  y_attribute(n_actor_,y_attribute_, type_y_, scale_y_){
  }
  
  void print() {
    Rcout << "X Attribute" << std::endl;
    x_attribute.print();
    Rcout << "Y Attribute" << std::endl;
    y_attribute.print();
    Rcout << "Z Network" << std::endl;
    Rcout << map_to_mat(z_network.adj_list, n_actor) << std::endl;
    Rcout << "Neighborhood Matrix" << std::endl;
    Rcout << map_to_mat(neighborhood, n_actor) << std::endl;
  }
  
  void set_info(arma::vec x_attribute_,arma::vec y_attribute_,
                           std::unordered_map< int, std::unordered_set<int>> z_network_) {
    x_attribute.attribute = x_attribute_;
    y_attribute.attribute = y_attribute_;
    z_network.adj_list = z_network_;
  }
  
  void set_info_arma(arma::vec x_attribute_, arma::vec y_attribute_, arma::mat z_network_) {
    x_attribute.attribute = x_attribute_;
    y_attribute.attribute = y_attribute_;
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
  
  
  void copy_from(XYZ_class obj) {
    x_attribute = obj.x_attribute;
    y_attribute = obj.y_attribute;
    z_network = obj.z_network;
    neighborhood = obj.neighborhood;
    n_actor = obj.n_actor;
  }
  
};
