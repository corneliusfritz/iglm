#pragma once

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include "xz_class.h"
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
  
  XYZ_class(int n_actor_, bool directed_, std::vector<std::vector<int>> neighborhood_,
                       std::vector<std::vector<int>> overlap_, 
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
    Rcout << "XYZ Class: Implemented as dense flat structures" << std::endl;
  }
  
  void set_info(arma::vec x_attribute_,arma::vec y_attribute_,
                           std::vector<std::vector<int>> z_network_) {
    // Deprecated for vector refactor
  }
  
  void set_info_arma(arma::vec x_attribute_, arma::vec y_attribute_, arma::mat z_network_) {
    x_attribute.attribute = x_attribute_;
    y_attribute.attribute = y_attribute_;
    mat_to_map_vec(z_network_, n_actor, z_network.directed, z_network.adj_list, z_network.adj_list_in, z_network.adj_mat);
    for (int i = 1; i <= n_actor; i++){
      adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
      if(z_network.directed){
        adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
      } 
    }
    initialize_overlap_counts();
  }
  
  void copy_from(XYZ_class obj) {
    XZ_class::copy_from(obj);
    y_attribute = obj.y_attribute;
  }
  
};

