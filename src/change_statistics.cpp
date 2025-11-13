#include <RcppArmadillo.h>
#include "iglm/extension_api.hpp"
#include <random>
#include <set>
#include <unordered_map>
#include "iglm/xyz_class.h"


// The KEY MACRO: Define the signature wrapper (the lambda capture list is empty)
#define CHANGESTAT [](const XYZ_class &object,const int &actor_i,const int &actor_j, const arma::mat &data,const double &type,const std::string &mode,const bool &is_full_neighborhood) -> double


// Alias for a function pointer used to calculate a validation metric or score.
// This signature defines a **mapping** from a complex state space (captured by the seven
// reference arguments) to a single **double** metric. All parameters are passed by 
// constant reference, indicating the function may modify the state variables 
// (`actor_i`, `data`, `mode`, etc.) during its execution.
//  @param object A reference to the state container, `XYZ_class`.
//  @param actor_i, actor_j References to integer indices, likely representing **agents** or **nodes**.
//  @param data A reference to an `arma::mat` (Armadillo matrix), likely a **system matrix** (e.g., adjacency, feature set).
//  @param type A reference to a `double`, likely a **model parameter** or **hyperparameter**.
//  @param mode A reference to a `std::string`, specifying the **operational regime** (e.g., "fast", "full").
//  @param is_full_neighborhood A reference to a `bool` flag, controlling the scope or **boundary condition**.
//  @returns double The calculated validation score or fitness value.

auto xyz_stat_repetition = CHANGESTAT{
  if(mode == "z"){
    return(object.z_network.get_val(actor_j, actor_i));
  }  
  else{
    return(0);
  }  
}; 
EFFECT_REGISTER("mutual_global", ::xyz_stat_repetition, "mutual_global", 0);

auto xyz_stat_edges= CHANGESTAT
{
  
  if (mode == "z")
    return 1.0;
  else
    return 0.0;
}; 
// Register: name, function pointer, short name, double
EFFECT_REGISTER("edges_global", ::xyz_stat_edges, "edges_global", 0);

auto xyz_stat_repetition_nonb= CHANGESTAT{
  if(mode == "z"){
    return(object.z_network.get_val(actor_j, actor_i)*(1-object.get_val_overlap(actor_i,actor_j)));
  } 
  else{
    return(0);
  }  
};
EFFECT_REGISTER("mutual_alocal", ::xyz_stat_repetition_nonb, "mutual_alocal", 0);
auto xyz_stat_repetition_nb= CHANGESTAT{
  if(mode == "z"){
    return(object.z_network.get_val(actor_j, actor_i)*object.get_val_overlap(actor_i,actor_j));
  }
  else{
    return(0);
  } 
};
EFFECT_REGISTER("mutual_local", ::xyz_stat_repetition_nb, "mutual_local", 0);


auto xyz_stat_cov_z_out_nb= CHANGESTAT{
  if(mode == "z"){ 
    return(data.at(0, actor_i-1)*object.get_val_overlap(actor_i,actor_j));
  }else {
    return(0);
  } 
};
EFFECT_REGISTER("cov_z_out_local", ::xyz_stat_cov_z_out_nb, "cov_z_out_local", 0);


auto xyz_stat_cov_z_in_nb= CHANGESTAT{
  if(mode == "z"){ 
    return(data.at(0, actor_j-1)*object.get_val_overlap(actor_i,actor_j));
  }else {
    return(0);
  } 
};
EFFECT_REGISTER("cov_z_in_local", ::xyz_stat_cov_z_in_nb, "cov_z_in_local", 0);



auto xyz_stat_cov_z_out_nonb= CHANGESTAT{
  if(mode == "z"){ 
    return(data.at(0, actor_i-1)*(1-object.get_val_overlap(actor_i,actor_j)));
  }else {
    return(0);
  }  
};
EFFECT_REGISTER("cov_z_out_alocal", ::xyz_stat_cov_z_out_nonb, "cov_z_out_alocal", 0);

auto xyz_stat_cov_z_in_nonb= CHANGESTAT{
  if(mode == "z"){ 
    return(data.at(0, actor_j-1)*(1-object.get_val_overlap(actor_i,actor_j)));
  }else {
    return(0);
  } 
};
EFFECT_REGISTER("cov_z_in_alocal", ::xyz_stat_cov_z_in_nonb, "cov_z_in_alocal", 0);

auto xyz_stat_cov_z_out= CHANGESTAT{
  if(mode == "z"){ 
    return(data.at(0, actor_i-1));
  }else {
    return(0);
  }  
};
EFFECT_REGISTER("cov_z_out_global", ::xyz_stat_cov_z_out, "cov_z_out_global", 0);

auto xyz_stat_cov_z_in= CHANGESTAT{
  if(mode == "z"){ 
    return(data.at(0, actor_j-1));
  }else {
    return(0);
  } 
};
EFFECT_REGISTER("cov_z_in_global", ::xyz_stat_cov_z_in, "cov_z_in_global", 0);


auto xyz_stat_cov_z_nb= CHANGESTAT{
  if(mode == "z"){
    if(object.get_val_overlap(actor_i,actor_j)){
      return(data.at(actor_i-1, actor_j-1));  
    } else {
      return(0);
    }
    
  }
  else{
    return(0);
  }
};
EFFECT_REGISTER("cov_z_local", ::xyz_stat_cov_z_nb, "cov_z_local", 0);


auto xyz_stat_cov_z_nonb= CHANGESTAT{
  if(mode == "z"){
    if((1-object.get_val_overlap(actor_i,actor_j))){
      return(data.at(actor_i-1, actor_j-1));  
    } else {
      return(0);
    }
    
  }
  else{
    return(0);
  }
};
EFFECT_REGISTER("cov_z_alocal", ::xyz_stat_cov_z_nonb, "cov_z_alocal", 0);

auto xyz_stat_cov_z= CHANGESTAT{
  if(mode == "z"){
    return(data.at(actor_i-1, actor_j-1));
  }
  else{
    return(0);
  }
};
EFFECT_REGISTER("cov_z_global", ::xyz_stat_cov_z, "cov_z_global", 0);


auto xyz_stat_cov_x= CHANGESTAT{
  if(mode == "x"){
    return(data.at(0, actor_i-1));
  }
  else{
    return(0);
  }
};
EFFECT_REGISTER("cov_x", ::xyz_stat_cov_x, "cov_x", 0);

auto xyz_stat_cov_y= CHANGESTAT{
  if(mode == "y"){
    return(data.at(0, actor_i-1));
  }
  else{
    return(0);
  }
};
EFFECT_REGISTER("cov_y", ::xyz_stat_cov_y, "cov_y", 0);



auto xyz_stat_edges_nonb= CHANGESTAT{
  if(mode == "z"){
    if(is_full_neighborhood) {
      return(0);
    } else { 
      return(1-object.get_val_overlap(actor_i,actor_j));
      // bool same_group;
      // std::unordered_set<int> intersect_group;
      // std::set_intersection(std::begin(object.overlap.at(actor_i)),
      //                       std::end(object.overlap.at(actor_i)),
      //                       std::begin(object.overlap.at(actor_j)),
      //                       std::end(object.overlap.at(actor_j)),
      //                       std::inserter(intersect_group, std::begin(intersect_group)));
      // // The union of the two neighborhoods is saved in intersect_group
      // same_group = intersect_group.size()>0;
      // if(object.overlap.at(actor_i).count(actor_j) != object.get_val_overlap(actor_i,actor_j)){
      //   Rcout << "There is an issue between actors:" ;
      //   Rcout << actor_i ;
      //   Rcout << " and " ;
      //   Rcout << actor_j << std::endl;
      // }
      // if(same_group){
      //   return(0);
      // } else {
      //   return(1);
      // }
      // if(object.get_val_overlap(actor_i,actor_j)){
      //   return(0);
      // } else {
      //   return(1);
      // }
    }
    
  }
  else{
    return(0);
  } 
};
EFFECT_REGISTER("edges_alocal", ::xyz_stat_edges_nonb, "edges_alocal", 0);


auto xyz_stat_edges_nb= CHANGESTAT{
  if(mode == "z"){
    if(is_full_neighborhood) {
      return(1);
    } else {
      return(object.get_val_overlap(actor_i,actor_j));
      // bool same_group;
      // std::unordered_set<int> intersect_group;
      // std::set_intersection(std::begin(object.overlap.at(actor_i)),
      //                       std::end(object.overlap.at(actor_i)),
      //                       std::begin(object.overlap.at(actor_j)),
      //                       std::end(object.overlap.at(actor_j)),
      //                       std::inserter(intersect_group, std::begin(intersect_group)));
      // // The union of the two neighborhoods is saved in intersect_group
      // same_group = intersect_group.size()>0;
      // if(same_group){
      //   return(1);
      // } else {
      //   return(0);
      // }
    }
    
  }
  else{
    return(0);
  }
};
EFFECT_REGISTER("edges_local", ::xyz_stat_edges_nb, "edges_local", 0);


auto xyz_stat_attribute_xy_nb= CHANGESTAT{
  if(mode == "y"){
    int res = 0;
    for (auto itr = object.overlap.at(actor_i).begin(); itr != object.overlap.at(actor_i).end(); itr++) {
      res+= object.y_attribute.get_val(*itr);
    }
    return(res);
  } else if(mode == "x"){
    int res = 0;
    for (auto itr = object.overlap.at(actor_i).begin(); itr != object.overlap.at(actor_i).end(); itr++) {
      res+= object.x_attribute.get_val(*itr);
    }
    return(res);
  } else {
    return(0);
  }
};
EFFECT_REGISTER("attribute_xy_local", ::xyz_stat_attribute_xy_nb, "attribute_xy_local", 0);

auto xyz_stat_attribute_xy_nonb= CHANGESTAT{
  if(mode == "y"){
    
    std::unordered_set<int> difference_result =
      get_set_difference(object.all_actors, object.overlap.at(actor_i));;
    
    
    int res = 0;
    for (auto itr = difference_result.begin(); itr != difference_result.end(); itr++) {
      res+= object.y_attribute.get_val(*itr);
    }
    return(res);
  } else if(mode == "x"){ 
    std::unordered_set<int> difference_result = 
      get_set_difference(object.all_actors, object.overlap.at(actor_i));
    int res = 0;
    for (auto itr = difference_result.begin(); itr != difference_result.end(); itr++) {
      res+= object.x_attribute.get_val(*itr);
    } 
    return(res);
  } else { 
    return(0);
  }
};
EFFECT_REGISTER("attribute_xy_alocal", ::xyz_stat_attribute_xy_nonb, "attribute_xy_alocal", 0);


auto xyz_stat_attribute_yz_nb= CHANGESTAT{
  if(mode == "y"){
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i)); 
    } else {   
      connections_of_i = connections_of_i_all;
    }    
    res = connections_of_i.size();
    return(res);
  } else if(mode == "z"){  
    return(object.y_attribute.get_val(actor_i) + object.y_attribute.get_val(actor_j));
  } else { 
    return(0);
  } 
};
EFFECT_REGISTER("attribute_yz_local", ::xyz_stat_attribute_yz_nb, "attribute_yz_local", 0);

auto xyz_stat_attribute_xz_nb= CHANGESTAT{
  if(mode == "x"){
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i= get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {   
      connections_of_i = connections_of_i_all;
    }   
    res = connections_of_i.size();
    return(res);
  } else if(mode == "z"){ 
    return(object.x_attribute.get_val(actor_i) + object.x_attribute.get_val(actor_j));
  } else { 
    return(0);
  } 
};
EFFECT_REGISTER("attribute_xz_local", ::xyz_stat_attribute_xz_nb, "attribute_xz_local", 0);

// 
// double xyz_stat_attribute_x_nb(XYZ_class &object,
//                             int &actor_i,
//                             int &actor_j,
//                             arma::mat &data,
//                             double &type,
//                             std::string &mode,
//                             bool &is_full_neighborhood){
//   if(mode == "x"){
//     int res = 0;
//     std::unordered_set<int>::iterator itr;
//     for (itr = object.overlap.at(actor_i).begin(); itr != object.overlap.at(actor_i).end(); itr++) {
//       res+= object.x_attribute.get_val(*itr);
//     }
//     return(res);
//   } else {
//     return(0);
//   }
// }
// 
// double xyz_stat_attribute_x_z_nb(XYZ_class &object,
//                                  int &actor_i,
//                                  int &actor_j,
//                                  arma::mat &data,
//                                  double &type,
//                                  std::string &mode,
//                                  bool &is_full_neighborhood){
//   if(mode == "x"){
//     int res = 0;
//     std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
//     std::unordered_set<int> connections_of_i;
//     
//     // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//     if(!is_full_neighborhood){
//       // Next we only want to get the connections within the same group
//       std::set_intersection(std::begin(connections_of_i_all), std::end(connections_of_i_all),
//                             std::begin(object.overlap.at(actor_i)), std::end(object.overlap.at(actor_i)),
//                             std::inserter(connections_of_i, std::begin(connections_of_i)));
//     } else {  
//       connections_of_i = connections_of_i_all;
//     }  
//     res = connections_of_i.size();
//     return(res);
//   } else if(mode == "z"){
//     return(object.x_attribute.get_val(actor_i) + object.x_attribute.get_val(actor_j));
//   } else {
//     return(0);
//   }
// } 

auto xyz_stat_edges_x_out_nb= CHANGESTAT{
  if(mode == "x"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {   
      connections_of_i = connections_of_i_all;
    }    
    return(connections_of_i.size());
  } else if(mode == "z"){  
    return(object.x_attribute.get_val(actor_i)*object.get_val_overlap(actor_i,actor_j));
  }else { 
    return(0);
  }
}; 
EFFECT_REGISTER("outedges_x_local", ::xyz_stat_edges_x_out_nb, "outedges_x_local", 0);

auto xyz_stat_edges_x_out_nonb= CHANGESTAT{
  if(mode == "x"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_set_difference(connections_of_i_all, object.overlap.at(actor_i));
    } else {   
      connections_of_i = connections_of_i_all;
    }    
    return(connections_of_i.size());
  } else if(mode == "z"){  
    return(object.x_attribute.get_val(actor_i)*(1-object.get_val_overlap(actor_i,actor_j)));
  }else { 
    return(0);
  }
};
EFFECT_REGISTER("outedges_x_alocal", ::xyz_stat_edges_x_out_nonb, "outedges_x_alocal", 0);


auto xyz_stat_edges_x_out= CHANGESTAT{
  if(mode == "x"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    return(connections_of_i_all.size());
  } else if(mode == "z"){ 
    return(object.x_attribute.get_val(actor_i));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("outedges_x_global", ::xyz_stat_edges_x_out, "outedges_x_global", 0);

auto xyz_stat_edges_x_in_nb= CHANGESTAT{
  if(mode == "x"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
      
    } else {   
      connections_of_i = connections_of_i_all;
    }   
    return(connections_of_i.size());
  } else if(mode == "z"){ 
    return(object.x_attribute.get_val(actor_j)*object.overlap.at(actor_i).count(actor_j));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("inedges_x_local", ::xyz_stat_edges_x_in_nb, "inedges_x_local", 0);

auto xyz_stat_edges_x_in_nonb= CHANGESTAT{
  if(mode == "x"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_set_difference(connections_of_i_all, object.overlap.at(actor_i));
    } else {   
      return(0.0);
    }   
    return(connections_of_i.size());
  } else if(mode == "z"){ 
    return(object.x_attribute.get_val(actor_j)*(1-object.overlap.at(actor_i).count(actor_j)));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("inedges_x_alocal", ::xyz_stat_edges_x_in_nonb, "inedges_x_alocal", 0);

auto xyz_stat_edges_x_in= CHANGESTAT{
  if(mode == "x"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    return(connections_of_i_all.size());
  } else if(mode == "z"){ 
    return(object.x_attribute.get_val(actor_j));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("inedges_x_global", ::xyz_stat_edges_x_in, "inedges_x_global", 0);


auto xyz_stat_edges_y_out_nb= CHANGESTAT{
  if(mode == "y"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i)); 
    } else {   
      connections_of_i = connections_of_i_all;
    }   
    return(connections_of_i.size());
  } else if(mode == "z"){ 
    return(object.y_attribute.get_val(actor_i)*object.overlap.at(actor_i).count(actor_j));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("outedges_y_local", ::xyz_stat_edges_y_out_nb, "outedges_y_local", 0);

auto xyz_stat_edges_y_out_nonb= CHANGESTAT{
  if(mode == "y"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_set_difference(connections_of_i_all, object.overlap.at(actor_i));
      return(0.0);
    }   
    return(connections_of_i.size());
  } else if(mode == "z"){ 
    return(object.y_attribute.get_val(actor_i)*(1-object.overlap.at(actor_i).count(actor_j)));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("outedges_y_alocal", ::xyz_stat_edges_y_out_nonb, "outedges_y_alocal", 0);


auto xyz_stat_edges_y_out = CHANGESTAT{
  if(mode == "y"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    return(connections_of_i_all.size());
  } else if(mode == "z"){ 
    return(object.y_attribute.get_val(actor_i));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("outedges_y_global", ::xyz_stat_edges_y_out, "outedges_y_global", 0);


auto xyz_stat_edges_y_in_nb= CHANGESTAT{
  if(mode == "y"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all,object.overlap.at(actor_i));
    } else {    
      connections_of_i = connections_of_i_all;
    }    
    return(connections_of_i.size());
  } else if(mode == "z"){ 
    return(object.y_attribute.get_val(actor_j)*object.overlap.at(actor_i).count(actor_j));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("inedges_y_local", ::xyz_stat_edges_y_in_nb, "inedges_y_local", 0);

auto xyz_stat_edges_y_in_nonb= CHANGESTAT{
  if(mode == "y"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_set_difference(connections_of_i_all, object.overlap.at(actor_i));
    } else {    
      return(0.0);
    }    
    return(connections_of_i.size());
  } else if(mode == "z"){ 
    return(object.y_attribute.get_val(actor_j)*(1-object.overlap.at(actor_i).count(actor_j)));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("inedges_y_alocal", ::xyz_stat_edges_y_in_nonb, "inedges_y_alocal", 0);

auto xyz_stat_edges_y_in= CHANGESTAT{
  if(mode == "y"){
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    return(connections_of_i_all.size());
  } else if(mode == "z"){ 
    return(object.y_attribute.get_val(actor_j));
  }else {
    return(0);
  }
};
EFFECT_REGISTER("inedges_y_global", ::xyz_stat_edges_y_in, "inedges_y_global", 0);


auto xyz_stat_attribute_x= CHANGESTAT{
  if(mode == "x"){
    return(1);
  } else {
    return(0);
  }
};
EFFECT_REGISTER("attribute_x", ::xyz_stat_attribute_x, "attribute_x", 0);

auto xyz_stat_attribute_y= CHANGESTAT{
  if(mode == "y"){
    return(1);
  } else {
    return(0);
  }
};
EFFECT_REGISTER("attribute_y", ::xyz_stat_attribute_y, "attribute_y", 0);


auto xyz_stat_interaction_edges= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    if(object.get_val_overlap(actor_i, actor_j)){
      return(object.x_attribute.get_val(actor_i)*object.y_attribute.get_val(actor_j)+
             object.x_attribute.get_val(actor_j)*object.y_attribute.get_val(actor_i));  
    } else {
      return(0);
    }
    
  } else if (mode == "x"){
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections with overlap between them
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
      
    } else {
      connections_of_i = connections_of_i_all;
    }
    
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      res+= object.y_attribute.get_val(*itr);
      // }
    } 
    // Rcout << res << std::endl;
    return(res);
  } else{
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all,object.overlap.at(actor_i)); 
      
    } else {
      connections_of_i = connections_of_i_all;
    }
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // if(*itr != actor_j){
      res+= object.x_attribute.get_val(*itr);  
      // }
    }
    // Rcout << res << std::endl;
    return(res);
  }
};
EFFECT_REGISTER("spillover_xy_symm", ::xyz_stat_interaction_edges, "spillover_xy_symm", 0);

// cov_i *cov_j * z_ij*c_ij
auto xyz_stat_interaction_edges_cov= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    if(object.get_val_overlap(actor_i, actor_j)){
      
      return(data.at(0, actor_i-1)*object.y_attribute.get_val(actor_j)+
             data.at(0, actor_j-1)*object.y_attribute.get_val(actor_i));  
    } else {
      return(0);
    }
    
  } else if (mode == "x"){
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    int res = 0;
    return(res);
  } else{
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all,object.overlap.at(actor_i));
      
    } else { 
      connections_of_i = connections_of_i_all;
    }
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // if(*itr != actor_j){
      res+= data.at(0, *itr-1);  
      // }
    } 
    // Rcout << res << std::endl;
    return(res);
  }
};
EFFECT_REGISTER("spillover_yc_symm", ::xyz_stat_interaction_edges_cov, "spillover_yc_symm", 0);

auto xyz_stat_interaction_edges_xy= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    if(object.get_val_overlap(actor_i, actor_j)){
      return(object.x_attribute.get_val(actor_i)*object.y_attribute.get_val(actor_j));  
    } else { 
      return(0);
    } 
    
  } else if (mode == "x"){ 
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections with overlap between them
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {
      connections_of_i = connections_of_i_all;
    }
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      res+= object.y_attribute.get_val(*itr); 
    } 
    // Rcout <<  res << std::endl;
    return(res);
  } else{ 
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else { 
      connections_of_i = connections_of_i_all;
    } 
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // if(*itr != actor_j){
      res+= object.x_attribute.get_val(*itr);   
      // }
    } 
    // Rcout << res << std::endl;
    return(res);
  }
};
EFFECT_REGISTER("spillover_xy", ::xyz_stat_interaction_edges_xy, "spillover_xy", 0);
// 
// auto xyz_stat_interaction_edges_directed= CHANGESTAT{
//   if(mode == "z"){
//     // What to do if the network change stat is desired
//     // z_ij from 0 -> 1
//     if(object.get_val_overlap(actor_i, actor_j)){
//       return(object.x_attribute.get_val(actor_i)*object.y_attribute.get_val(actor_j));  
//     } else { 
//       return(0);
//     }  
//     
//   } else if (mode == "x"){  
//     // What to do if the attribute change stat is wanted
//     // x_i from 0 -> 1
//     int res = 0;
//     std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
//     std::unordered_set<int> connections_of_i;
//     
//     // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//     if(!is_full_neighborhood){
//       // Next we only want to get the connections with overlap between them
//       std::set_intersection(std::begin(connections_of_i_all), std::end(connections_of_i_all),
//                             std::begin(object.overlap.at(actor_i)), std::end(object.overlap.at(actor_i)),
//                             std::inserter(connections_of_i, std::begin(connections_of_i)));
//     } else { 
//       connections_of_i = connections_of_i_all;
//     } 
//     
//     std::unordered_set<int>::iterator itr;
//     for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
//       res+= object.y_attribute.get_val(*itr); 
//     }  
//     
//     std::unordered_set<int> connections_of_j_all =  object.z_network.adj_list.at(actor_j);
//     std::unordered_set<int> connections_of_j;
//     
//     // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//     if(!is_full_neighborhood){
//       // Next we only want to get the connections with overlap between them
//       std::set_intersection(std::begin(connections_of_j_all), std::end(connections_of_j_all),
//                             std::begin(object.overlap.at(actor_j)), std::end(object.overlap.at(actor_j)),
//                             std::inserter(connections_of_j, std::begin(connections_of_j)));
//     } else { 
//       connections_of_j = connections_of_j_all;
//     } 
//     
//     for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
//       res+= object.y_attribute.get_val(*itr); 
//     }  
//     
//     return(res);
//   } else{  
//     // What to do if the attribute change stat is wanted
//     // y_i from 0 -> 1
//     int res = 0;
//     std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
//     std::unordered_set<int> connections_of_i;
//     
//     // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//     if(!is_full_neighborhood){
//       // Next we only want to get the connections within the same group
//       std::set_intersection(std::begin(connections_of_i_all), std::end(connections_of_i_all),
//                             std::begin(object.overlap.at(actor_i)), std::end(object.overlap.at(actor_i)),
//                             std::inserter(connections_of_i, std::begin(connections_of_i)));
//     } else {  
//       connections_of_i = connections_of_i_all;
//     }  
//     std::unordered_set<int>::iterator itr;
//     for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
//       // if(*itr != actor_j){
//       res+= object.x_attribute.get_val(*itr);   
//       // }
//     }  
//     
//     std::unordered_set<int> connections_of_j_all =  object.z_network.adj_list_in.at(actor_j);
//     std::unordered_set<int> connections_of_j;
//     
//     // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//     if(!is_full_neighborhood){
//       // Next we only want to get the connections within the same group
//       std::set_intersection(std::begin(connections_of_j_all), std::end(connections_of_j_all),
//                             std::begin(object.overlap.at(actor_j)), std::end(object.overlap.at(actor_j)),
//                             std::inserter(connections_of_j, std::begin(connections_of_j)));
//     } else {  
//       connections_of_j = connections_of_j_all;
//     }  
//     for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
//       // if(*itr != actor_j){
//       res+= object.x_attribute.get_val(*itr);   
//       // }
//     }  
//     return(res);
//   }
// };
// EFFECT_REGISTER("interaction_edges_directed", ::xyz_stat_interaction_edges_directed, "interaction_edges_directed", 0);


auto xyz_stat_interaction_edges_y_cov= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    if(object.get_val_overlap(actor_i, actor_j)){
      return(data.at(0,actor_i-1)*object.y_attribute.get_val(actor_j));  
    } else { 
      return(0);
    }  
  } else if (mode == "x"){  
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    // Do nothing!
    int res = 0;
    return(res);
  } else{  
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {  
      connections_of_i = connections_of_i_all;
    }  
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // if(*itr != actor_j){
      res+= data.at(0,*itr-1);   
      // }
    }  
    // Rcout << res << std::endl;
    return(res);
  }
};
EFFECT_REGISTER("spillover_yc", ::xyz_stat_interaction_edges_y_cov, "spillover_yc", 0);



auto xyz_stat_interaction_edges_xy_scaled= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    
    
    if(object.get_val_overlap(actor_i, actor_j)){
      if(object.neighborhood.at(actor_i).size()>0){
        return(object.x_attribute.get_val(actor_i)*object.y_attribute.get_val(actor_j)/object.neighborhood.at(actor_i).size());  
      } else {
        return(0);
      }
    } else {  
      return(0);
    }  
    
  } else if (mode == "x"){  
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    double res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections with overlap between them
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {
      connections_of_i = connections_of_i_all;
    }
    
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      res+= object.y_attribute.get_val(*itr); 
    } 
    // Rcout <<  res << std::endl;
    // Rcout <<  res/object.neighborhood.at(actor_i).size() << std::endl;
    if(object.neighborhood.at(actor_i).size()>0){
      return(res/object.neighborhood.at(actor_i).size());  
    } else {
      return(0);
    }
  } else{ 
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    double res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i, overlap;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {  
      connections_of_i = connections_of_i_all;
    }  
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // if(*itr != actor_j){
      res+= object.x_attribute.get_val(*itr);   
      // }
    }  
    // Rcout <<  res << std::endl;
    // Rcout <<  res/object.neighborhood.at(actor_i).size() << std::endl;
    if(object.neighborhood.at(actor_i).size()>0){
      return(res/object.neighborhood.at(actor_i).size());  
    } else {
      return(0);
    }
    
  }
};
EFFECT_REGISTER("spillover_xy_scaled", ::xyz_stat_interaction_edges_xy_scaled, "spillover_xy_scaled", 0);


auto xyz_stat_interaction_edges_yx= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    if(object.get_val_overlap(actor_i, actor_j)){
      return(object.x_attribute.get_val(actor_j)*object.y_attribute.get_val(actor_i));  
    } else { 
      return(0);
    } 
    
  } else if (mode == "x"){ 
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections with overlap between them
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {
      connections_of_i = connections_of_i_all;
    }
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // Rcout << "Attribute of actor";
      // Rcout << *itr << std::endl;
      // Rcout << object.attribute.get_val(*itr) << std::endl;
      // if(*itr != actor_j){
      res+= object.y_attribute.get_val(*itr); 
      // }
    } 
    // Rcout <<  res << std::endl;
    return(res);
  } else{ 
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {  
      connections_of_i = connections_of_i_all;
    } 
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // if(*itr != actor_j){
      res+= object.x_attribute.get_val(*itr);   
      // }
    } 
    // Rcout << res << std::endl;
    return(res);
  }
};
EFFECT_REGISTER("spillover_yx", ::xyz_stat_interaction_edges_yx, "spillover_yx", 0);



// auto xyz_stat_interaction_dyads= CHANGESTAT{
//   if (mode == "x"){
//     // What to do if the attribute change stat is wanted
//     // x_i from 0 -> 1
//     int res = 0;
//     std::unordered_set<int>::iterator itr;
//     for (itr = object.overlap.at(actor_i).begin(); itr != object.overlap.at(actor_i).end(); itr++) {
//       res+= object.y_attribute.get_val(*itr);
//     }  
//     // Rcout << res << std::endl;
//     return(res);
//   } else if (mode == "y"){ 
//     // What to do if the attribute change stat is wanted
//     // y_i from 0 -> 1
//     int res = 0;
//     std::unordered_set<int>::iterator itr;
//     for (itr = object.overlap.at(actor_i).begin(); itr != object.overlap.at(actor_i).end(); itr++) {
//       res+= object.x_attribute.get_val(*itr);
//     } 
//     // Rcout << res << std::endl;
//     return(res);
//   }
//   // Rcout << "Here we are!" << std::endl;
//   return(0);
// };
// EFFECT_REGISTER("interaction_dyads", ::xyz_stat_interaction_dyads, "interaction_dyads", 0);

auto xyz_stat_attribute_xy= CHANGESTAT{
  // Rcout <<  "Starting" << std::endl;
  if (mode == "x"){
    // Rcout <<  actor_i << std::endl;
    // Rcout <<  object.y_attribute.get_val(actor_i) << std::endl;
    return(object.y_attribute.get_val(actor_i));
  } else if (mode == "y") {  
    // Rcout <<  actor_i << std::endl;
    // Rcout <<  object.x_attribute.attribute.size() << std::endl;
    return(object.x_attribute.get_val(actor_i));
  } else {
    return(0);
  }
};
EFFECT_REGISTER("attribute_xy_global", ::xyz_stat_attribute_xy, "attribute_xy_global", 0);

// x_i*x_j*z_ij * c_ij
auto xyz_stat_matching_edges_x= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // z_ij from 0 -> 1
    // The change statistic will be x_i*x_j if actors i and j are in the same neighborhood
    bool same_group = object.get_val_overlap(actor_i, actor_j);
    // If the neighborhood is a full graph all actors are within the same group
    if(same_group) {
      // Rcout << "Change of an network!" << std::endl;
      // Rcout << actor_i << std::endl;
      // Rcout << actor_j << std::endl;
      // Rcout << object.attribute.get_val(actor_i)*object.attribute.get_val(actor_j) << std::endl;
      return(object.x_attribute.get_val(actor_i)*object.x_attribute.get_val(actor_j));
    } else {
      return(0);
    }
  } else if (mode == "x"){
    // What to do if the attribute change stat is wanted
    // x_i from 0 -> 1
    // The change statistic will be sum_{h with h and i being in the same neighborhood}x_h z_{h,i} 
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {
      connections_of_i = connections_of_i_all;
    }
    for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // Rcout << "Attribute of actor";
      // Rcout << *itr << std::endl;
      // Rcout << object.attribute.get_val(*itr) << std::endl;
      // if(*itr != actor_j){
      res+= object.x_attribute.get_val(*itr);
      // }
    }
    if(object.z_network.directed){
      std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
      std::unordered_set<int> connections_of_i;
      
      // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
      if(!is_full_neighborhood){
        // Next we only want to get the connections within the same group
        connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
      } else {
        connections_of_i = connections_of_i_all;
      }
      for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
        // Rcout << "Attribute of actor";
        // Rcout << *itr << std::endl;
        // Rcout << object.attribute.get_val(*itr) << std::endl;
        // if(*itr != actor_j){
        res+= object.x_attribute.get_val(*itr);
        // }
      }
    }
    // Rcout << res << std::endl;
    return(res);
  } else {
    return(0);
  }
};
EFFECT_REGISTER("spillover_xx", ::xyz_stat_matching_edges_x, "spillover_xx", 0);

// y_i*y_j*z_ij * c_ij
auto xyz_stat_matching_edges_y= CHANGESTAT{
  if(mode == "z"){
    // What to do if the network change stat is desired
    // x_ij from 0 -> 1
    bool same_group = object.get_val_overlap(actor_i, actor_j);
    // If the neighborhood is a full graph all actors are within the same group
    if(same_group) {
      // Rcout << "Change of an network!" << std::endl;
      // Rcout << actor_i << std::endl;
      // Rcout << actor_j << std::endl;
      // Rcout << object.attribute.get_val(actor_i)*object.attribute.get_val(actor_j) << std::endl;
      return(object.y_attribute.get_val(actor_i)*object.y_attribute.get_val(actor_j));
    } else {
      return(0);
    }
  } else if (mode == "y"){
    // What to do if the attribute change stat is wanted
    // y_i from 0 -> 1
    int res = 0;
    std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
    std::unordered_set<int> connections_of_i;
    
    // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
    if(!is_full_neighborhood){
      // Next we only want to get the connections within the same group
      connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
    } else {
      connections_of_i = connections_of_i_all;
    }
    for ( auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
      // Rcout << "Attribute of actor";
      // Rcout << *itr << std::endl;
      // Rcout << object.attribute.get_val(*itr) << std::endl;
      // if(*itr != actor_j){
      res+= object.y_attribute.get_val(*itr);
      // }
      
    }
    if(object.z_network.directed){
      std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
      std::unordered_set<int> connections_of_i;
      
      // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
      if(!is_full_neighborhood){
        // Next we only want to get the connections within the same group
        connections_of_i = get_intersection(connections_of_i_all, object.overlap.at(actor_i));
      } else {
        connections_of_i = connections_of_i_all;
      }
      for (auto itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
        // Rcout << "Attribute of actor";
        // Rcout << *itr << std::endl;
        // Rcout << object.attribute.get_val(*itr) << std::endl;
        // if(*itr != actor_j){
        res+= object.y_attribute.get_val(*itr);
        // }
        
      }
    }
    
    // Rcout << res << std::endl;
    return(res);
  } else {
    return(0);
  }
};
EFFECT_REGISTER("spillover_yy", ::xyz_stat_matching_edges_y, "spillover_yy", 0);


// // Are there common partners between actor_i and actor_j? -> now this should be the number of transitive edges
// auto xyz_stat_transitive_edges= CHANGESTAT{
//   if(mode == "z"){
//     // If the two actors are not within the same group the change stat is always 0
//     int res = 0;
//     
//     // object.print();
//     // Rcout <<actor_j << std::endl;
//     bool same_group = object.get_val_overlap(actor_i, actor_j);
//     if(same_group){
//       // First we look at the case for directed connections 
//       // (a transitive connection is here z_i,j such that a z_i,h and z_h,j exists with h in N_i and in N_j)
//       if(object.z_network.directed){
//         std::unordered_set<int> intersect_group, intersect_nb;
//         // Get the intersection group -> this encompasses all actors that are in the overlapping neighborhoods of i and j 
//         std::set_intersection(std::begin(object.overlap.at(actor_i)),
//                               std::end(object.overlap.at(actor_i)),
//                               std::begin(object.overlap.at(actor_j)),
//                               std::end(object.overlap.at(actor_j)),
//                               std::inserter(intersect_group, std::begin(intersect_group)));
//         
//         std::set_intersection(std::begin(object.neighborhood.at(actor_i)),
//                               std::end(object.neighborhood.at(actor_i)),
//                               std::begin(object.neighborhood.at(actor_j)),
//                               std::end(object.neighborhood.at(actor_j)),
//                               std::inserter(intersect_nb, std::begin(intersect_nb)));
//         // Step 1: Get all connections of actor_i
//         std::unordered_set<int> out_connections_of_i_all =  object.z_network.adj_list.at(actor_i);
//         std::unordered_set<int> in_connections_of_i_all =  object.z_network.adj_list_in.at(actor_i);
//         std::unordered_set<int> out_connections_of_i, out_connections_of_i_nb;
//         std::unordered_set<int> in_connections_of_i, in_connections_of_i_nb;
//         
//         // Step 2: Get all connections of actor_j
//         std::unordered_set<int> in_connections_of_j_all = object.z_network.adj_list_in.at(actor_j);
//         std::unordered_set<int> out_connections_of_j_all = object.z_network.adj_list.at(actor_j);
//         std::unordered_set<int> in_connections_of_j, in_connections_of_j_nb;
//         std::unordered_set<int> out_connections_of_j, out_connections_of_j_nb;
//         
//         // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//         if(!is_full_neighborhood){
//           // Next we only want to get the connections of i within the same group
//           std::set_intersection(std::begin(out_connections_of_i_all), std::end(out_connections_of_i_all),
//                                 std::begin(intersect_group), std::end(intersect_group),
//                                 std::inserter(out_connections_of_i, std::begin(out_connections_of_i)));
//           std::set_intersection(std::begin(in_connections_of_i_all), std::end(in_connections_of_i_all),
//                                 std::begin(intersect_group), std::end(intersect_group),
//                                 std::inserter(in_connections_of_i, std::begin(in_connections_of_i)));
//           std::set_intersection(std::begin(out_connections_of_j_all), std::end(out_connections_of_j_all),
//                                 std::begin(intersect_group), std::end(intersect_group),
//                                 std::inserter(out_connections_of_j, std::begin(out_connections_of_j)));
//           std::set_intersection(std::begin(in_connections_of_j_all), std::end(in_connections_of_j_all),
//                                 std::begin(intersect_group), std::end(intersect_group),
//                                 std::inserter(in_connections_of_j, std::begin(in_connections_of_j)));
//           
//           std::set_intersection(std::begin(out_connections_of_i_all), std::end(out_connections_of_i_all),
//                                 std::begin(object.neighborhood.at(actor_i)), std::end(object.neighborhood.at(actor_i)),
//                                 std::inserter(out_connections_of_i_nb, std::begin(out_connections_of_i_nb)));
//           std::set_intersection(std::begin(in_connections_of_i_all), std::end(in_connections_of_i_all),
//                                 std::begin(object.neighborhood.at(actor_i)), std::end(object.neighborhood.at(actor_i)),
//                                 std::inserter(in_connections_of_i_nb, std::begin(in_connections_of_i_nb)));
//           std::set_intersection(std::begin(out_connections_of_j_all), std::end(out_connections_of_j_all),
//                                 std::begin(object.neighborhood.at(actor_j)), std::end(object.neighborhood.at(actor_j)),
//                                 std::inserter(out_connections_of_j_nb, std::begin(out_connections_of_j_nb)));
//           std::set_intersection(std::begin(in_connections_of_j_all), std::end(in_connections_of_j_all),
//                                 std::begin(object.neighborhood.at(actor_j)), std::end(object.neighborhood.at(actor_j)),
//                                 std::inserter(in_connections_of_j_nb, std::begin(in_connections_of_j_nb)));
//           
//         } else {
//           out_connections_of_i = out_connections_of_i_all;
//           in_connections_of_i = in_connections_of_i_all;
//           out_connections_of_j = out_connections_of_j_all;
//           in_connections_of_j = in_connections_of_j_all;
//           
//           out_connections_of_i_nb = out_connections_of_i_all;
//           in_connections_of_i_nb = in_connections_of_i_all;
//           out_connections_of_j_nb = out_connections_of_j_all;
//           in_connections_of_j_nb = in_connections_of_j_all;
//         }
//         
//         // If i and j are already connected exclude them from the respective sets
//         in_connections_of_j.erase(actor_i);
//         out_connections_of_i.erase(actor_j);
//         // intersection_ij is the set of shared connections of actors i and j
//         // intersection is the set of shared connections of actors *itr and j or i
//         // tmpis the set of connections or *itr
//         std::unordered_set<int> intersection_i_out_j_in_nb, intersection_i_in_j_in, 
//         intersection_i_out_j_out, intersection_i_out_j_in,intersection;
//         // Step 3: Get the intersection of connections to i and j -> case a)
//         std::set_intersection(std::begin(in_connections_of_i), std::end(in_connections_of_i),
//                               std::begin(in_connections_of_j), std::end(in_connections_of_j),
//                               std::inserter(intersection_i_in_j_in, std::begin(intersection_i_in_j_in)));
//         std::set_intersection(std::begin(out_connections_of_i), std::end(out_connections_of_i),
//                               std::begin(out_connections_of_j), std::end(out_connections_of_j),
//                               std::inserter(intersection_i_out_j_out, std::begin(intersection_i_out_j_out)));
//         std::set_intersection(std::begin(out_connections_of_i), std::end(out_connections_of_i),
//                               std::begin(in_connections_of_j), std::end(in_connections_of_j),
//                               std::inserter(intersection_i_out_j_in, std::begin(intersection_i_out_j_in)));
//         std::set_intersection(std::begin(intersection_i_out_j_in), std::end(intersection_i_out_j_in),
//                               std::begin(intersect_nb), std::end(intersect_nb),
//                               std::inserter(intersection_i_out_j_in_nb, std::begin(intersection_i_out_j_in_nb)));
//         
//         if(intersection_i_out_j_in_nb.size()>0){
//           // If actors i and j have common partners in the same neighborhood then the change statistic is at least 1
//           // (since 0->1 would make z_ij a transitive edge)
//           res += 1;
//         }
//         //  Going through all actors h that send to i and j
//         // the change is one if there is no other two path (with z_{h,k} = z_{k,j}) between h and j than the one through i
//         if(object.get_val_neighborhood(actor_j,actor_i) && (intersection_i_in_j_in.size()>0)){
//           // Rcout << "The intersection is larger than 0 of size " + std::to_string(intersection_ij.size())<< std::endl;
//           // For all the partners of i and j we need to check whether they have a connection
//           // Step 4.1: Get all connections from actors in this intersection
//           std::unordered_set<int>::iterator itr;
//           for (itr = intersection_i_in_j_in.begin(); itr != intersection_i_in_j_in.end(); itr++) {
//             // Save the connections of *itr (which is an iterator going through all common partners within the same neighborhood of i and j)
//             if(object.get_val_neighborhood(*itr,actor_i)){
//               std::unordered_set<int> tmp_out, tmp_out2, tmp_out3;
//               // Possible other connections
//               tmp_out = object.z_network.adj_list.at(*itr);
//               // Those are all the outgoing ties within the neighbohood of *itr
//               std::set_intersection(std::begin(tmp_out), std::end(tmp_out),
//                                     std::begin(object.neighborhood.at(*itr)), std::end(object.neighborhood.at(*itr)),
//                                     std::inserter(tmp_out3, std::begin(tmp_out3)));
//               // Transitive edges can be observed by any two path going through the common neighborhood of h and j
//               // Only look at the connections of *itr that are in the neighborhood of h and j
//               std::set_intersection(std::begin(tmp_out3), std::end(tmp_out3),
//                                     std::begin(object.neighborhood.at(actor_j)), std::end(object.neighborhood.at(actor_j)),
//                                     std::inserter(tmp_out2, std::begin(tmp_out2)));
//               // i should not be included here
//               tmp_out2.erase(actor_i);
//               intersection.clear();
//               std::set_intersection(std::begin(tmp_out2), std::end(tmp_out2),
//                                     std::begin(in_connections_of_j_nb), std::end(in_connections_of_j_nb),
//                                     std::inserter(intersection, std::begin(intersection)));
//               // If i and itr* have no common partner, the edge between i and *itr would have become a transitive edge
//               if(intersection.size()==0){
//                 res += 1;
//               }  
//             }
//           }
//         }
//         //  Going through all actors h that got sent to from i and j
//         // the change is one if there is no other two path (with z_{i,k} = z_{k,h} = 1) between h and j than the one through j
//         if(object.get_val_neighborhood(actor_i,actor_j) &&(intersection_i_out_j_out.size()>0)){
//           // Rcout << "The intersection is larger than 0 of size " + std::to_string(intersection_ij.size())<< std::endl;
//           // For all the partners of i and j we need to check whether they have a connection
//           // Step 4.1: Get all connections from actors in this intersection
//           std::unordered_set<int>::iterator itr;
//           for (itr = intersection_i_out_j_out.begin(); itr != intersection_i_out_j_out.end(); itr++) {
//             if(object.get_val_neighborhood(*itr,actor_j)){
//               std::unordered_set<int> tmp_in, tmp_in2, tmp_in3;
//               tmp_in = object.z_network.adj_list_in.at(*itr);
//               // Only look at the connections of *itr that are in its neighborhood
//               std::set_intersection(std::begin(tmp_in), std::end(tmp_in),
//                                     std::begin(object.neighborhood.at(*itr)), std::end(object.neighborhood.at(*itr)),
//                                     std::inserter(tmp_in3, std::begin(tmp_in3)));
//               // Transitive edges can be observed by any two path going through the common neighborhood of h and j
//               // Only look at the connections of *itr that are in the neighborhood of h and j
//               std::set_intersection(std::begin(tmp_in3), std::end(tmp_in3),
//                                     std::begin(object.neighborhood.at(actor_i)), std::end(object.neighborhood.at(actor_i)),
//                                     std::inserter(tmp_in2, std::begin(tmp_in2)));
//               tmp_in2.erase(actor_j);
//               intersection.clear();
//               std::set_intersection(std::begin(tmp_in2), std::end(tmp_in2),
//                                     std::begin(out_connections_of_i_nb), std::end(out_connections_of_i_nb),
//                                     std::inserter(intersection, std::begin(intersection)));
//               // If i and itr* have no common partner, the edge between i and *itr would have become a transitive edge
//               if(intersection.size()==0){
//                 res += 1;
//               }  
//             }
//           }
//         }
//       } else {
//         // Here we specify what happends in the undirected case
//         std::unordered_set<int> intersect_group, common_nb;
//         // Step 1: Get all connections of actor_i and actor_j
//         std::unordered_set<int> connections_of_i_all =  object.z_network.adj_list.at(actor_i);
//         std::unordered_set<int> connections_of_j_all = object.z_network.adj_list.at(actor_j);
//         std::unordered_set<int> connections_of_i;
//         std::unordered_set<int> connections_of_j;
//         // Get the intersection group of actors h that are in the overlapping neighborhoods of actor i (N_i \cap N_h \neq empty) and j (N_j \cap N_h \neq empty)
//         std::set_intersection(std::begin(object.overlap.at(actor_i)),
//                               std::end(object.overlap.at(actor_i)),
//                               std::begin(object.overlap.at(actor_j)),
//                               std::end(object.overlap.at(actor_j)),
//                               std::inserter(intersect_group, std::begin(intersect_group)));
//         
//         // We are only interesten in the connection within intersect_group
//         // If there is no full neighborhood we need to cut the connections of i to only include other actors within the same neighborhood
//         if(!is_full_neighborhood){
//           // Next we only want to get the connections of i within the same group
//           std::set_intersection(std::begin(connections_of_i_all), std::end(connections_of_i_all),
//                                 std::begin(intersect_group), std::end(intersect_group),
//                                 std::inserter(connections_of_i, std::begin(connections_of_i)));
//           // Next we only want to get the connections within the same group
//           std::set_intersection(std::begin(connections_of_j_all), std::end(connections_of_j_all),
//                                 std::begin(intersect_group), std::end(intersect_group),
//                                 std::inserter(connections_of_j, std::begin(connections_of_j)));
//         } else {
//           connections_of_i = connections_of_i_all;
//           connections_of_j = connections_of_j_all;
//           
//         }
//         connections_of_j.erase(actor_i);
//         connections_of_i.erase(actor_j);
//         
//         // Get the intersection group of actors h that in the neighborhoods of i and j
//         std::set_intersection(std::begin(object.neighborhood.at(actor_i)),
//                               std::end(object.neighborhood.at(actor_i)),
//                               std::begin(object.neighborhood.at(actor_j)),
//                               std::end(object.neighborhood.at(actor_j)),
//                               std::inserter(common_nb, std::begin(common_nb)));
//         std::unordered_set<int> intersection_ij, intersection_ij_nb;
//         std::set_intersection(std::begin(connections_of_i), std::end(connections_of_i),
//                               std::begin(connections_of_j), std::end(connections_of_j),
//                               std::inserter(intersection_ij, std::begin(intersection_ij)));
//         std::set_intersection(std::begin(intersection_ij), std::end(intersection_ij),
//                               std::begin(common_nb), std::end(common_nb),
//                               std::inserter(intersection_ij_nb, std::begin(intersection_ij_nb)));
//         if(intersection_ij_nb.size()>0){
//           res += 1;
//         }
//         
//         // Step 3: Get the intersection of connections to i and j -> case a)
//         
//         if(object.get_val_neighborhood(actor_j,actor_i) && (intersection_ij.size()>0)){
//           // Step 4.1: Get all connections from actors in this intersection
//           std::unordered_set<int>::iterator itr;
//           for (itr = intersection_ij.begin(); itr != intersection_ij.end(); itr++) {
//             std::unordered_set<int> intersection, tmp, tmp3;
//             
//             // Save the connections of *itr (which is an iterator going through all common partners within the same neighborhood of i and j)
//             tmp = object.z_network.adj_list.at(*itr);
//             // Only look at the connections of *itr that are in its neighborhood
//             std::set_intersection(std::begin(tmp), std::end(tmp),
//                                   std::begin(object.neighborhood.at(*itr)), std::end(object.neighborhood.at(*itr)),
//                                   std::inserter(tmp3, std::begin(tmp3)));
//             // Go through the common partners and check if they have no other common partners
//             // Get the size of the intersection of tmp and connections_of_i
//             // -> case b) where *itr is connected to people that are connected to i
//             intersection.clear();
//             // If i and itr* have no common partner, the edge between i and *itr would have become a transitive edge
//             if(object.get_val_neighborhood(*itr,actor_i)){
//               std::unordered_set<int> tmp2;
//               std::set_intersection(std::begin(tmp3), std::end(tmp3),
//                                     std::begin(object.neighborhood.at(actor_j)), std::end(object.neighborhood.at(actor_j)),
//                                     std::inserter(tmp2, std::begin(tmp2)));
//               // tmp2 is the set of actors that is connected to *itr and in the neighborhood of actor j and *itr
//               tmp2.erase(actor_i);
//               std::set_intersection(std::begin(tmp2), std::end(tmp2),
//                                     std::begin(connections_of_j), std::end(connections_of_j),
//                                     std::inserter(intersection, std::begin(intersection)));
//               if(intersection.size()==0){
//                 res += 1;
//               }
//             }
//             
//             // Get the size of the intersection of tmp and connections_of_j -> case c)
//             intersection.clear();
//             
//             // If j and itr* have no common partner, the edge between j and *itr would have become a transitive edge
//             if(object.get_val_neighborhood(*itr,actor_j)){
//               std::unordered_set<int> tmp1;
//               std::set_intersection(std::begin(tmp3), std::end(tmp3),
//                                     std::begin(object.neighborhood.at(actor_i)), std::end(object.neighborhood.at(actor_i)),
//                                     std::inserter(tmp1, std::begin(tmp1)));
//               // tmp1 is the set of actors that is connected to *itr and in the neighborhood of actor i and *itr
//               tmp1.erase(actor_j);
//               std::set_intersection(std::begin(tmp1), std::end(tmp1),
//                                     std::begin(connections_of_i), std::end(connections_of_i),
//                                     std::inserter(intersection, std::begin(intersection)));
//               if(intersection.size()==0){
//                 res += 1;
//               }
//               intersection.clear();
//             }
//           }
//         }
//       }
//     }
//     // Rcout <<"Actor i: " + std::to_string(actor_i) + " Actor j: " + std::to_string(actor_j) + " Change: "+ std::to_string(res)<< std::endl;
//     return(res);
//   }
//   else{
//     // Rcout << "Triangle on attribute" << std::endl;
//     return(0);
//   }
// };
// EFFECT_REGISTER("transitive", ::xyz_stat_transitive_edges, "transitive", 0);


static inline bool has_alternative_h_to_j(
    int h, int excluded, int j,
    const XYZ_class &object,
    const std::unordered_set<int> &in_connections_of_j_nb,
    const std::unordered_set<int> &neighborhood_h,
    const std::unordered_set<int> &neighborhood_j
) {
  auto it_out_h = object.z_network.adj_list.find(h);
  if (it_out_h == object.z_network.adj_list.end()) return false;
  const auto &out_h = it_out_h->second;
  
  // iterate smaller: out_h vs neighborhood_h or neighborhood_j isn't required -
  // we iterate out_h (likely smaller than neighborhoods) and test membership
  for (int k : out_h) {
    if (k == excluded) continue;
    // must be in both neighborhoods and be an in-connection to j within j's nb
    if (neighborhood_h.find(k) != neighborhood_h.end() &&
        neighborhood_j.find(k) != neighborhood_j.end() &&
        in_connections_of_j_nb.find(k) != in_connections_of_j_nb.end()) {
      return true;
    }
  }
  return false;
}

static inline bool has_alternative_i_to_h(
    int i, int excluded, int h,
    const XYZ_class &object,
    const std::unordered_set<int> &out_connections_of_i_nb,
    const std::unordered_set<int> &neighborhood_h,
    const std::unordered_set<int> &neighborhood_i
) {
  // iterate in-connections of h (who -> h)
  auto it_in_h = object.z_network.adj_list_in.find(h);
  if (it_in_h == object.z_network.adj_list_in.end()) return false;
  const auto &in_h = it_in_h->second;
  
  for (int k : in_h) {
    if (k == excluded) continue;
    if (neighborhood_h.find(k) != neighborhood_h.end() &&
        neighborhood_i.find(k) != neighborhood_i.end() &&
        out_connections_of_i_nb.find(k) != out_connections_of_i_nb.end()) {
      return true;
    }
  }
  return false;
}


auto xyz_stat_transitive_edges = CHANGESTAT {
  if (mode != "z") return 0.0;
  // if(object.overlap.at(actor_i).count(actor_j) == 0.0) {
  //   // Rcout << "Actors " << actor_i << " and " << actor_j << object.overlap.at(actor_i).count(actor_j)<< std::endl;
  //   return 0.0;
  //   
  // }
  // Quick presence checks
  if ((object.z_network.directed && (object.z_network.adj_list_in.find(actor_i) == object.z_network.adj_list_in.end() ||
      object.z_network.adj_list_in.find(actor_j) == object.z_network.adj_list_in.end())))
  {
    return 0.0;
  }
  
  int res = 0;
  const auto &overlap_i = object.overlap.at(actor_i);
  if (!overlap_i.count(actor_j)) return 0.0; // same_group check
  
  const auto &neighborhood_i = object.neighborhood.at(actor_i);
  const auto &neighborhood_j = object.neighborhood.at(actor_j);
  const auto &out_i_all = object.z_network.adj_list.at(actor_i);
  const auto &out_j_all = object.z_network.adj_list.at(actor_j);
  
  // Precompute intersect_group only if needed as unordered_set for O(1) lookups
  std::unordered_set<int> intersect_group_nb; // only used if !is_full_neighborhood
  if (!is_full_neighborhood) {
    // iterate smaller overlap set
    const auto &overlap_j = object.overlap.at(actor_j);
    const auto *small_ov = (&overlap_i);
    const auto *large_ov = (&overlap_j);
    if (overlap_j.size() < overlap_i.size()) { small_ov = &overlap_j; large_ov = &overlap_i; }
    intersect_group_nb.reserve(std::min<size_t>(small_ov->size(), large_ov->size()));
    for (int h : *small_ov) {
      if (large_ov->find(h) != large_ov->end()) intersect_group_nb.insert(h);
    }
  }
  
  if (object.z_network.directed) {
    const auto &in_i_all = object.z_network.adj_list_in.at(actor_i);
    const auto &in_j_all = object.z_network.adj_list_in.at(actor_j);
    
    // Build fast membership sets for the neighborhood-filtered in/out connections
    std::unordered_set<int> in_connections_of_j_nb;
    in_connections_of_j_nb.reserve(in_j_all.size());
    for (int n : in_j_all) if (neighborhood_j.find(n) != neighborhood_j.end()) in_connections_of_j_nb.insert(n);
    
    std::unordered_set<int> out_connections_of_i_nb;
    out_connections_of_i_nb.reserve(out_i_all.size());
    for (int n : out_i_all) if (neighborhood_i.find(n) != neighborhood_i.end()) out_connections_of_i_nb.insert(n);
    
    // --- Logic Part 1: Simple Transitivity i -> h -> j ---
    // iterate smaller of neighborhood intersection candidates: neighborhood_i vs neighborhood_j
    // and count presence of i->h and h->j
    int simple_transitivity_count = 0;
    const auto *small_nb = (&neighborhood_i);
    const auto *large_nb = (&neighborhood_j);
    if (neighborhood_j.size() < neighborhood_i.size()) { small_nb = &neighborhood_j; large_nb = &neighborhood_i; }
    for (int h : *small_nb) {
      if (h == actor_i || h == actor_j) continue;
      if (large_nb->find(h) == large_nb->end()) continue; // not common neighbor
      // check i->h and h->j using adjacency sets
      if (out_i_all.find(h) != out_i_all.end() && in_j_all.find(h) != in_j_all.end()) {
        simple_transitivity_count = 1; // we only need >0
        break;
      }
    }
    if (simple_transitivity_count) res += 1;
    
    // --- Logic Part 2: uniqueness check for h -> i -> j path (h != j) ---
    bool check_cond_part2 = neighborhood_j.find(actor_i) != neighborhood_j.end();
    
    if (check_cond_part2) {
      // iterate smaller of in_i_all and in_j_all to find common in-neighbors h
      const auto *small_in = &in_i_all;
      const auto *large_in = &in_j_all;
      if (in_j_all.size() < in_i_all.size()) { small_in = &in_j_all; large_in = &in_i_all; }
      
      for (int h : *small_in) {
        if (h == actor_i || h == actor_j) continue;
        if (large_in->find(h) == large_in->end()) continue;
        if (!is_full_neighborhood && intersect_group_nb.find(h) == intersect_group_nb.end()) continue;
        // check object.get_val_neighborhood(h, actor_i) placeholder:
        if (object.neighborhood.at(h).find(actor_i) == object.neighborhood.at(h).end()) continue;
        
        // check uniqueness: no alternative path h -> k -> j where k != i
        if (!has_alternative_h_to_j(h, actor_i, actor_j, object, in_connections_of_j_nb,
                                    object.neighborhood.at(h), neighborhood_j)) {
          res += 1;
        }
      }
    }
    
    // --- Logic Part 3: uniqueness check for i -> j -> h path (h != i) ---
    bool check_cond_part3 = neighborhood_i.find(actor_j) != neighborhood_i.end();
    
    if (check_cond_part3) {
      // iterate smaller of out_i_all and out_j_all to find common out-neighbors h
      const auto *small_out = &out_i_all;
      const auto *large_out = &out_j_all;
      if (out_j_all.size() < out_i_all.size()) { small_out = &out_j_all; large_out = &out_i_all; }
      
      for (int h : *small_out) {
        if (h == actor_i || h == actor_j) continue;
        if (large_out->find(h) == large_out->end()) continue;
        if (!is_full_neighborhood && intersect_group_nb.find(h) == intersect_group_nb.end()) continue;
        if (object.neighborhood.at(h).find(actor_j) == object.neighborhood.at(h).end()) continue;
        
        // check uniqueness: no alternative path i -> k -> h where k != j
        if (!has_alternative_i_to_h(actor_i, actor_j, h, object, out_connections_of_i_nb,
                                    object.neighborhood.at(h), neighborhood_i)) {
          res += 1;
        }
      }
    }
  }
  else { // undirected
    // Find common neighbors h (i - h and j - h)
    const auto *small_conn = &out_i_all;
    const auto *large_conn = &out_j_all;
    if (out_j_all.size() < out_i_all.size()) { small_conn = &out_j_all; large_conn = &out_i_all; }
    
    int simple_triangle_count = 0;
    // We'll collect the relevant common neighbors into a small vector to reuse for uniqueness checks
    std::vector<int> common_neighbors;
    common_neighbors.reserve(std::min(small_conn->size(), large_conn->size()));
    
    for (int h : *small_conn) {
      if (h == actor_i || h == actor_j) continue;
      if (large_conn->find(h) == large_conn->end()) continue;
      if (!is_full_neighborhood && intersect_group_nb.find(h) == intersect_group_nb.end()) continue;
      common_neighbors.push_back(h);
      // check if h also in N_i intersect N_j (common_direct_neighborhood)
      if (neighborhood_i.find(h) != neighborhood_i.end() && neighborhood_j.find(h) != neighborhood_j.end()) {
        simple_triangle_count = 1;
      }
    }
    if (simple_triangle_count) res += 1;
    // Uniqueness checks if object.get_val_neighborhood(actor_j, actor_i)
    bool check_cond_part2 = neighborhood_j.find(actor_i) != neighborhood_j.end();
    if (check_cond_part2 && !common_neighbors.empty()) {
      // Build fast neighborhood-filtered adjacency sets for i and j
      std::unordered_set<int> connections_of_i_nb;
      connections_of_i_nb.reserve(out_i_all.size());
      for (int n : out_i_all) if (neighborhood_i.find(n) != neighborhood_i.end()) connections_of_i_nb.insert(n);
      std::unordered_set<int> connections_of_j_nb;
      connections_of_j_nb.reserve(out_j_all.size());
      for (int n : out_j_all) if (neighborhood_j.find(n) != neighborhood_j.end()) connections_of_j_nb.insert(n);
      
      for (int h : common_neighbors) {
        // check object.get_val_neighborhood(h, actor_i)
        if (object.neighborhood.at(h).find(actor_i) != object.neighborhood.at(h).end()) {
          if (!has_alternative_h_to_j(h, actor_i, actor_j, object, connections_of_j_nb,
                                      object.neighborhood.at(h), neighborhood_j)) {
            res += 1; 
          }
        }
        // now symmetric: check object.get_val_neighborhood(h, actor_j)
        if (object.neighborhood.at(h).find(actor_j) != object.neighborhood.at(h).end()) {
          if (!has_alternative_h_to_j(h, actor_j, actor_i, object, connections_of_i_nb,
                                      object.neighborhood.at(h), neighborhood_i)) {
            res += 1;
          }
        }
      }
    }
  } // end undirected
  
  return static_cast<double>(res);
};

EFFECT_REGISTER("transitive", ::xyz_stat_transitive_edges, "transitive", 0);

auto xyz_stat_nonisolates= CHANGESTAT{
  if(mode == "z"){ 
    arma::vec res(3);
    int degree_i,degree_j;  
    if(object.z_network.directed){
      degree_i = object.z_network.adj_list.at(actor_i).size() + 
        object.z_network.adj_list_in.at(actor_i).size();
      degree_j = object.z_network.adj_list.at(actor_j).size() + 
        object.z_network.adj_list_in.at(actor_j).size();
    } else {
      degree_i = object.z_network.adj_list.at(actor_i).size();
      degree_j = object.z_network.adj_list.at(actor_j).size();
    }
    // If the edge is already there, we need to substract one of the degrees
    if(object.z_network.get_val(actor_i, actor_j)){
      degree_i -= 1;
      degree_j -= 1;
    }
    return((degree_i == 0) +(degree_j == 0));
  }else {
    return(0);
  }  
};
EFFECT_REGISTER("nonisolates", ::xyz_stat_nonisolates, "nonisolates", 0);

auto xyz_stat_isolates= CHANGESTAT{
  if(mode == "z"){ 
    arma::vec res(3);
    int degree_i,degree_j;  
    if(object.z_network.directed){
      degree_i = object.z_network.adj_list.at(actor_i).size() + 
        object.z_network.adj_list_in.at(actor_i).size();
      degree_j = object.z_network.adj_list.at(actor_j).size() + 
        object.z_network.adj_list_in.at(actor_j).size();
    } else {
      degree_i = object.z_network.adj_list.at(actor_i).size();
      degree_j = object.z_network.adj_list.at(actor_j).size();
    } 
    // If the edge is already there, we need to substract one of the degrees
    if(object.z_network.get_val(actor_i, actor_j)){
      degree_i -= 1;
      degree_j -= 1;
    }
    return(-(degree_i == 0)-(degree_j == 0));
  }else {
    return(0);
  }  
};
EFFECT_REGISTER("isolates", ::xyz_stat_isolates, "isolates", 1.0);

auto xyz_stat_gwesp_local_ITP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all ISP of i and j 
    std::unordered_set<int> itp_ij = object.get_common_partners_nb(actor_i, actor_j, "ITP");
    double res = expo_pos*(1- pow(expo_min, 
                                  itp_ij.size()));
    // 2. Step: For all h in ITP of i and j check their ISP between j and h 
    
    std::unordered_set<int>::iterator itr;
    for (itr = itp_ij.begin(); itr != itp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(actor_j, *itr, "ITP")); 
      res += pow(expo_min, object.count_common_partners_nb(*itr,actor_i,  "ITP")); 
    }
    return(res);
  }else {  
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_local_ITP", ::xyz_stat_gwesp_local_ITP, "gwesp_local_ITP",0.0);

auto xyz_stat_gwesp_local_ISP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all ISP of i and j 
    double res = expo_pos*(1- pow(expo_min, 
                                  object.count_common_partners_nb(actor_i, actor_j, "ISP")));
    // 2. Step: For all h in OSP of i and j check their ISP between j and h 
    std::unordered_set<int> osp_ij = object.get_common_partners_nb(actor_i, actor_j, "OSP");
    std::unordered_set<int>::iterator itr;
    for (itr = osp_ij.begin(); itr != osp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(actor_j, *itr, "ISP")); 
    }
    // 3. Step: For all h in OTP of i and j check their ISP between h and j
    std::unordered_set<int> otp_ij = object.get_common_partners_nb(actor_i, actor_j, "OTP");
    for (itr = otp_ij.begin(); itr != otp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(*itr,actor_j,  "ISP")); 
    }
    return(res);
  }else {
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_local_ISP", ::xyz_stat_gwesp_local_ISP, "gwesp_local_ISP",0.0);

auto xyz_stat_gwesp_local_OTP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all OTP of i and j 
    
    double res = expo_pos*(1- pow(expo_min, 
                                  object.count_common_partners_nb(actor_i, actor_j, "OTP")));
    // 2. Step: 
    std::unordered_set<int> osp_ij = object.get_common_partners_nb(actor_i, actor_j, "OSP");
    std::unordered_set<int>::iterator itr;
    for (itr = osp_ij.begin(); itr != osp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(actor_i, *itr, "OTP")); 
    }
    // 3. Step:
    std::unordered_set<int> isp_ij = object.get_common_partners_nb(actor_i, actor_j, "ISP");
    for (itr = isp_ij.begin(); itr != isp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(*itr,actor_j,  "OTP")); 
    }
    return(res);
  }else { 
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_local_OTP", ::xyz_stat_gwesp_local_OTP, "gwesp_local_OTP",0.0);

auto xyz_stat_gwesp_local_OSP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all OSP of i and j 
    double res = expo_pos*(1- pow(expo_min, 
                                  object.count_common_partners_nb(actor_i, actor_j, "OSP")));
    // 2. Step: 
    std::unordered_set<int>::iterator itr;
    std::unordered_set<int> otp_ij = object.get_common_partners_nb(actor_i, actor_j, "OTP");
    for (itr = otp_ij.begin(); itr != otp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(actor_i, *itr, "OSP")); 
    }
    // 3. Step:
    std::unordered_set<int> isp_ij = object.get_common_partners_nb(actor_i, actor_j, "ISP");
    for (itr = isp_ij.begin(); itr != isp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners_nb(*itr,actor_i,  "OSP")); 
    }
    return(res);
  }else { 
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_local_OSP", ::xyz_stat_gwesp_local_OSP, "gwesp_local_OSP",0.0);

auto xyz_stat_gwesp_ITP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all ISP of i and j 
    std::unordered_set<int> itp_ij = object.get_common_partners(actor_i, actor_j, "ITP");
    double res = expo_pos*(1- pow(expo_min, 
                                  itp_ij.size()));
    // 2. Step: For all h in ITP of i and j check their ISP between j and h 
    
    std::unordered_set<int>::iterator itr;
    for (itr = itp_ij.begin(); itr != itp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(actor_j, *itr, "ITP")); 
      res += pow(expo_min, object.count_common_partners(*itr,actor_i,  "ITP")); 
    }
    return(res);
  }else { 
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_global_ITP", ::xyz_stat_gwesp_ITP, "gwesp_global_ITP",0.0);

auto xyz_stat_gwesp_ISP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all ISP of i and j 
    double res = expo_pos*(1- pow(expo_min, 
                                  object.count_common_partners(actor_i, actor_j, "ISP")));
    // 2. Step: For all h in OSP of i and j check their ISP between j and h 
    std::unordered_set<int> osp_ij = object.get_common_partners(actor_i, actor_j, "OSP");
    std::unordered_set<int>::iterator itr;
    for (itr = osp_ij.begin(); itr != osp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(actor_j, *itr, "ISP")); 
    }
    // 3. Step: For all h in OTP of i and j check their ISP between h and j
    std::unordered_set<int> otp_ij = object.get_common_partners(actor_i, actor_j, "OTP");
    for (itr = otp_ij.begin(); itr != otp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(*itr,actor_j,  "ISP")); 
    }
    return(res);
  }else {
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_global_ISP", ::xyz_stat_gwesp_ISP, "gwesp_global_ISP",0.0);

auto xyz_stat_gwesp_OTP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all OTP of i and j 
    
    double res = expo_pos*(1- pow(expo_min, 
                                  object.count_common_partners(actor_i, actor_j, "OTP")));
    // 2. Step: 
    std::unordered_set<int> osp_ij = object.get_common_partners(actor_i, actor_j, "OSP");
    std::unordered_set<int>::iterator itr;
    for (itr = osp_ij.begin(); itr != osp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(actor_i, *itr, "OTP")); 
    }
    // 3. Step:
    std::unordered_set<int> isp_ij = object.get_common_partners(actor_i, actor_j, "ISP");
    for (itr = isp_ij.begin(); itr != isp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(*itr,actor_j,  "OTP")); 
    }
    return(res);
  }else {  
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_global_OTP", ::xyz_stat_gwesp_OTP, "gwesp_global_OTP",0.0);

auto xyz_stat_gwesp_OSP= CHANGESTAT{
  if(mode == "z"){
    double expo_min = (1-exp(-data.at(0,0)));  
    double expo_pos = exp(data.at(0,0));
    // 1. Step: For all OSP of i and j 
    double res = expo_pos*(1- pow(expo_min, 
                                  object.count_common_partners(actor_i, actor_j, "OSP")));
    // 2. Step: 
    std::unordered_set<int>::iterator itr;
    std::unordered_set<int> otp_ij = object.get_common_partners(actor_i, actor_j, "OTP");
    for (itr = otp_ij.begin(); itr != otp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(actor_i, *itr, "OSP")); 
    }
    // 3. Step:
    std::unordered_set<int> isp_ij = object.get_common_partners(actor_i, actor_j, "ISP");
    for (itr = isp_ij.begin(); itr != isp_ij.end(); itr++) {
      res += pow(expo_min, object.count_common_partners(*itr,actor_i,  "OSP")); 
    }
    return(res);
  }else { 
    return(0);
  }
}; 
EFFECT_REGISTER("gwesp_global_OSP", ::xyz_stat_gwesp_OSP, "gwesp_global_OSP",0.0);


// arma::vec stat_two_star_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
//                             const std::unordered_map< int, std::set<int>> &edges_neg,
//                             int &n_actors,
//                             int &tmp_random_i,
//                             int &tmp_random_j, 
//                             arma::mat &data, 
//                             int &type){
//   arma::vec res(3);
//   int degree_i = edges_pos.at(tmp_random_i).size();
//   int degree_j = edges_pos.at(tmp_random_j).size();
//   if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
//     res.at(0) = degree_i + degree_j-2;
//     res.at(1) = 0;
//     res.at(2) = -res.at(0);
//   } else {
//     res.at(0) = degree_i + degree_j;
//     res.at(1) = 0;
//     res.at(2) = -res.at(0);
//   }
//   return(res);
// }

// arma::vec stat_gwdsp(const std::unordered_map< int, std::set<int>> &edges_pos,
//                          const  std::unordered_map< int, std::set<int>> &edges_neg,
//                          int &n_actors,
//                          int &tmp_random_i,
//                          int &tmp_random_j, 
//                          arma::mat &data, 
//                          int &type){
//   arma::vec res(3);
//   res.fill(0);
//   double expo_min = (1-exp(-data.at(0,0)));  
//   // double expo_pos = exp(data.at(0,0));  
//   
//   std::set<int> intersection, tmp,intersection_ij;
//   // Step 1: Get all connections of tmp_random_i
//   std::set<int> connections_of_i = edges_pos.at(tmp_random_i);
//   // Step 2: Get all connections of tmp_random_j
//   std::set<int> connections_of_j = edges_pos.at(tmp_random_j);
//   connections_of_j.erase(tmp_random_i);
//   connections_of_i.erase(tmp_random_j);
//   // Step 3: Go through all connections of j and count its two-paths with i (call them h)
//   std::set<int>::iterator itr;
//   
//   for (itr = connections_of_j.begin(); itr != connections_of_j.end(); itr++) {
//     // For each connection of j we go through the positive edges 
//     // and count how many of them are intersecting with connections of i 
//     tmp = edges_pos.at(*itr);
//     // Get the size of the intersection of tmp and connections_of_i
//     intersection.clear();    
//     std::set_intersection(std::begin(tmp), std::end(tmp),
//                           std::begin(connections_of_i), std::end(connections_of_i),
//                           std::inserter(intersection, std::begin(intersection)));
//     res.at(0) +=  pow(expo_min, intersection.size());
//   }
//   // Step 4: Go through all connections of i and count its two-paths with j
//   for (itr = connections_of_i.begin(); itr != connections_of_i.end(); itr++) {
//     // For each connection of j we go through the positive edges 
//     tmp = edges_pos.at(*itr);
//     // Get the size of the intersection of tmp and connections_of_i
//     intersection.clear();    
//     std::set_intersection(std::begin(tmp), std::end(tmp),
//                           std::begin(connections_of_j), std::end(connections_of_j),
//                           std::inserter(intersection, std::begin(intersection)));
//     res.at(0) +=  pow(expo_min, intersection.size());
//   }
//   res.at(1) = 0;
//   res.at(2) = -res.at(0);
//   
//   return(res);
// }
// arma::vec stat_gwese_pos(const std::unordered_map< int, std::set<int>> &edges_pos,
//                          const std::unordered_map< int, std::set<int>> &edges_neg,
//                          int &n_actors,
//                          int &tmp_random_i,
//                          int &tmp_random_j, 
//                          arma::mat &data, 
//                          int &type){
//   arma::vec res(3);
//   res.fill(0);
//   double expo_min = (1-exp(-data.at(0,0)));  
//   double expo_pos = exp(data.at(0,0));  
//   
//   std::set<int>::iterator itr;
//   std::set<int> intersection_ij,intersection, tmp;
//   std::set<int> connections_of_i_neg = edges_neg.at(tmp_random_i);
//   std::set<int> connections_of_j_neg = edges_neg.at(tmp_random_j);
//   std::set<int> connections_of_i_pos = edges_pos.at(tmp_random_i);
//   std::set<int> connections_of_j_pos = edges_pos.at(tmp_random_j);
// 
//   // Step 3: Get the intersection of connections to i and j -> case a)
//   // Rcout << "Positive connections of i" << std::endl;
//   std::set_intersection(std::begin(connections_of_i_neg), std::end(connections_of_i_neg),
//                         std::begin(connections_of_j_neg), std::end(connections_of_j_neg),
//                         std::inserter(intersection_ij, std::begin(intersection_ij)));
//   res.at(0) += expo_pos*(1- pow(expo_min, intersection_ij.size()));
//   
//   // If tmp_i and tmp_j are already connected, we need to delete them from the tmp graph 
//   connections_of_i_neg.erase(tmp_random_j);
//   connections_of_j_neg.erase(tmp_random_i);
//   // Check for each positive connection of i that has a negative connection to j how many negative common partner it has with i
//   // How many common enemies do we have?  
//   for (itr = connections_of_i_pos.begin(); itr != connections_of_i_pos.end(); itr++) {
//     tmp = edges_neg.at(*itr);
//     if(tmp.count(tmp_random_j)){
//       intersection.clear();    
//       std::set_intersection(std::begin(tmp), std::end(tmp),
//                             std::begin(connections_of_i_neg), std::end(connections_of_i_neg),
//                             std::inserter(intersection, std::begin(intersection)));
//       res.at(1) +=  pow(expo_min, intersection.size());
//       // Rcout << intersection.size() << std::endl;
//       
//     }
//   }
//   // Check for each positive connection of j that has a negative connection to i how many negative common partner it has with j
//   // How many common enemies do we have?  
//   for (itr = connections_of_j_pos.begin(); itr != connections_of_j_pos.end(); itr++) {
//     tmp = edges_neg.at(*itr);
//     if(tmp.count(tmp_random_i)){
//       intersection.clear();    
//       std::set_intersection(std::begin(tmp), std::end(tmp),
//                             std::begin(connections_of_j_neg), std::end(connections_of_j_neg),
//                             std::inserter(intersection, std::begin(intersection)));
//       res.at(1) +=  pow(expo_min, intersection.size());
//       }
//   }
//   res.at(2) = res.at(1) - res.at(0);
//   return(res);
// }
// arma::vec stat_gwdegree(const std::unordered_map< int, std::set<int>> &edges_pos,
//                         const std::unordered_map< int, std::set<int>> &edges_neg,
//                         int &n_actors,
//                         int &tmp_random_i,
//                         int &tmp_random_j, 
//                         arma::mat &data, 
//                         int &type){
//   arma::vec res(3);
//   int degree_i = edges_pos.at(tmp_random_i).size() + edges_neg.at(tmp_random_i).size();
//   int degree_j = edges_pos.at(tmp_random_j).size() + edges_neg.at(tmp_random_j).size();
//   double expo_min = (1-exp(-data.at(0,0)));  
//   // double expo_pos = exp(data.at(0,0)); 
//   // If the edge is already there, we need to substract one of the degrees
//   if(edges_pos.at(tmp_random_i).count(tmp_random_j)){
//     degree_i -= 1;
//     degree_j -= 1;
//   }
//   if(edges_neg.at(tmp_random_i).count(tmp_random_j)){
//     degree_i -= 1;
//     degree_j -= 1;
//   } 
// 
//   res.at(0) = pow(expo_min, degree_i) + pow(expo_min, degree_j);
//   res.at(1) = res.at(0);
//   res.at(2) = 0;
//   return(res);
// }
