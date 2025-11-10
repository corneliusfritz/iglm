#include "iglm/xz_class.h"
// Constructors
XZ_class::XZ_class(int n_actor_, bool directed_, std::string type_, double scale_): x_attribute(n_actor_, type_, scale_), z_network(n_actor_,directed_){
  n_actor = n_actor_;
  for (int i = 1; i <= n_actor; i++){
    // Rcout << i  << std::endl;
    neighborhood[i] = std::unordered_set<int>();
    overlap[i] = std::unordered_set<int>();
  }
  overlap_mat = arma::zeros<arma::mat>(2, 0);
}

XZ_class::XZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_): x_attribute(n_actor_, type_, scale_), z_network(n_actor_,directed_){
  n_actor = n_actor_;
  mat_to_map_neighborhood(neighborhood_,overlap_, n_actor_, directed_,
                          neighborhood,
                          overlap);
  overlap_mat = overlap_;
}


XZ_class::XZ_class(int n_actor_, bool directed_, std::unordered_map< int, std::unordered_set<int>> neighborhood_,
                   std::unordered_map< int, std::unordered_set<int>> overlap_,
                   arma::mat overlap_mat_, std::string type_, double scale_): x_attribute(n_actor_, type_, scale_), z_network(n_actor_,directed_){

  neighborhood = neighborhood_;
  overlap = overlap_;
  overlap_mat = overlap_mat_;
  n_actor = n_actor_;
}


XZ_class::XZ_class(int n_actor_, bool directed_, arma::mat z_network_, arma::vec x_attribute_, 
                   arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
  x_attribute(n_actor_, x_attribute_, type_, scale_), z_network(n_actor_,directed_, z_network_){
  n_actor = n_actor_;
  mat_to_map_neighborhood(neighborhood_,overlap_, n_actor, directed_,
                          neighborhood,
                          overlap);
  overlap_mat = overlap_;
  
  // Rcout <<  map_to_mat(overlap, n_actor) << std::endl;
  
}


bool XZ_class::get_val_overlap(int from, int to )const {
  if(overlap.at(from).count(to)){
    return(true);
  } else if(overlap.at(to).count(from)) {
    return(true);
  } else {
    return(false);
  }
}


bool XZ_class::get_val_neighborhood(int from, int to ) const{
  if(neighborhood.at(from).count(to)){
    return(true);
  } else {
    return(false);
  }
}

// This function checks if the neighborhood is full or not
bool XZ_class::check_if_full_neighborhood() {
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
void XZ_class::print() {
  Rcout << "Network" << std::endl;
  Rcout << map_to_mat(z_network.adj_list, n_actor) << std::endl;
  Rcout << "Attribute" << std::endl;
  x_attribute.print();
  Rcout << "Neighborhood Matrix" << std::endl;
  Rcout << map_to_mat(neighborhood, n_actor) << std::endl;
  // print_set(neighborhood,n_ac);
}

void XZ_class::copy_from(XZ_class obj) {
  z_network = obj.z_network;
  x_attribute = obj.x_attribute;
  neighborhood = obj.neighborhood;
  n_actor = obj.n_actor;
}


void XZ_class::set_neighborhood_from_mat(arma::mat mat) {
  arma::rowvec tmp_row;
  for (int i = 1; i <= n_actor; i++){
    tmp_row = mat.row(i-1);
    arma::uvec ids = find(tmp_row == 1) + 1; // Find indices
    // Rcout << ids << std::endl;
    neighborhood.at(i)= std::unordered_set<int>(ids.begin(),ids.end());
  }
}

void XZ_class::neighborhood_initialize() {
  for (int i = 1; i <= n_actor; i++){
    // Rcout << i  << std::endl;
    neighborhood[i] = std::unordered_set<int>();
  }
}

// void XZ_class::neighborhood_randomization(int seed, int K) {
//   arma::vec memberships =  simulate_numbers(1, K,n_actor,seed);
//   // Rcout << memberships  << std::endl;
//   arma::uvec tmp_member;
//   for(int i = 1; i<= n_actor; i ++){
//     // Rcout << i  << std::endl;
//     tmp_member = find(memberships == memberships.at(i-1)) +1;
//     // Rcout << tmp_member  << std::endl;
//     // std::vector<int> tmp(tmp_member.begin(), tmp_member.end());
//     // tmp = arma::conv_to<std::vector<int>>::from(tmp_member);
//     // std::unordered_set<int> s(tmp_member.begin(), tmp_member.end());
//     // print_set(s);
//     std::vector<int> tmp = arma::conv_to< std::vector<int> >::from(tmp_member);
//     std::unordered_set<int> s(tmp.begin(), tmp.end());
//     neighborhood.at(i)=  s;
//     // print_set(s);
//     // neighborhood.at(i) = std::unordered_set<int> s(tmp_member.begin(), tmp_member.end());
//   }
// }

void XZ_class::assign_neighborhood(std::unordered_map< int, std::unordered_set<int>> new_neighborhood) {
  neighborhood = new_neighborhood;
}
void XZ_class::set_info(arma::vec x_attribute_, std::unordered_map< int, std::unordered_set<int>> z_network_) {
  x_attribute.attribute = x_attribute_;
  z_network.adj_list = z_network_;
}

void XZ_class::set_info_arma(arma::vec x_attribute_, arma::mat z_network_) {
  x_attribute.attribute = x_attribute_;
  mat_to_map(z_network_,n_actor, z_network.directed, z_network.adj_list,z_network.adj_list_in);
}

void XZ_class::change_neighborhood(int actor, std::unordered_set<int> new_neighborhood) {
  neighborhood.at(actor -1 ) = new_neighborhood;
}
// void XZ_class::initialize(arma::mat network_tmp, arma::vec attribute_tmp, arma::mat neighborhood_tmp) {
//   z_network.initialize();
//   x_attribute.initialize();
//   neighborhood_initialize();
//   z_network.add_edges_from_mat(network_tmp);
//   set_neighborhood_from_mat(neighborhood_tmp);
//   x_attribute.set_attr_from_armavec(attribute_tmp);
// }
//
// void XZ_class::initialize_from_zero(arma::mat neighborhood_tmp) {
//   z_network.initialize();
//   x_attribute.initialize();
//   neighborhood_initialize();
//   set_neighborhood_from_mat(neighborhood_tmp);
// }
