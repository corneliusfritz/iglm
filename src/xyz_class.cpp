#include "iglm/xyz_class.h"
// #include "xyz_sampling.h" 

// Constructors
XYZ_class::XYZ_class(int n_actor_, bool directed_, std::string type_x_,std::string type_y_, 
                     double scale_x_, double scale_y_): XZ_class(n_actor_,directed_, type_x_, scale_x_), y_attribute(n_actor_, type_y_, scale_y_){

} 
XYZ_class::XYZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_x_,std::string type_y_, double scale_x_, double scale_y_): 
  XZ_class(n_actor_,directed_, neighborhood_, overlap_, type_x_, scale_x_), y_attribute(n_actor_, type_y_, scale_y_){
  n_actor = n_actor_;
}

XYZ_class::XYZ_class(int n_actor_, bool directed_,std::unordered_map< int, std::unordered_set<int>> neighborhood_,
                     std::unordered_map< int, std::unordered_set<int>> overlap_, 
                     arma::mat overlap_mat_, 
                     std::string type_x_,std::string type_y_,
                     double scale_x_, double scale_y_): XZ_class(n_actor_,directed_, neighborhood_, overlap_, overlap_mat_, type_x_, scale_x_),  y_attribute(n_actor_, type_y_, scale_y_){
  n_actor = n_actor_;
}


XYZ_class::XYZ_class(int n_actor_, bool directed_,  arma::vec x_attribute_, arma::vec y_attribute_, arma::mat z_network_,arma::mat neighborhood_, 
                     arma::mat overlap_, std::string type_x_,std::string type_y_, double scale_x_, double scale_y_):
 XZ_class(n_actor_,directed_, z_network_,x_attribute_, neighborhood_,overlap_, type_x_, scale_x_),  y_attribute(n_actor_,y_attribute_, type_y_, scale_y_){
}


// XYZ_class::XYZ_class(int b, bool a): y_attribute(b), XZ_class(b,a){
// }

void XYZ_class::print() {
  Rcout << "X Attribute" << std::endl;
  x_attribute.print();
  Rcout << "Y Attribute" << std::endl;
  y_attribute.print();
  Rcout << "Z Network" << std::endl;
  Rcout << map_to_mat(z_network.adj_list, n_actor) << std::endl;
  Rcout << "Neighborhood Matrix" << std::endl;
  Rcout << map_to_mat(neighborhood, n_actor) << std::endl;
}

void XYZ_class::set_info(arma::vec x_attribute_,arma::vec y_attribute_,
                         std::unordered_map< int, std::unordered_set<int>> z_network_) {
  x_attribute.attribute = x_attribute_;
  y_attribute.attribute = y_attribute_;
  z_network.adj_list = z_network_;
}

void XYZ_class::set_info_arma(arma::vec x_attribute_, arma::vec y_attribute_, arma::mat z_network_) {
  x_attribute.attribute = x_attribute_;
  y_attribute.attribute = y_attribute_;
  mat_to_map(z_network_,n_actor, z_network.directed, z_network.adj_list,z_network.adj_list_in);
}


void XYZ_class::copy_from(XYZ_class obj) {
  x_attribute = obj.x_attribute;
  y_attribute = obj.y_attribute;
  z_network = obj.z_network;
  neighborhood = obj.neighborhood;
  n_actor = obj.n_actor;
}




//
// void XYZ_class::initialize(arma::mat network_tmp, arma::vec x_attribute_tmp,arma::vec y_attribute_tmp, arma::mat neighborhood_tmp) {
//   x_attribute.initialize();
//   y_attribute.initialize();
//   z_network.initialize();
//   neighborhood_initialize();
//   z_network.add_edges_from_mat(network_tmp);
//   set_neighborhood_from_mat(neighborhood_tmp);
//   x_attribute.set_attr_from_armavec(x_attribute_tmp);
//   y_attribute.set_attr_from_armavec(y_attribute_tmp);
// }
//
// void XYZ_class::initialize_from_zero(arma::mat neighborhood_tmp) {
//   x_attribute.initialize();
//   y_attribute.initialize();
//   z_network.initialize();
//
//   neighborhood_initialize();
//   set_neighborhood_from_mat(neighborhood_tmp);
// }

