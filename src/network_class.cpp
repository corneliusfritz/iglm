// #include <Rcpp.h>
#include "iglm/network_class.h"
#include "iglm/helper_functions.h"
using namespace Rcpp;
Network::Network (int n_actors_, bool directed_) {
  n_actors = n_actors_;
  directed = directed_;
  for (int i = 1; i <= n_actors; i++){
    adj_list[i] = std::unordered_set<int>();
    adj_list_in[i] = std::unordered_set<int>();
  }
}

Network::Network (int n_actors_, bool directed_, arma::mat mat) {
  n_actors = n_actors_;
  directed = directed_;
  mat_to_map(mat,n_actors, directed, adj_list,adj_list_in);
  // Rcout << "Finished" << std::endl;
  
}
// This is a member function to change the state of a particular pai in the network
// (either from 1->0 or 0->1)
void Network::change_edge(int from, int to) {
  if(directed){
    if(adj_list.at(from).count(to)){
      adj_list.at(from).erase(to);
      adj_list_in.at(to).erase(from);
    } else {
      adj_list.at(from).insert(to);
      adj_list_in.at(to).insert(from);
      
    }
  } else {
    if(adj_list.at(from).count(to)){
      adj_list.at(from).erase(to);
      adj_list.at(to).erase(from);
    } else {
      adj_list.at(from).insert(to);
      adj_list.at(to).insert(from);
    }
  }
}


arma::uvec Network::get_common_partners(unsigned int from,unsigned int to, std::string type = "OSP")const {
  arma::uvec res;
  if(type == "OTP"){
    std::unordered_set<unsigned int> intersect_group;
    std::set_intersection(std::begin(adj_list.at(from)),
                          std::end(adj_list.at(from)),
                          std::begin(adj_list_in.at(to)),
                          std::end(adj_list_in.at(to)),
                          std::inserter(intersect_group, std::begin(intersect_group)));
    std::vector<int> output;
    std::copy(intersect_group.begin(), intersect_group.end(), std::back_inserter(output));
    res = arma::conv_to<arma::uvec>::from(output);
  } else  if(type == "ISP"){
    std::unordered_set<unsigned int> intersect_group;
    std::set_intersection(std::begin(adj_list_in.at(from)),
                          std::end(adj_list_in.at(from)),
                          std::begin(adj_list_in.at(to)),
                          std::end(adj_list_in.at(to)),
                          std::inserter(intersect_group, std::begin(intersect_group)));
    std::vector<int> output;
    std::copy(intersect_group.begin(), intersect_group.end(), std::back_inserter(output));
    res = arma::conv_to<arma::uvec>::from(output);
  }else  if(type == "OSP"){
    std::unordered_set<unsigned int> intersect_group;
    std::set_intersection(std::begin(adj_list.at(from)),
                          std::end(adj_list.at(from)),
                          std::begin(adj_list.at(to)),
                          std::end(adj_list.at(to)),
                          std::inserter(intersect_group, std::begin(intersect_group)));
    std::vector<int> output;
    std::copy(intersect_group.begin(), intersect_group.end(), std::back_inserter(output));
    res = arma::conv_to<arma::uvec>::from(output);
  }
  else  if(type == "ITP"){
    std::unordered_set<unsigned int> intersect_group;
    std::set_intersection(std::begin(adj_list_in.at(from)),
                          std::end(adj_list_in.at(from)),
                          std::begin(adj_list.at(to)),
                          std::end(adj_list.at(to)),
                          std::inserter(intersect_group, std::begin(intersect_group)));
    std::vector<int> output;
    std::copy(intersect_group.begin(), intersect_group.end(), std::back_inserter(output));
    res = arma::conv_to<arma::uvec>::from(output);
  }
  return(res);
}

double Network::number_edges() const{
  double count = 0.0;
  for(int i=1; i<=n_actors; i++){
    count += adj_list.at(i).size();
  }
  return(count);
}

double Network::get_val(int from, int to)const {
  if(adj_list.at(from).count(to)){
    return(1.0);
  } else {
    return(0.0);
  }
}

// This is a member function to add a particular matrix
void Network::add_edge(int from, int to) {
  if(directed){
    adj_list.at(from).insert(to);
    adj_list_in.at(to).insert(from);
  } else{
    adj_list.at(from).insert(to);
    adj_list.at(to).insert(from);
  }

}
// This is a member function to delete a particular matrix
void Network::delete_edge(int from, int to) {
  if(directed){
    adj_list.at(from).erase(to);
    adj_list_in.at(to).erase(from);
  } else{
    adj_list.at(from).erase(to);
    adj_list.at(to).erase(from);
  }
}
// This is a member function to initialize the matrix
// void Network::initialize() {
//   for (int i = 1; i <= n_actors; i++){
//     adj_list[i] = std::unordered_set<int>();
//   }
// }
// This is a member function to add edges to a network through an adjacency matrix
void Network::add_edges_from_mat(arma::mat mat) {
  mat_to_map(mat,n_actors, directed, adj_list,adj_list_in);
}

