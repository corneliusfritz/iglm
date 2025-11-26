
#ifndef helper_H
#define helper_H
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <random>
#include <set>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include "iglm/extension_api.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

template<typename T>
void print_vector(std::vector<T> tmp){
  // std::vector<T>::iterator itr;
  for (std::vector<int>::iterator it = tmp.begin() ; it != tmp.end(); ++it) {
    Rcout << *it << std::endl;
  }
  // for (itr = tmp.begin(); itr != tmp.end(); itr++) {
  //
  // }
}



template<typename T>
std::vector<int> findItems(std::vector<T> const &v, int target) {
  std::vector<int> indices;
  auto it = v.begin();
  while ((it = std::find_if(it, v.end(), [&] (T const &e) { return e == target; }))
           != v.end())
  {
    indices.push_back(std::distance(v.begin(), it));
    it++;
  }
  return indices;
}

inline void print_set(std::unordered_set<int> tmp){
  std::unordered_set<int>::iterator itr;
  for (itr = tmp.begin(); itr != tmp.end(); itr++) {
    Rcout << *itr << std::endl;
  }
}

inline void set_seed(int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

// This is a function to transform a map to one matrix
inline arma::mat map_to_mat(std::unordered_map< int, std::unordered_set<int>> &adj_list,
                     int & n_actors) {
  arma::mat res(n_actors,n_actors);
  res.fill(0);
  arma::mat tmp_row(n_actors,1);
  NumericVector tmp_vec(n_actors);
  for (int i = 1; i <= n_actors; i++){
    // Set the vector to 0
    tmp_row.fill(0);
    tmp_vec = adj_list[i];
    // Go through each row and write the matrix
    for (int j = 0; j < tmp_vec.length() ; j++){
      tmp_row.at(tmp_vec.at(j)-1) = 1;
    }
    res.col(i-1) = tmp_row.as_col();
  }
  return(res);
}

// This is a function to transform a matrix to a map
inline void mat_to_map(arma::mat mat, int n_actors, bool directed,
                       std::unordered_map< int, std::unordered_set<int>> &adj_list,
                       std::unordered_map< int, std::unordered_set<int>> &adj_list_in) {
  if (mat.is_empty() || mat.n_elem == 0) {
    return; 
  }
  for (int i = 1; i <= n_actors; i++){ 
    adj_list[i] = std::unordered_set<int>();
  } 
  if(directed){
    for (int i = 1; i <= n_actors; i++){
      adj_list_in[i] = std::unordered_set<int>();
    }  
  }
  //  This checks whether an edge list or adjacency matrix is provided
  if(mat.n_cols==2){
    if(directed){
      arma::vec tmp_row1 = mat.col(0);
      arma::vec tmp_row2 = mat.col(1);
      for (arma::uword i = 1; i <= n_actors; i++){
        arma::vec ids1 = tmp_row1.elem(find(mat.col(1) == i)); // Find indices
        arma::vec ids2 = tmp_row2.elem(find(mat.col(0) == i)); // Find indices
        adj_list.at(i)= std::unordered_set<int>(ids2.begin(), ids2.end());
        adj_list_in.at(i)= std::unordered_set<int>(ids1.begin(), ids1.end());
      }
    } else { 
      arma::vec tmp_row1 = mat.col(0);
      arma::vec tmp_row2 = mat.col(1);

      // double u = 0.0;
      for (arma::uword i = 1; i <= n_actors; i++){
        std::unordered_set<int> part1,part2, res;
        arma::vec ids1 = tmp_row1.elem(find(mat.col(1) == i)); // Find indices from
        arma::vec ids2 = tmp_row2.elem(find(mat.col(0) == i)); // Find indices to
        part1= std::unordered_set<int>(ids1.begin(), ids1.end());
        part2= std::unordered_set<int>(ids2.begin(), ids2.end());
        std::set_union(std::begin(part1), std::end(part1),
                              std::begin(part2), std::end(part2),
                              std::inserter(res, std::begin(res)));
        
        adj_list.at(i)= std::unordered_set<int>(res.begin(), res.end());
        // u += res.size();
      }
      // Rcout << "The undirected network has " << u/2 << " edges." << std::endl;
    } 
  } else {
    arma::rowvec tmp_row;
    for (arma::uword i = 1; i <= n_actors; i++){
      tmp_row = mat.row(i-1);
      arma::uvec ids = find(tmp_row == 1) + 1; // Find indices
      // Rcout << ids << std::endl;
      adj_list.at(i)= std::unordered_set<int>(ids.begin(),ids.end());
    } 
    if(directed){
      // Rcout << "Starting with the in degree network" << std::endl;
      for (arma::uword i = 1; i <= n_actors; i++){
        tmp_row = mat.col(i-1).as_row();
        // Rcout << tmp_col << std::endl;
        arma::uvec ids = find(tmp_row == 1) + 1; // Find indices
        // Rcout << ids << std::endl;
        adj_list_in.at(i)= std::unordered_set<int>(ids.begin(),ids.end());
      } 
    } 
  }
}

inline void mat_to_map_neighborhood(arma::mat neighborhood_,arma::mat overlap_, int n_actors, bool directed,
                       std::unordered_map< int, std::unordered_set<int>> &neighborhood,
                       std::unordered_map< int, std::unordered_set<int>> &overlap) {
  for (int i = 1; i <= n_actors; i++){ 
    neighborhood[i] = std::unordered_set<int>();
    overlap[i] = std::unordered_set<int>();
  } 

  if(neighborhood_.n_cols==2){
    arma::vec tmp_row1 = neighborhood_.col(0);
    arma::vec tmp_row2 = neighborhood_.col(1);
    arma::vec tmp_row1_ov = overlap_.col(0);
    arma::vec tmp_row2_ov = overlap_.col(1);
    for (int i = 1; i <= n_actors; i++){
      std::unordered_set<int> part1,part2, res;
      arma::vec ids1 = tmp_row1.elem(find(neighborhood_.col(1) == i)); // Find indices
      arma::vec ids2 = tmp_row2.elem(find(neighborhood_.col(0) == i)); // Find indices
      part1= std::unordered_set<int>(ids1.begin(), ids1.end());
      part2= std::unordered_set<int>(ids2.begin(), ids2.end());
      std::set_union(std::begin(part1), std::end(part1),
                     std::begin(part2), std::end(part2),
                     std::inserter(res, std::begin(res)));
      neighborhood.at(i)= std::unordered_set<int>(res.begin(), res.end());
      
      std::unordered_set<int> res_nb;
      ids1 = tmp_row1_ov.elem(find(overlap_.col(1) == i)); // Find indices
      ids2 = tmp_row2_ov.elem(find(overlap_.col(0) == i)); // Find indices
      part1= std::unordered_set<int>(ids1.begin(), ids1.end());
      part2= std::unordered_set<int>(ids2.begin(), ids2.end());
      std::set_union(std::begin(part1), std::end(part1),
                     std::begin(part2), std::end(part2),
                     std::inserter(res_nb, std::begin(res_nb)));
      overlap.at(i)= std::unordered_set<int>(res_nb.begin(), res_nb.end());
    }
  } else {
    arma::rowvec tmp_row;
    arma::uvec ids;
    for (int i = 1; i <= n_actors; i++){
      tmp_row = neighborhood_.row(i-1);
      ids = find(tmp_row == 1) + 1; // Find indices
      neighborhood.at(i)= std::unordered_set<int>(ids.begin(),ids.end());
      
      tmp_row = overlap_.row(i-1);
      ids = find(tmp_row == 1) + 1; // Find indices
      overlap.at(i)= std::unordered_set<int>(ids.begin(),ids.end());
    } 
  } 
}

// This is a function to transform an arma::vec to a std::vector<int>
inline std::unordered_set< int> armavec_to_set(arma::vec vec, int type) {
  std::unordered_set< int> res;
  arma::uvec ids = find(vec == type) + 1; // Find indices
  for (unsigned int i = 0; i < ids.size(); i++){
    res.insert(ids.at(i));
  }
  // print_vector(res);
  return(res);
}


// This is a function to transform an arma::vec to a std::vector<int>
inline std::vector< int> armavec_to_vector(arma::vec vec) {
  int n_actors = vec.size();
  std::vector< int> res;

  for (int i = 0; i < n_actors; i++){
    res.push_back(vec.at(i));
  }
  // print_vector(res);
  return(res);
}

inline std::unordered_set<int> get_set_difference(
    const std::unordered_set<int>& set1,
    const std::unordered_set<int>& set2)
{
  std::unordered_set<int> difference;
  difference.reserve(set1.size()); 
   
  for (const int& element : set1) {
    if (!set2.count(element)) {
      difference.insert(element);
    }
  } 
  
  return difference;
}

inline std::unordered_set<int> get_intersection(
    const std::unordered_set<int>& set1,
    const std::unordered_set<int>& set2)
{
  const auto& smaller_set = (set1.size() < set2.size()) ? set1 : set2;
  const auto& larger_set  = (set1.size() < set2.size()) ? set2 : set1;
  
  // Allocate the result set
  std::unordered_set<int> intersection;
  intersection.reserve(smaller_set.size());
  
  for (const int& element : smaller_set) {
    if (larger_set.count(element)) {
      intersection.insert(element);
    }
  }
  
  return intersection;
}


inline size_t count_difference(
    const std::unordered_set<int>& set1,
    const std::unordered_set<int>& set2)
{
  size_t count = 0;
  for (const int& element : set1) {
    if (!set2.count(element)) {
      count++;
    }
  } 
  
  return count;
} 

inline size_t count_intersection(
    const std::unordered_set<int>& set1,
    const std::unordered_set<int>& set2)
{ 
  const auto& smaller_set = (set1.size() < set2.size()) ? set1 : set2;
  const auto& larger_set  = (set1.size() < set2.size()) ? set2 : set1;
  size_t count = 0;
  // Allocate the result set
  
  for (const int& element : smaller_set) {
    if (larger_set.count(element)) {
      count++;
    }
  }
  
  return count;
}


#endif

