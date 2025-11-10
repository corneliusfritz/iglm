// Defines a header file containing function signatures for functions in src/

// Protect signatures using an inclusion guard.
#ifndef network_class_H
#define network_class_H
#define DARMA_USE_CURRENT
#include <RcppArmadillo.h>
#include <set>
#include <unordered_map>

class Network {
public:
  // Members
  bool directed;
  // Note that adj_list refers to the outgoing and adj_list_in to the ingoing ties
  // for undirected networks adj_list_in is not initialized and never used 
  std::unordered_map< int, std::unordered_set<int>> adj_list;
  std::unordered_map< int, std::unordered_set<int>> adj_list_in;
  // Constructors
  Network(int, bool); // Empty Constructor
  Network(int, bool, arma::mat mat); // Matrix Constructor
  // Member functions
  // void initialize();
  arma::uvec get_common_partners(unsigned int,unsigned int, std::string) const;
  void add_edge(int, int);
  void delete_edge(int, int);
  double get_val(int, int) const;
  void change_edge(int, int);
  void add_edges_from_mat(arma::mat mat);
  double number_edges() const;
private:
  // Private members
  int n_actors;
};
#endif

