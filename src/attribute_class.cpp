#include "iglm/attribute_class.h"

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_NEWARP
#define ARMA_DONT_USE_ARPACK

Attribute::Attribute (int a, std::string type_, double scale_) {
  n_actors = a;
  scale = scale_;
 
  
  arma::vec tmp(a);
  tmp.fill(0);
  attribute = tmp;
  if(type_ == "binomial" || type_ == "poisson" || type_ == "normal"){
    type = type_;   
  } else {
    Rcpp::Rcout << "Invalid type, we assume that it is binomial (only binomial, poisson, and normal are implemented.)\n";
    type = "binomial";
  }
  // std::cout << "Size";
  // std::cout << attribute.size();
}
Attribute::Attribute (int a, arma::vec attribute_tmp, std::string type_, double scale_) {
  n_actors = a;
  scale = scale_;
  attribute = attribute_tmp;
  if(type_ == "binomial" || type_ == "poisson" || type_ == "normal"){
    type = type_;   
  } else {
    Rcpp::Rcout << "Invalid type, we assume that it is binomial (only binomial, poisson, and normal are implemented.)\n";
    type = "binomial";
  }
}

// function object that acts as a specific function
// we overload the operator "IsMoreThan" for that reason
struct IsMoreThan {
  IsMoreThan(int i) : i_{i} {}
  bool operator()(int i) { return i > i_; }
  int i_;
};


double Attribute::get_scale() const{
  return(scale);
}

// This is just a function to check whether the attribute is in line with expectations (all entries should be below the number of actors)
bool Attribute::check () const {
  // int number_attribute = attribute.size();
  // How many entries in the set are larger than n_actors?
  int mycount = std::count_if(attribute.begin(),
                              attribute.end(),
                              IsMoreThan(n_actors));
  // If this count is larger than 0 return false else true
  return(mycount==0);
}

// How many attributes are set to 1?
int Attribute::get_sum()const{
  return(attribute.size());
}

double Attribute::get_val(int from) const {
  if(from > n_actors) {
    // std::cout << "Error: trying to access attribute of non-existing actor.\n";
    return(false);
  }
  return(attribute.at(from-1)/scale);
}

void Attribute::print() {
  attribute.print();
}


// void Attribute::initialize() {
//   std::unordered_set<int> tmp;
//
//   attribute = tmp;
// }

void Attribute::set_attr_0(int from) {
  attribute(from-1) = 0;
}

void Attribute::set_attr_value(int from, double val) {
  attribute(from-1) = val;
}

// void Attribute::set_attr(std::unordered_set<int> attribute_tmp) {
//   attribute = attribute_tmp;
// }


void Attribute::set_attr_from_armavec(arma::vec attribute_tmp) {
  attribute = attribute_tmp;
}


void Attribute::set_attr_1(int from) {
  attribute(from-1) = 1;
}
// // [[Rcpp::export]]
// int tmp_function() {
// int a = 3;
// Attribute attribute(a);
// attribute.initialize();
// attribute.print();
//   attribute.set_attr_1(a-1);
//   attribute.print();
//   // Rcout << attribute.vector(a-1) << std::endl;
//   return (1);
// }
