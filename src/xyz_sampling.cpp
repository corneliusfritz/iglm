#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <random>
#include <set>
#include <unordered_map>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <limits>
#include "iglm/xyz_class.h"
#include "iglm/extension_api.hpp"

using xyz_ValidateFunction = double(*)(const XYZ_class &object,
                                    const int &actor_i,
                                    const int &actor_j,
                                    const arma::mat &data,
                                    const double &type,
                                    const std::string &mode, const bool &is_full_neighborhood);

//[[Rcpp::depends(RcppProgress)]]

void show_info() {
  auto& reg = iglm::Registry::instance();
  for (auto& name : reg.names()) {
    auto m = reg.info(name);
    Rcpp::Rcout << "Function: " << name
                << " | short: " << m.short_name
                << " | value: " << m.value << "\n";
  }
}

// [[Rcpp::export]]
void iglm_print_registered_functions() {
  auto meta = iglm::Registry::instance().all_meta();
  Rcpp::Rcout << "Registered functions:\n";
  if (meta.empty()) {
    Rcpp::Rcout << "  <none>\n";
    return;
  } 
  
  for (auto& kv : meta) {
    const std::string& name = kv.short_name;
    const iglm::FUN& m = kv;
    Rcpp::Rcout << "  " << name
                << "  (" << m.short_name << ", " << m.value << ")\n";
  }
}


std::vector<xyz_ValidateFunction> xyz_change_statistics_generate_new(std::vector<std::string> terms) {
  // Res
  std::vector<xyz_ValidateFunction> fns;
  fns.reserve(terms.size());
  // Get all the registered functions 
  auto& reg = iglm::Registry::instance();
  
  for (size_t i = 0; i < terms.size(); ++i) {
    const std::string& name = terms[i];
    if (!reg.has(name)) {
      throw std::invalid_argument("The statistic " + name + " does not exist");
    }    fns.push_back(reg.get(name));
  }
  return fns;
}


arma::vec xyz_eval_at_empty_network_new(std::vector<std::string> terms, const XYZ_class& object) {
  arma::vec res(terms.size());
  
  auto& reg = iglm::Registry::instance();
  for (size_t i = 0; i < terms.size(); ++i) {
    const std::string& name = terms[i];
    if (!reg.has(name)) {
      throw std::invalid_argument("The statistic " + name + " does not exist");
    }
    const auto meta = reg.info(name);
    if(meta.value == 1.0){
      res.at(i) = object.n_actor;
    } else {
      res.at(i) = meta.value;
    }
    
  }
  return res;
}

// Function to call all functions in the vector functions
inline void xyz_calculate_change_stats(arma::vec & change_stat,
                                            const int actor_i,
                                      const int actor_j,
                                      const XYZ_class &object,
                                      const std::vector<arma::mat> &data_list,
                                      const std::vector<double> &type_list,
                                      const std::string &mode,
                                      const bool &is_full_neighborhood,
                                      const std::vector<xyz_ValidateFunction> functions){

  // Rcout << functions.size() << std::endl;
  for(unsigned int a = 0; a < (functions.size()); a += 1 ) {
    change_stat.at(a) = functions.at(a)(object,actor_i,actor_j, data_list.at(a), type_list.at(a),mode, is_full_neighborhood);
  }
}



arma::vec xyz_count_global_statistic( XYZ_class &object,
                                      std::vector<arma::mat> &data_list,
                                      std::vector<double> &type_list,
                                      std::vector<xyz_ValidateFunction> functions, 
                                      std::string type_x, 
                                      std::string type_y, 
                                      double attr_x_scale, 
                                      double attr_y_scale) {
  // Generate empty network that we will fill as we go through all observed edges in the network
  XYZ_class alt_object(object.n_actor,object.z_network.directed, 
                       object.neighborhood, 
                       object.overlap,
                       object.overlap_mat,  
                       type_x, type_y,attr_x_scale, attr_y_scale);
  // XYZ_class alt_object(object.n_actor, object.z_network.directed, neighborhood);
  bool is_full_neighborhood = object.check_if_full_neighborhood();
  arma::vec res(functions.size());
  arma::vec change_stat(functions.size());
  // arma::vec tmp_row;
  std::unordered_set<int> tmp_js;
  std::string z = "z", x = "x", y = "y";
  // Go through all actors i and switch them incrementally from 0 to 1
  for (int i = 1; i <= object.n_actor; i++){
    tmp_js = object.z_network.adj_list.at(i);
    
    if(tmp_js.size()>0){
      std::unordered_set<int>::iterator it = tmp_js.begin();
      while (it != tmp_js.end()) {
        if(*it == i){ 
        } else if(!object.z_network.directed){
          if(*it>i){
            xyz_calculate_change_stats(change_stat, i,
                                                     *it,
                                                     alt_object,
                                                     data_list,
                                                     type_list,
                                                     z,
                                                     is_full_neighborhood,
                                                     functions);
            alt_object.z_network.add_edge(i,*it);
            res +=change_stat;
          }
        } else {
          // Rcout << "directed" << std::endl;
          xyz_calculate_change_stats(change_stat, i,
                                                   *it,
                                                   alt_object,
                                                   data_list,
                                                   type_list,
                                                   z,
                                                   is_full_neighborhood,
                                                   functions);
          alt_object.z_network.add_edge(i,*it);
          res +=change_stat;
        }
        it++;
        
      }
    }
  }
  // Rcout << "X Attr" << std::endl;
  for(int i = 1; i <= object.x_attribute.attribute.size(); i++){
    // Rcout << i << std::endl;
    xyz_calculate_change_stats(change_stat, i,
                                             i,
                                             alt_object,
                                             data_list,
                                             type_list,
                                             x,
                                             is_full_neighborhood,
                                             functions);
    // Rcout << change_stat << std::endl;
    // Rcout << "Here" << std::endl;
    // Rcout << object.x_attribute.attribute.size() << std::endl;
    // Rcout << object.x_attribute.attribute(i-1) << std::endl;
    
    alt_object.x_attribute.set_attr_value(i, object.x_attribute.attribute.at(i-1));
    res +=change_stat*object.x_attribute.attribute.at(i-1);
  }
  // Rcout << "Y Attr" << std::endl;
  
  for(int i = 1; i <= object.y_attribute.attribute.size(); i++){
    // Rcout << i << std::endl;
    xyz_calculate_change_stats(change_stat, i,
                                             i,
                                             alt_object,
                                             data_list,
                                             type_list,
                                             y,
                                             is_full_neighborhood,
                                             functions);
    // Rcout << change_stat << std::endl;
    // Rcout << "Here" << std::endl;
    alt_object.y_attribute.set_attr_value(i, object.y_attribute.attribute.at(i-1));
    res +=change_stat*object.y_attribute.attribute.at(i-1);
  }
  return(res);
}



// [[Rcpp::export]]
arma::vec xyz_count_global(arma::mat z_network,
                           arma::vec x_attribute,
                           arma::vec y_attribute,
                           arma::mat neighborhood,
                           arma::mat overlap,
                           bool directed,
                           std::vector<std::string> terms,
                           int n_actor,
                           std::vector<arma::mat> &data_list,
                           std::vector<double> &type_list, 
                           std::string type_x, 
                           std::string type_y, 
                           double attr_x_scale, 
                           double attr_y_scale) {
  // std::unordered_map< int, std::unordered_set<int>> edges;
  // // Convert the matrix to two unordered_map objects
  // edges = mat_to_map(network,1, n_actor);
  XYZ_class object(n_actor,directed, x_attribute,y_attribute,z_network, neighborhood, overlap, type_x, type_y,attr_x_scale, attr_y_scale);
  // object.initialize(z_network, x_attribute,y_attribute,neighborhood );
  // object.print();
  std::vector<xyz_ValidateFunction> functions, functions_new;
  // functions = xyz_change_statistics_generate(terms);
  functions = xyz_change_statistics_generate_new(terms);
  
  arma::vec at_zero;
  at_zero = xyz_eval_at_empty_network_new(terms, object);
  // Rcout << xyz_change_statistics_generate_new(terms) << std::endl;
  // Rcout << functions << std::endl;
  // Rcout << at_zero << std::endl;
  // Rcout <<at_zero_new << std::endl;
  
  arma::vec global_stats(xyz_count_global_statistic(object,
                                                    data_list,
                                                    type_list,
                                                    functions, 
                                                    type_x, 
                                                    type_y, 
                                                    attr_x_scale, 
                                                    attr_y_scale));
  global_stats = global_stats + at_zero;
  
  // arma::vec global_stats(xyz_count_global_statistic(object,
  //                                                   data_list,
  //                                                   type_list,
  //                                                   functions, 
  //                                                   type_x, 
  //                                                   type_y, 
  //                                                   attr_x_scale, 
  //                                                   attr_y_scale));
  // global_stats = global_stats + at_zero;
  // arma::vec global_stats;
  // Rcout << global_stats_new << std::endl;
  // Rcout <<global_stats << std::endl;
  
  return(global_stats);
}

arma::vec xyz_count_global_internal(XYZ_class object,
                                    std::vector<std::string> terms,
                                    int n_actor,
                                    std::vector<arma::mat> &data_list,
                                    std::vector<double> &type_list, 
                                    std::string type_x, 
                                    std::string type_y, 
                                    double attr_x_scale, 
                                    double attr_y_scale) {
  std::vector<xyz_ValidateFunction> functions;
  functions = xyz_change_statistics_generate_new(terms);
  arma::vec at_zero;
  // Rcout << "at_zero" << std::endl;
  at_zero = xyz_eval_at_empty_network_new(terms, object);
  // Rcout << "at_zero" << std::endl;
  arma::vec global_stats(xyz_count_global_statistic(object,
                                                    data_list,
                                                    type_list,
                                                    functions, type_x, 
                                                    type_y, 
                                                    attr_x_scale, 
                                                    attr_y_scale));
  global_stats = global_stats + at_zero;
  // Rcout << "at_zero" << std::endl;
  // arma::vec global_stats;
  return(global_stats);
}



void xyz_simulate_network_consecutive_mh( const arma::vec &coef,
                                          XYZ_class &object,
                                          const int seed,
                                          const std::vector<arma::mat> &data_list,
                                          const std::vector<double> &type_list,
                                          const bool &is_full_neighborhood,
                                          const std::vector<xyz_ValidateFunction> &functions,
                                          arma::vec &global_stats, 
                                          const double offset_nonoverlap) {
  int n_proposals = object.n_actor * (object.n_actor -1)/(object.z_network.directed ? 1 : 2);
  
  std::string z = "z";
  // double MR;
  arma::mat HR;
  // change_staCt = arma::vec(n_elem);
  arma::vec change_stat(functions.size());
  set_seed(seed);
  
  NumericVector random_accept= runif(n_proposals,0,1);
  
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int>  distr(1, object.n_actor);
  
  arma::vec tmp_vec, tmp_stat;
  
  // Go through a loop for all actor changes
  int a = 0; 
  if(object.z_network.directed){
    for(int i = 1; i <=(object.n_actor); ++i) {
      for(int j = 1; j <=(object.n_actor); ++j) {
        if(object.overlap.at(i).count(j)){
          continue;
        } 
        if(i == j){
          continue;
        } 
        // Calculate the change stat for actor_i, actor_j from 0 to 1
        xyz_calculate_change_stats(change_stat, i,
                                                 j,
                                                 object,
                                                 data_list,
                                                 type_list,
                                                 z,
                                                 is_full_neighborhood,
                                                 functions);
        // 3. Calculate the Hastings Ratios by exp(delta(tmp_entry)*coef)
        tmp_stat=change_stat;
        
        HR= exp(coef.t()*tmp_stat +offset_nonoverlap)/(exp(coef.t()*tmp_stat +offset_nonoverlap)+1);
        // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
        if(random_accept(a)<HR.at(0)){
          if(object.z_network.get_val(i,j) == 0){
            object.z_network.add_edge(i,j);
            global_stats += tmp_stat;
          }
        } else {
          if(object.z_network.get_val(i,j)){
            global_stats -= tmp_stat;  
            object.z_network.delete_edge(i,j);
          }
        }
        ++ a; 
      }
    }  
  } else {
    for(int i = 1; i <=(object.n_actor-1); ++i) {
      for(int j = i+1; j <=(object.n_actor); ++j) {
        if(object.overlap.at(i).count(j)){
          continue;
        } 
        // Calculate the change stat for actor_i, actor_j from 0 to 1
        xyz_calculate_change_stats(change_stat, i,
                                                 j,
                                                 object,
                                                 data_list,
                                                 type_list,
                                                 z,
                                                 is_full_neighborhood,
                                                 functions);
        // 3. Calculate the Hastings Ratios by exp(delta(tmp_entry)*coef)
        tmp_stat=change_stat;
        HR= exp(coef.t()*tmp_stat +offset_nonoverlap)/(exp(coef.t()*tmp_stat +offset_nonoverlap)+1);
        // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
        if(random_accept(a)<HR.at(0)){
          if(object.z_network.get_val(i,j) == 0){
            object.z_network.add_edge(i,j);
            global_stats += tmp_stat;
          } 
        } else {
          if(object.z_network.get_val(i,j)){
            global_stats -= tmp_stat;  
            object.z_network.delete_edge(i,j);
          } 
        }
        ++ a; 
      }
    }
  }
  
}


void xyz_simulate_network_consecutive_popularity_mh( const arma::vec &coef_nonpopularity,
                                                     const arma::vec &coef_popularity,
                                                     XYZ_class &object,
                                                     const int seed,
                                                     const std::vector<arma::mat> &data_list,
                                                     const std::vector<double> &type_list,
                                                     const bool &is_full_neighborhood,
                                                     const std::vector<xyz_ValidateFunction> &functions,
                                                     arma::vec &global_stats, 
                                                     const double offset_nonoverlap) {
  int n_proposals = object.n_actor * (object.n_actor -1)/(object.z_network.directed ? 1 : 2);
  std::string z = "z";
  // double MR;
  arma::mat HR;
  // change_staCt = arma::vec(n_elem);
  arma::vec change_stat(functions.size());
  set_seed(seed);
  
  NumericVector random_accept= runif(n_proposals,0,1);
  // Go through a loop for all actor changes
  int a = 0; 
  if(object.z_network.directed){
    for(int i = 1; i <=(object.n_actor); ++i) {
      double coef_popularity_i = coef_popularity.at(i-1); 
      for(int j = 1; j <=(object.n_actor); ++j) {

        if(object.overlap.at(i).count(j)){
          continue;
        }  
        if(i==j){
          continue;
        }  
        // Calculate the change stat for actor_i, actor_j from 0 to 1
        xyz_calculate_change_stats(change_stat, i,
                                                 j,
                                                 object,
                                                 data_list,
                                                 type_list,
                                                 z,
                                                 is_full_neighborhood,
                                                 functions);
        // Rcout << "Got CS" << std::endl;
        
        // 3. Calculate the Hastings Ratios by exp(delta(tmp_entry)*coef)
        HR= exp(coef_nonpopularity.t()*change_stat + offset_nonoverlap +
          (coef_popularity_i + coef_popularity.at(j-1+object.n_actor)));
        HR= HR/(HR+1);
        // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
        if(random_accept(a)<HR.at(0)){
          if(object.z_network.get_val(i,j) == 0){
            object.z_network.add_edge(i,j);
            global_stats += change_stat;
          }
        } else {
          if(object.z_network.get_val(i,j)){
            global_stats -= change_stat;  
            object.z_network.delete_edge(i,j);
          }
        }
        ++ a; 
      }
    }
  } else {
    for(int i = 1; i <=(object.n_actor-1); ++i) {
      double coef_popularity_i = coef_popularity.at(i-1); 
      for(int j = i+1; j <=(object.n_actor); ++j) {
        if(object.overlap.at(i).count(j)){
          continue;
        }  
        // Calculate the change stat for actor_i, actor_j from 0 to 1
        xyz_calculate_change_stats(change_stat, i,
                                                 j,
                                                 object,
                                                 data_list,
                                                 type_list,
                                                 z,
                                                 is_full_neighborhood,
                                                 functions);
        // 3. Calculate the Hastings Ratios by exp(delta(tmp_entry)*coef)
        HR= exp(coef_nonpopularity.t()*change_stat + offset_nonoverlap +
          (coef_popularity_i + coef_popularity.at(j-1)));
        
        HR= HR/(HR+1);
        // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
        if(random_accept(a)<HR.at(0)){
          if(object.z_network.get_val(i,j) == 0){
            object.z_network.add_edge(i,j);
            global_stats += change_stat;
          }
        } else {
          if(object.z_network.get_val(i,j)){
            global_stats -= change_stat;  
            object.z_network.delete_edge(i,j);
          }
        }
        ++ a; 
      }
    }
  }
  
}




void xyz_simulate_network_mh( const arma::vec coef,
                              XYZ_class &object,
                              const int &n_proposals,
                              const int seed,
                              const std::vector<arma::mat> &data_list,
                              const std::vector<double> &type_list,
                              const bool &is_full_neighborhood,
                              const std::vector<xyz_ValidateFunction> &functions,
                              arma::vec &global_stats, 
                              const double offset_nonoverlap) {
  // Set up objects
  // Plan is to implement TNT sampling in the future ...
  // where edges are selected with prob. 0.5
  // std::unordered_set<char> nonzero_edges;
  if(n_proposals == 0){
    return;
  }
  int proposed_change;
  std::string z = "z";
  arma::mat HR;
  arma::vec change_stat(functions.size());
  
  set_seed(seed);
  
  NumericVector random_accept= runif(n_proposals,0,1);
  // NumericVector random_propose= runif(n_proposals,0,1);
  arma::vec tmp_vec, tmp_stat;
  arma::umat comp_vec;
  int multiplier = 1;
  int tmp_i, tmp_j, proposal_idx, tmp_switch; 
  // bool is_wrong = true;
  
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int>  proposal_nnn(0, object.overlap_mat.n_rows-1);
  // std::bernoulli_distribution proposal_within_nb(prob_nb);
  
  
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals-1); a = a + 1 ) {
    proposal_idx = proposal_nnn(generator);
    tmp_i = object.overlap_mat(proposal_idx,0);
    tmp_j = object.overlap_mat(proposal_idx,1);
    if(!object.z_network.directed){
      if(tmp_i>tmp_j){
        tmp_switch = tmp_j; 
        tmp_j = tmp_i; 
        tmp_i = tmp_switch;
      }
    }
    
    // If the entries are the same (implying loops) we skip the proposal
    // (actor_i==actor_j) ? continue;
    if(object.z_network.get_val(tmp_i,tmp_j)){
      proposed_change = 0;
      multiplier = -1;
    }  else {
      proposed_change = 1;
      multiplier = 1;
    }
    // Calculate the change stat for actor_i, actor_j from 0 to 1
    xyz_calculate_change_stats(change_stat, tmp_i,
                                             tmp_j,
                                             object,
                                             data_list,
                                             type_list,
                                             z,
                                             is_full_neighborhood,
                                             functions);
    
    // 3. Calculate the Hastings Ratios by exp(delta(tmp_entry)*coef)
    tmp_stat=change_stat*multiplier;
    
    if(object.overlap.at(tmp_i).count(tmp_j)){
      HR= exp(coef.t()*tmp_stat);
    } else {
      HR= exp(coef.t()*tmp_stat +offset_nonoverlap*multiplier);
    }
    
    // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
    // clock.tick("Change val");
    if(random_accept(a)<HR.at(0)){
      global_stats += tmp_stat;
      // Here we modify the network
      if(proposed_change == 0){
        object.z_network.delete_edge(tmp_i,tmp_j);
      } 
      if(proposed_change == 1){
        object.z_network.add_edge(tmp_i,tmp_j);
      }
    }
  }
}

void xyz_simulate_network_mh_popularity( const arma::vec coef_nonpopularity,
                                         const arma::vec coef_popularity,
                                         XYZ_class &object,
                                         const int &n_proposals,
                                         const int seed,
                                         const std::vector<arma::mat> &data_list,
                                         const std::vector<double> &type_list,
                                         const bool &is_full_neighborhood,
                                         const std::vector<xyz_ValidateFunction> &functions,
                                         arma::vec &global_stats, 
                                         const double offset_nonoverlap) {
  // Set up objects
  int proposed_change;
  if(n_proposals == 0){
    return;
  }
  std::string z = "z";
  arma::mat HR;
  arma::vec change_stat(functions.size());
  set_seed(seed);
  
  NumericVector random_accept= runif(n_proposals,0,1);
  arma::vec tmp_vec, random_i(n_proposals), random_j(n_proposals), tmp_stat;
  arma::umat comp_vec;
  
  int multiplier = 1;
  int tmp_i, tmp_j, proposal_idx, tmp_switch; 
  // bool is_wrong = true;
  
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int>  distr(0, object.overlap_mat.n_rows-1);
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals-1); a = a + 1 ) {
    // Rcout << "Proposal " << a+1 << " out of " << n_proposals << "\n";
    proposal_idx = distr(generator);
    // Rcout << "Proposal index: " << proposal_idx << "\n";
    // Rcout << "Proposal index: " << object.overlap_mat.n_rows << "\n";
    
    tmp_i = object.overlap_mat(proposal_idx,0);
    tmp_j = object.overlap_mat(proposal_idx,1);
    if(!object.z_network.directed){
      if(tmp_i>tmp_j){
        tmp_switch = tmp_j; 
        tmp_j = tmp_i; 
        tmp_i = tmp_switch;
      }
    }
    // Rcout << "Proposed edge: " << tmp_i << " " << tmp_j << "\n";
    // If the entries are the same (implying loops) we skip the proposal
    // (actor_i==actor_j) ? continue;
    if(object.z_network.get_val(tmp_i,tmp_j)){
      proposed_change = 0;
      multiplier = -1;
    }  else { 
      proposed_change = 1;
      multiplier = 1;
    } 
    // Calculate the change stat for actor_i, actor_j from 0 to 1
    xyz_calculate_change_stats(change_stat, tmp_i,
                                             tmp_j,
                                             object,
                                             data_list,
                                             type_list,
                                             z,
                                             is_full_neighborhood,
                                             functions);
    // 3. Calculate the Hastings Ratios by exp(delta(tmp_entry)*coef)
    tmp_stat=change_stat*multiplier;
    
    if(object.overlap.at(tmp_i).count(tmp_j)){
      if(object.z_network.directed){
        HR= exp(coef_nonpopularity.t()*tmp_stat + 
          multiplier*(coef_popularity.at(tmp_i-1) + coef_popularity.at(tmp_j-1 + object.n_actor)));
      } else {
        HR= exp(coef_nonpopularity.t()*tmp_stat + 
          multiplier*(coef_popularity.at(tmp_i-1) + coef_popularity.at(tmp_j-1)));  
      }
    } else {
      if(object.z_network.directed){
        HR= exp(coef_nonpopularity.t()*tmp_stat + offset_nonoverlap*multiplier +
          multiplier*(coef_popularity.at(tmp_i-1) + coef_popularity.at(tmp_j-1+ object.n_actor)));
      } else {
        HR= exp(coef_nonpopularity.t()*tmp_stat + offset_nonoverlap*multiplier +
          multiplier*(coef_popularity.at(tmp_i-1) + coef_popularity.at(tmp_j-1)));
      }
      
    }
    // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
    if(random_accept(a)<HR.at(0)){
      global_stats += tmp_stat;
      // Here we modify the network
      if(proposed_change == 0){
        object.z_network.delete_edge(tmp_i,tmp_j);
      } 
      if(proposed_change == 1){
        object.z_network.add_edge(tmp_i,tmp_j);
      }
    }
  }
}

// void xyz_simulate_network_mh_popularity(
//     const arma::vec& coef_nonpopularity, // Use const&
//     const arma::vec& coef_popularity,   // Use const&
//     XYZ_class &object,
//     int n_proposals,                   // Pass by value
//     int seed,
//     std::vector<arma::mat>& data_list, // Use const&
//     std::vector<double>& type_list,   // Use const&
//     bool is_full_neighborhood,             // Pass by value
//     const std::vector<xyz_ValidateFunction>& functions, // Use const&
//     arma::vec &global_stats,
//     double offset_nonoverlap)
// {
//   if (n_proposals <= 0) {
//     return;
//   }
//   
//   std::string z = "z";
//   arma::vec change_stat;
//   const int n_actor = object.n_actor; // Cache n_actor
//   
//   // --- C++ Random Number Generation ---
//   std::mt19937 generator(seed);
//   std::uniform_int_distribution<int> node_distr(1, );
//   std::uniform_real_distribution<double> accept_distr(0.0, 1.0);
//   
//   // --- Main MH Loop ---
//   for (int a = 0; a < n_proposals; ++a) {
//     int tmp_i, tmp_j;
//     bool is_within_proposal = false; // For NNN sampling
//     
//     // --- Dyad Proposal Loop ---
//     while (true) { // Loop until a valid dyad is found
//       tmp_i = node_distr(generator);
//       tmp_j = node_distr(generator);
//       
//       if (tmp_i == tmp_j) continue; // Skip self-loops
//       
//       // Ensure i < j for undirected graphs to avoid proposing same dyad twice
//       if (!object.z_network.directed && tmp_i > tmp_j) {
//         std::swap(tmp_i, tmp_j);
//       }
//       
//       // Check NNN sampling constraints if applicable
//       if (NNN_sampling) {
//         // Check if actors exist in overlap map before accessing .at()
//         if (object.overlap.find(tmp_i) == object.overlap.end()){
//           // Rcpp::Rcout << "Warning: Actor " << tmp_i << " not found in overlap map. Skipping proposal." << std::endl;
//           continue; // Skip if actor not in map
//         }
//         
//         is_within_proposal = proposal_within_nb_distr(generator);
//         bool is_in_overlap = object.overlap.at(tmp_i).count(tmp_j); // Assuming overlap check is efficient
//         
//         if (is_within_proposal && is_in_overlap) break; // Valid within-NB proposal
//         if (!is_within_proposal && !is_in_overlap) break; // Valid outside-NB proposal
//         // Otherwise, continue loop to find a matching proposal type
//       } else {
//         break; // Valid proposal if not NNN sampling
//       }
//     } // End dyad proposal loop
//     
//     // --- Calculate Change Statistics ---
//     bool current_edge_exists = object.z_network.get_val(tmp_i, tmp_j);
//     int multiplier = current_edge_exists ? -1 : 1; // -1 for removal, +1 for addition proposal
//     
//     // Calculate change_stat assuming the toggle happens (stat for ADDITION)
//     change_stat = xyz_calculate_change_stats(tmp_i, tmp_j, object, data_list, type_list, z, is_full_neighborhood, functions);
//     
//     // This is the change stat vector for the *proposed* move (add or delete)
//     arma::vec delta_stats = change_stat * multiplier;
//     
//     // --- Calculate Log Hastings Ratio ---
//     double log_HR = arma::dot(coef_nonpopularity, delta_stats);
//     
//     // Add popularity terms (multiplied by multiplier)
//     if (object.z_network.directed) {
//       log_HR += multiplier * (coef_popularity.at(tmp_i - 1) + coef_popularity.at(tmp_j - 1 + n_actor));
//     } else {
//       log_HR += multiplier * (coef_popularity.at(tmp_i - 1) + coef_popularity.at(tmp_j - 1));
//     }
//     
//     // Add offset term if the dyad is NOT in the overlap set
//     // Check if actors exist in overlap map before accessing .at()
//     bool dyad_in_overlap = false;
//     if (object.overlap.find(tmp_i) != object.overlap.end()){
//       dyad_in_overlap = object.overlap.at(tmp_i).count(tmp_j);
//     }
//     
//     if (!dyad_in_overlap) {
//       log_HR += offset_nonoverlap * multiplier;
//     }
//     
//     // --- Acceptance Check ---
//     // Accept if log_HR >= 0 (HR >= 1) or exp(log_HR) > U(0,1)
//     // Using log comparison is often more stable and avoids exp when HR >= 1
//     double log_u = std::log(accept_distr(generator));
//     if (log_HR >= log_u) { // Accept the move
//       global_stats += delta_stats; // Update global stats
//       
//       // Apply the change to the network object
//       if (multiplier == 1) { // We proposed adding
//         object.z_network.add_edge(tmp_i, tmp_j);
//       } else { // We proposed deleting
//         object.z_network.delete_edge(tmp_i, tmp_j);
//       }
//     }
//     // If rejected, do nothing to global_stats or object.z_network
//     
//   } // End main MH loop
// }



void xyz_simulate_attribute_mh( const arma::vec coef,
                                XYZ_class &object,
                                const int &n_proposals,
                                const int seed,
                                const  std::vector<arma::mat> &data_list,
                                const std::vector<double> &type_list,
                                const bool &is_full_neighborhood,
                                const std::vector<xyz_ValidateFunction> &functions,
                                arma::vec &global_stats, 
                                const std::string type) {
  if(n_proposals == 0){
    return;
  }
  int proposed_change;
  arma::vec HR;
  set_seed(seed);
  std::mt19937 gen(seed); // Mersenne Twister random generator
  // std::uniform_real_distribution<> dist(0.0, 1.0);  // Uniform distribution [0, 1)
  arma::vec change_stat(functions.size());
  
  NumericVector random_accept= runif(n_proposals,0,1);
  arma::vec tmp_stat(functions.size());
  int multiplier;
  int tmp_i; 
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int>  distr(1, object.n_actor);
  const double MAX_LOG_RATE = 100.0;
  arma::vec tmp;
  arma::vec tmp_row;
  // Go through a loop for each proposed change
  for(int a = 0; a <=(n_proposals-1); a ++ ) {
    // Here we pick the random entry
    tmp_i = distr(generator);
    // Here we calculate the change stat from turning y_i from 0 to 1
    xyz_calculate_change_stats(change_stat, tmp_i,
                                             tmp_i,
                                             object,
                                             data_list,
                                             type_list,
                                             type,
                                             is_full_neighborhood,
                                             functions);
    // Rcpp::Rcout << tmp_stat << std::endl;
    // Rcpp::Rcout << global_stats <<std::endl;
    if(type == "x"){
      if(object.x_attribute.type == "binomial"){
        if(object.x_attribute.get_val(tmp_i)){
          proposed_change = 0;
          multiplier = -1;
        } else {
          proposed_change = 1;
          multiplier = 1;
        }  
        // 3. Step: Calculate the Hastings Ratios
        tmp_stat=change_stat*multiplier;
        HR= exp(coef.t()*tmp_stat);
        
        // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
        if(random_accept(a)<HR.at(0)){
          global_stats += tmp_stat;
          // Here we modify the network
          if(proposed_change == 0){
            object.x_attribute.set_attr_0(tmp_i);  
          } else {
            object.x_attribute.set_attr_1(tmp_i);
          }
        }
      }
      if(object.x_attribute.type == "poisson"){
        tmp = coef.t()*change_stat;
        // Rcout << tmp << std::endl;
        double safe_eta = std::min(tmp.at(0), MAX_LOG_RATE);
        double tmp_val = R::rpois(exp(safe_eta)); 
        global_stats += (tmp_val- object.x_attribute.get_val(tmp_i))*change_stat;
        object.x_attribute.set_attr_value(tmp_i, tmp_val);  
        
        
      }
      if(object.x_attribute.type == "normal"){
        HR= coef.t()*change_stat;
        double tmp_val = R::rnorm(HR.at(0), 1); 
        global_stats += (tmp_val- object.x_attribute.get_val(tmp_i))*change_stat;
        object.x_attribute.set_attr_value(tmp_i, tmp_val);  
      }
    }
    if(type == "y"){
      // Rcpp::Rcout << object.y_attribute.get_val(tmp_i) <<std::endl;
      
      if(object.y_attribute.type == "binomial"){
        if(object.y_attribute.get_val(tmp_i)){
          proposed_change = 0;
          multiplier = -1;
        } else {
          proposed_change = 1;
          multiplier = 1;
        }
        
        // 3. Step: Calculate the Hastings Ratios
        tmp_stat=change_stat*multiplier;
        HR= exp(coef.t()*tmp_stat);
        
        // 4. Step: Sample a random number between 0 and 1, accept if it is > HR
        if(random_accept(a)<HR.at(0)){
          global_stats += tmp_stat;
          // Here we modify the network
          if(proposed_change == 0){
            object.y_attribute.set_attr_0(tmp_i);  
          } else {
            object.y_attribute.set_attr_1(tmp_i);
          }
        }
      }
      if(object.y_attribute.type == "poisson"){
        tmp = coef.t()*change_stat;
        // Rcout << tmp << std::endl;
        double safe_eta = std::min(tmp.at(0), MAX_LOG_RATE);
        double tmp_val = R::rpois(exp(safe_eta)); 
        global_stats +=  (tmp_val- object.y_attribute.get_val(tmp_i))*change_stat;
        object.y_attribute.set_attr_value(tmp_i, tmp_val);  
      }
      if(object.y_attribute.type == "normal"){
        HR= coef.t()*change_stat;
        double tmp_val = R::rnorm(HR.at(0), 1); 
        global_stats += (tmp_val- object.y_attribute.get_val(tmp_i))*change_stat;
        object.y_attribute.set_attr_value(tmp_i, tmp_val);  
      }
      
    }
    
    
    
    
  }
}

arma::mat xyz_simulate_internal(XYZ_class & object,
                                const arma::vec& coef,
                                const  arma::vec& coef_popularity,
                                const std::vector<arma::mat>& data_list,
                                const std::vector<double>& type_list,
                                arma::vec & global_stats,
                                const int n_proposals_x,
                                const int seed_x,
                                const int n_proposals_y,
                                const int seed_y,
                                const  int n_proposals_z,
                                const int seed_z,
                                const int n_burn_in,
                                const  int n_simulation,
                                std::vector<arma::vec>& res_x,
                                std::vector<arma::vec>& res_y,
                                std::vector<std::unordered_map< int, std::unordered_set<int>>>& res_z,
                                const bool only_stats,
                                const bool is_full_neighborhood,
                                const std::vector<xyz_ValidateFunction> functions,
                                const bool display_progress, 
                                const bool popularity, 
                                const double offset_nonoverlap, 
                                const bool fix_x = false){
  arma::mat stats(n_simulation,functions.size());
  stats.fill(0);
  std::string x, y; 
  x = "x";
  y = "y";
  Progress p(n_simulation + n_burn_in, display_progress);
  // Start for a burn in period with the normal number of proposals
  // Intialize global statistics and then adapt them peu a peu
  for(int i = 1; i <=(n_simulation + n_burn_in);i ++) {
    // Simulate the network
    // Rcout << global_stats << std::endl;
    
    p.increment(); // update progress
    if(!fix_x){
      // Rcout << "Sampling X| Y,Z" << std::endl;
      // Sample X| Y,Z
      xyz_simulate_attribute_mh(coef,object,
                                n_proposals_x,seed_x +i,
                                data_list, 
                                type_list,
                                is_full_neighborhood, 
                                functions,
                                global_stats, x);  
    }
    // Rcout << "Sampling Y| X,Z" << std::endl;
    // Sample Y| X,Z
    xyz_simulate_attribute_mh(coef,object,
                              n_proposals_y,seed_y +i,
                              data_list, type_list,
                              is_full_neighborhood, functions,
                              global_stats, y);
    
    // Rcout << global_stats << std::endl;
    
    
    // Sample Z|X,Y
    // Rcout << "Sampling Z| X,Y" << std::endl;
    if(popularity){
      xyz_simulate_network_mh_popularity(coef,
                                         coef_popularity,
                                         object,
                                         n_proposals_z, seed_z +i,
                                         data_list, type_list,
                                         is_full_neighborhood, functions,
                                         global_stats, offset_nonoverlap); 
    } else {
      xyz_simulate_network_mh(coef,object,
                              n_proposals_z, seed_z +i,
                              data_list, type_list,
                              is_full_neighborhood, functions,
                              global_stats,  offset_nonoverlap);  
    }
    // Rcout << "Sampling Z non-overlap | X,Y" << std::endl;
    if(popularity){
      xyz_simulate_network_consecutive_popularity_mh(coef,
                                                     coef_popularity,object,
                                                     seed_z*2 +i,
                                                     data_list, type_list,
                                                     is_full_neighborhood, functions,
                                                     global_stats, offset_nonoverlap);
    } else {
      xyz_simulate_network_consecutive_mh(coef,object,
                                          seed_z*2 +i,
                                          data_list, type_list,
                                          is_full_neighborhood, functions,
                                          global_stats, offset_nonoverlap);
    }
    // We throw the first n_burn_in samples away
    if(i>n_burn_in){
      if(only_stats){
        // Count global statistics
        stats.row(i - n_burn_in-1) = global_stats.as_row();
      } else{
        // Save the object
        // if(!fix_x){
        //   res_x.at(i - n_burn_in-1) =object.x_attribute.attribute; 
        // }
        res_x.at(i - n_burn_in-1) =object.x_attribute.attribute; 
        res_y.at(i - n_burn_in-1) =object.y_attribute.attribute;
        res_z.at(i - n_burn_in-1) =object.z_network.adj_list;
        // Count global statistics
        stats.row(i - n_burn_in-1) = global_stats.as_row();
      }
    }
    
  }
  return(stats);
}


// [[Rcpp::export]]
List xyz_simulate_cpp(arma::vec& coef,
                      arma::vec& coef_popularity,
                      std::vector<std::string>& terms,
                      int& n_actor,
                      arma::mat z_network,
                      arma::mat neighborhood,
                      arma::mat overlap,
                      arma::vec x_attribute,
                      arma::vec y_attribute,
                      bool init_empty,
                      bool directed,
                      bool popularity,
                      std::vector<arma::mat>& data_list,
                      std::vector<double>& type_list,
                      double offset_nonoverlap,
                      std::string type_x, 
                      std::string type_y, 
                      double attr_x_scale, 
                      double attr_y_scale,
                      int n_proposals_x = 100,
                      int seed_x = 123,
                      int n_proposals_y = 100,
                      int seed_y = 123,
                      int n_proposals_z = 100,
                      int seed_z = 123,
                      int n_proposals_z_nonoverlap = 100,
                      int seed_z_nonoverlap = 123,
                      int n_burn_in = 100,
                      int n_simulation = 1,
                      bool only_stats = false,
                      bool display_progress = false, 
                      bool fix_x = false){
  // res(n_simulation);
  // stats2.fill(0);
  XYZ_class object(n_actor,directed, neighborhood, overlap, type_x, type_y,attr_x_scale, attr_y_scale);
  if(!init_empty){
    object.set_info_arma(x_attribute,y_attribute, z_network);
  }
  bool is_full_neighborhood = object.check_if_full_neighborhood();
  // Rcout <<object.overlap.at(1)<< std::endl;
  
  // object.set_neighborhood_from_mat(neighborhood);
  // Generate change statistic function from the terms
  std::vector<xyz_ValidateFunction> functions;
  
  functions = xyz_change_statistics_generate_new(terms);
  // arma::vec at_zero;
  // at_zero = xyz_eval_at_empty_network_new(terms);
  // // Start for a burn in period with the normal number of proposals
  // // Intialize global statistics and then adapt them peu a peu
  // arma::vec global_stats(functions.size());
  // global_stats.fill(0);
  // global_stats = global_stats + at_zero;
  // Rcout << "Here"<< std::endl;
  arma::vec global_stats = xyz_count_global_internal( object,
                                                      terms,
                                                      n_actor,
                                                      data_list,
                                                      type_list,
                                                      type_x, type_y, 
                                                      attr_x_scale, 
                                                      attr_y_scale);
  // Rcout << global_stats << std::endl;
  std::vector<arma::vec> res_x(n_simulation);
  std::vector<arma::vec> res_y(n_simulation);
  std::vector<std::unordered_map< int, std::unordered_set<int>>> res_z(n_simulation);
  // Rcout << "B"<< std::endl;
  arma::mat stats = xyz_simulate_internal(object, coef,coef_popularity, data_list, type_list, global_stats,
                                          n_proposals_x, seed_x,
                                          n_proposals_y, seed_y,
                                          n_proposals_z, seed_z,
                                          n_burn_in, n_simulation,
                                          res_x,res_y,res_z,
                                          only_stats, 
                                          is_full_neighborhood, 
                                          functions, 
                                          display_progress, 
                                          popularity,
                                          offset_nonoverlap, 
                                          fix_x);
  if(only_stats){
    return(List::create(_["stats"] = stats));
  } else {
    return(List::create(_["simulation_attributes_x"] =res_x,_["simulation_attributes_y"] =res_y,
                        _["simulation_networks_z"] =res_z, _["stats"] = stats));  
  }
}
std::tuple<arma::mat, arma::vec> xyz_get_info_pl(XYZ_class object,
                                                 std::vector<std::string> terms,
                                                 std::vector<arma::mat> &data_list,
                                                 std::vector<double> &type_list, 
                                                 bool display_progress, 
                                                 arma::uvec &i_vec,
                                                 arma::uvec &j_vec, 
                                                 arma::uvec &overlap_vec, 
                                                 int n_actor, 
                                                 bool fix_x) {
  bool is_full_neighborhood = object.check_if_full_neighborhood();
  // Generate vector of functions that calculate the sufficient statistics 
  std::vector<xyz_ValidateFunction> functions;
  functions = xyz_change_statistics_generate_new(terms);
  // strings z, x, and y later needed to tell the sufficient statistics 
  // what type of change statistic is wanted
  std::string z = "z", x = "x", y = "y";
  // Just a temporary vector of the change statistics for one dyad of the network
  arma::vec change_stat_network(functions.size()),
  // The same thing but for the attributes i and j
  change_stat_attribute_i(functions.size()), change_stat_attribute_j(functions.size());
  int x_i, y_i, z_ij;
  // int ncores = 5;
  arma::mat res_covs;
  arma::vec res_target;
  Progress p(n_actor*(n_actor-1) + 2*n_actor, display_progress);
  
  // Rcpp::Rcout << "Rows: " << n_actor*(n_actor-1) + n_actor << std::endl;
  // Rcpp::Rcout << "Colum: " <<terms.size() << std::endl;
  if(object.z_network.directed){
    // Differentiate between the case that x and y or just y is assumed to be random
    if(fix_x){
      res_covs.reshape(n_actor*(n_actor-1) + n_actor,terms.size());
      res_target.reshape(n_actor*(n_actor-1) +n_actor,1);
    } else {
      res_covs.reshape(n_actor*(n_actor-1) + n_actor*2,terms.size());
      res_target.reshape(n_actor*(n_actor-1) +n_actor*2,1);
    }
  } else {
    if(fix_x){
      res_covs.reshape(n_actor*(n_actor-1)/2 + n_actor,terms.size());
      res_target.reshape(n_actor*(n_actor-1)/2 +n_actor,1);  
    } else {
      res_covs.reshape(n_actor*(n_actor-1)/2 + n_actor*2,terms.size());
      res_target.reshape(n_actor*(n_actor-1)/2 +n_actor*2,1);    
    }
    
  }
  int now = 0;
  if(object.z_network.directed){
    for(int i: seq(1,n_actor)){
      Rcpp::checkUserInterrupt();
      
      for(int j: seq(1,n_actor)){
        // Get present values of x_ij, y_i, y_j
        p.increment(); // update progress
        if(i == j){
          continue;
        }
        z_ij = object.z_network.get_val(i,j);
        xyz_calculate_change_stats(change_stat_network, i,
                                                         j,
                                                         object,
                                                         data_list,
                                                         type_list,
                                                         z,
                                                         is_full_neighborhood,
                                                         functions);
        // change_stat_network(change_stat_network.size()+1) =x_ij;
        res_covs.row(now)= (change_stat_network).as_row();
        res_target.at(now) = z_ij;
        i_vec.at(now) = i;
        j_vec.at(now) = j;
        overlap_vec.at(now) = object.overlap.at(i).count(j);
        now += 1;
      } 
    }
  } else {
    for(int i: seq(1,n_actor-1)){
      for(int j: seq(i+1,n_actor)){
        // Get present values of x_ij, y_i, y_j
        p.increment(); // update progress
        
        z_ij = object.z_network.get_val(i,j);
        xyz_calculate_change_stats(change_stat_network, i,
                                                         j,
                                                         object,
                                                         data_list,
                                                         type_list,
                                                         z,
                                                         is_full_neighborhood,
                                                         functions);
        // change_stat_network(change_stat_network.size()+1) =x_ij;
        res_covs.row(now)= (change_stat_network).as_row();
        res_target.at(now) = z_ij;
        i_vec.at(now) = i;
        j_vec.at(now) = j;
        overlap_vec.at(now) = object.overlap.at(i).count(j);
        now += 1;
      } 
    }
  }
  // End with 2 n_actor entries of the attributes
  for(int i: seq(1,n_actor)){
    if(!fix_x){
      p.increment(); 
      x_i = object.x_attribute.get_val(i);
      xyz_calculate_change_stats(change_stat_attribute_i, i,
                                                           i,
                                                           object,
                                                           data_list,
                                                           type_list,
                                                           x,
                                                           is_full_neighborhood,
                                                           functions);
      res_covs.row(now)= (change_stat_attribute_i).as_row();
      res_target.at(now) = x_i;
      now += 1;  
    }
    p.increment(); 
    y_i = object.y_attribute.get_val(i);
    xyz_calculate_change_stats(change_stat_attribute_i, i,
                                                         i,
                                                         object,
                                                         data_list,
                                                         type_list,
                                                         y,
                                                         is_full_neighborhood,
                                                         functions);
    res_covs.row(now)= (change_stat_attribute_i).as_row();
    res_target.at(now) = y_i;
    now += 1;
  } 
  // Rcpp::Rcout << res_target << std::endl;
  return(std::tuple<arma::mat, arma::vec> {res_covs, res_target});
} 

// Version where the argument of a XYZ_class object 
// and not the separate attributes and network as in the function xyz_prepare_composite_estimation_internal
// std::vector<arma::mat> xyz_get_info(XYZ_class object,
//                                     std::vector<std::string> terms,
//                                     std::vector<arma::mat> &data_list,
//                                     std::vector<double> &type_list, 
//                                     bool add_info, 
//                                     bool display_progress, 
//                                     arma::vec overlap_vec) {
//   // Set up objects
//   int n_actor = object.n_actor;
//   std::vector<xyz_ValidateFunction> functions;
//   functions = xyz_change_statistics_generate(terms);
//   std::string z = "z", x = "x", y = "y";
//   // arma::vec change_stat_x_i(functions.size()),
//   // change_stat_x_j(functions.size()), change_stat_y_i(functions.size()),
//   // change_stat_y_j(functions.size()), change_stat_z_ij(functions.size());
//   bool is_full_neighborhood = object.check_if_full_neighborhood();
//   
//   bool x_i, x_j, z_ij, y_i,y_j;
//   int now = -1;
//   
//   if(object.z_network.directed){
//     std::vector<arma::mat> res(n_actor*(n_actor-1));
//     Progress p_3(n_actor*(n_actor-1), display_progress);
//     
//     
//     for(int i: seq(1,n_actor)){
//       for(int j: seq(1,n_actor)){
//         if(i == j){
//           continue;
//         } else {
//           now += 1;
//           p_3.increment();
//         } 
//         // Get present values of x_ij, y_i, y_j
//         x_i = object.x_attribute.get_val(i);
//         x_j = object.x_attribute.get_val(j);
//         y_i = object.y_attribute.get_val(i);
//         y_j = object.y_attribute.get_val(j);
//         z_ij = object.z_network.get_val(i,j);
//         overlap_vec.at(now) = object.overlap.at(i).count(j);
//         if(add_info) {
//           res.at(now) = get_all_combinations_xyz(i,j,x_i, x_j,y_i,y_j,z_ij,now, 
//                  object,functions, data_list, type_list, is_full_neighborhood);
//         } else { 
//           res.at(now) = get_all_combinations_xyz(i,j,x_i, x_j,y_i,y_j,z_ij,now, 
//                  object,functions, data_list, type_list, is_full_neighborhood).submat(0,8,31,7+functions.size());
//         } 
//       } 
//     }
//     
//     return(res);
//     
//   } else {
//     std::vector<arma::mat> res(n_actor*(n_actor-1)/2);
//     Progress p_3(n_actor*(n_actor-1)/2, display_progress);
//     
//     for(int i: seq(1,n_actor-1)){
//       for(int j: seq(i+1,n_actor)){
//         
//         // 
//         // Get present values of x_ij, y_i, y_j
//         x_i = object.x_attribute.get_val(i);
//         x_j = object.x_attribute.get_val(j);
//         y_i = object.y_attribute.get_val(i);
//         y_j = object.y_attribute.get_val(j);
//         z_ij = object.z_network.get_val(i,j);
//         overlap_vec.at(now) = object.overlap.at(i).count(j);
//         
//         now += 1;
//         p_3.increment();
//         if(add_info) {
//           res.at(now) = get_all_combinations_xyz(i,j,x_i, x_j,y_i,y_j,z_ij,now, 
//                  object,functions, data_list, type_list, is_full_neighborhood);
//           
//         } else { 
//           res.at(now) = get_all_combinations_xyz(i,j,x_i, x_j,y_i,y_j,z_ij,now, 
//                  object,functions, data_list, type_list, is_full_neighborhood).submat(0,8,31,7+functions.size());
//         } 
//         
//       } 
//     }
//     p_3.cleanup();
//     return(res);
//   }
//   
// }

arma::mat get_A_exact(arma::uvec i_vec, 
                      arma::uvec  j_vec,
                      arma::uvec  overlap_vec,
                      std::tuple<arma::mat,arma::vec> &pseudo_lh, 
                      arma::vec &coef_popularity, 
                      arma::vec &coef_nonpopularity, 
                      double &offset_nonoverlap, 
                      bool directed, 
                      int n_actor){
  arma::vec res;
  arma::vec exp_tmp;
  arma::mat A(coef_popularity.size(),coef_popularity.size());
  unsigned int number_elements_network = i_vec.size();
  // For the calculation of B we only have to regard network information (relating to the first entries)
  
  for(unsigned int i = 0; i < number_elements_network; i++){
    if(overlap_vec.at(i)){
      if(directed) {
        exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1+ n_actor));  
      } else {
        exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1));  
      }
    } else {
      if(directed) {
        exp_tmp = arma::exp(offset_nonoverlap + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1 + n_actor));
      } else {
        exp_tmp = arma::exp(offset_nonoverlap + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1));
      }
    }
    res = exp_tmp.at(0)/(1+exp_tmp.at(0)) - pow(exp_tmp.at(0)/(1+exp_tmp.at(0)),2);
    
    A.at(i_vec.at(i)-1, i_vec.at(i)-1)+= res.at(0);
    if(directed) {
      A.at(j_vec.at(i)-1 + n_actor, j_vec.at(i)-1 + n_actor)+= res.at(0);
      A.at(i_vec.at(i)-1, j_vec.at(i)-1 + n_actor)+= res.at(0);
      A.at(j_vec.at(i)-1 + n_actor, i_vec.at(i)-1)+= res.at(0);
    } else {
      A.at(j_vec.at(i)-1, j_vec.at(i)-1)+= res.at(0);
      A.at(i_vec.at(i)-1, j_vec.at(i)-1)+= res.at(0);
      A.at(j_vec.at(i)-1, i_vec.at(i)-1)+= res.at(0);
    }
  }
  // Rcpp::Rcout << "Old" << std::endl;
  // Rcpp::Rcout << arma::accu(A) << std::endl;
  
  
  return(A);
}

std::tuple< arma::vec, arma::mat> get_B_pl(arma::uvec i_vec, 
                                           arma::uvec  j_vec,
                                           arma::uvec  overlap_vec,
                                           std::tuple<arma::mat,arma::vec> &pseudo_lh, 
                                           arma::vec &coef_popularity, 
                                           arma::vec &coef_nonpopularity, 
                                           double &offset_nonoverlap, 
                                           bool directed, 
                                           int n_actor) {
  arma::mat B_mat(coef_nonpopularity.size(),coef_popularity.size());
  arma::vec A_diag(coef_popularity.size()), res;
  B_mat.fill(0);
  arma::vec exp_tmp, cov_popularity(coef_popularity.size());
  cov_popularity.fill(0);
  unsigned int number_elements_network = i_vec.size();
  // For the calculation of B we only have to regard network information (relating to the first entries)
  
  for(unsigned int i = 0; i < number_elements_network; i++){
    if(overlap_vec.at(i)){
      if(directed) {
        exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1+ n_actor));  
        cov_popularity.at(j_vec.at(i)-1+ n_actor) = 1;
      } else {
        exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1));  
        cov_popularity.at(j_vec.at(i)-1) = 1;
      }
    } else {
      if(directed) {
        exp_tmp = arma::exp(offset_nonoverlap + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1 + n_actor));
        cov_popularity.at(j_vec.at(i)-1+ n_actor) = 1;
      } else {
        exp_tmp = arma::exp(offset_nonoverlap + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity +
          coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1));
        cov_popularity.at(j_vec.at(i)-1) = 1;
      }
    }
    cov_popularity.at(i_vec.at(i)-1) = 1;
    B_mat += exp_tmp.at(0)/(1+exp_tmp.at(0))*1/(1+exp_tmp.at(0))*
      std::get<0>(pseudo_lh).row(i).t()*cov_popularity.t();
    
    res = exp_tmp.at(0)/(1+exp_tmp.at(0)) - pow(exp_tmp.at(0)/(1+exp_tmp.at(0)),2);
    A_diag.at(i_vec.at(i)-1) += res.at(0);
    cov_popularity.at(i_vec.at(i)-1) = 0;
    if(directed) { 
      cov_popularity.at(j_vec.at(i)-1 + n_actor) = 0;
      A_diag.at(j_vec.at(i)-1 + n_actor) += res.at(0);
      
    } else {
      cov_popularity.at(j_vec.at(i)-1) = 0;
      A_diag.at(j_vec.at(i)-1) += res.at(0);
      
    }
    
  }
  return(std::tuple<arma::vec, arma::mat> {A_diag,B_mat});
}


arma::mat get_B(arma::uvec i_vec, arma::uvec  j_vec,
                arma::uvec  overlap_vec,
                std::tuple<arma::mat,arma::vec> &pseudo_lh, 
                arma::vec &coef_popularity,
                arma::vec &coef_nonpopularity, 
                double &offset_nonoverlap, 
                bool directed, 
                int n_actor) {
  arma::mat B_mat(coef_nonpopularity.size(),coef_popularity.size());
  B_mat.fill(0);
  arma::vec exp_tmp, cov_popularity(coef_popularity.size());
  cov_popularity.fill(0);
  unsigned int number_elements_network = i_vec.size();
  // For the calculation of B we only have to regard network information (relating to the first entries)
  
  arma::vec vec_network = arma::exp(offset_nonoverlap*(1-overlap_vec)+std::get<0>(pseudo_lh).rows(0,number_elements_network-1)*coef_nonpopularity +
    coef_popularity.elem(i_vec-1) + coef_popularity.elem(j_vec-1+ n_actor*directed));
  arma::vec vec_attribute = arma::exp(std::get<0>(pseudo_lh).rows(number_elements_network,std::get<0>(pseudo_lh).n_rows-1)*coef_nonpopularity); 
  arma::vec result = join_cols(vec_network, vec_attribute); 
  
  for(unsigned int i = 0; i < number_elements_network; i++){
  
    arma::uword i_idx = i_vec[i] - 1;
    arma::uword j_idx = j_vec[i] - 1 + (directed ? n_actor : 0);
    double weight = result.at(i) / std::pow(1.0 + result.at(i), 2);
    const arma::rowvec& x_row = std::get<0>(pseudo_lh).row(i);
    
    B_mat.col(i_idx) += weight * x_row.t();
    B_mat.col(j_idx) += weight * x_row.t();
    
  } 
  
  return(B_mat);
}



std::tuple<arma::vec, arma::vec, arma::mat, arma::mat> cond_estimation_nonpopularity_pl_old(arma::vec coef,
                                                                                            arma::uvec &i_vec,
                                                                                            arma::uvec  &j_vec,
                                                                                            arma::uvec  &overlap_vec,
                                                                                            bool directed,
                                                                                            std::tuple<arma::mat,arma::vec> &pseudo_lh,
                                                                                            int max_iteration,
                                                                                            double tol,
                                                                                            arma::vec &coef_popularity,
                                                                                            double offset_nonoverlap,
                                                                                            bool &non_stop, 
                                                                                            std::string attr_x_type, 
                                                                                            std::string attr_y_type, 
                                                                                            double x_scale, 
                                                                                            double y_scale, 
                                                                                            bool fix_x) {
  int n_coef = coef.size();
  List score_ind(std::get<0>(pseudo_lh).n_rows);
  List fisher_ind(std::get<0>(pseudo_lh).n_rows);
  
  arma::vec score(n_coef), score_tmp(n_coef), m(n_coef), exp_tmp(n_coef);
  arma::mat fisher(n_coef, n_coef),x_trans(n_coef, n_coef), coefs(max_iteration+1, n_coef);
  coefs.row(0)= coef.t();
  score.fill(0);
  fisher.fill(0);
  bool non_converged = true;
  int k = 1;
  unsigned int n_actor = coef_popularity.n_elem/2;
  unsigned int number_elements_network = i_vec.size();
  while(non_converged) {
    // Update the score and fisher info
    for(unsigned int i = 0; i < std::get<0>(pseudo_lh).n_rows; i++){
      // What to do if we are regarding network information (relating to the first entries)
      if(i < number_elements_network){
        if(directed){
          if(overlap_vec.at(i)){
            exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef +
              coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1 + n_actor));
          } else {
            exp_tmp = arma::exp(offset_nonoverlap+std::get<0>(pseudo_lh).row(i)*coef +
              coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1+ n_actor));
          }
        } else {
          if(overlap_vec.at(i)){
            exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef +
              coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1));
          } else {
            exp_tmp = arma::exp(offset_nonoverlap+std::get<0>(pseudo_lh).row(i)*coef +
              coef_popularity.at(i_vec.at(i)-1) + coef_popularity.at(j_vec.at(i)-1));
          }  
        }
        score_tmp = std::get<0>(pseudo_lh).row(i).t()*
          std::get<1>(pseudo_lh).at(i) -
          exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
        score += score_tmp;
        fisher += exp_tmp.at(0)/(pow(1+exp_tmp.at(0), 2))*
          std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
        
      } else if((i >= number_elements_network) && (i < number_elements_network +n_actor) && (fix_x == false)){
        // Rcpp::Rcout << "X" << std::endl;
        // Rcpp::Rcout << i << std::endl;
        // Rcpp::Rcout << (i >= number_elements_network) << std::endl;
        // Rcpp::Rcout <<  (i < number_elements_network +n_actor)  << std::endl;
        // Rcpp::Rcout << (fix_x == false) << std::endl;
        // 
        // What to do if we are regarding attribute information (relating to the last entries)
        if(attr_x_type == "binomial"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score_tmp = std::get<0>(pseudo_lh).row(i).t()*
            std::get<1>(pseudo_lh).at(i) -
            exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          score += score_tmp;
          fisher += exp_tmp.at(0)/(pow(1+exp_tmp.at(0), 2))*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
          
        }
        if(attr_x_type == "poisson"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          fisher += exp_tmp.at(0)*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
        } else if(attr_x_type == "normal"){
          exp_tmp = std::get<0>(pseudo_lh).row(i)*coef;
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t()/x_scale;
          fisher += std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i)/x_scale;
        }
        
        
      } else {
        // Rcpp::Rcout << "Y" << std::endl;
        // // Rcpp::Rcout << i << std::endl;
        // Rcpp::Rcout << i << std::endl;
        // Rcpp::Rcout << (i >= number_elements_network) << std::endl;
        // Rcpp::Rcout <<  (i < number_elements_network +n_actor)  << std::endl;
        // Rcpp::Rcout << (fix_x == false) << std::endl;
        
        
        if(attr_y_type == "binomial"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score_tmp = std::get<0>(pseudo_lh).row(i).t()*
            std::get<1>(pseudo_lh).at(i) -
            exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          score += score_tmp;
          fisher += exp_tmp.at(0)/(pow(1+exp_tmp.at(0), 2))*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
          
        }
        if(attr_y_type == "poisson"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          fisher += exp_tmp.at(0)*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
        } else if(attr_y_type == "normal"){
          exp_tmp = std::get<0>(pseudo_lh).row(i)*coef;
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t()/y_scale;
          fisher += std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i)/y_scale;
        }
        
      }
      
    }
    coef += (arma::inv(fisher)*score);
    coefs.row(k) = coef.t();
    // If the maximal value of iterations is met end the estimation also if the convergence criteria is met
    // (otherwise start another iteration and reset score and info)
    // TODO there should be a check for identifiabiliy -> popularity
    if(k == max_iteration){
      non_converged = false;
    } else if ((sqrt(sum(arma::pow(coefs.row(k)- coefs.row(k-1),2)))<tol) & !non_stop){
      non_converged = false;
    } else {
      // Reset the score and fisher info
      score.fill(0);
      fisher.fill(0);
    }
    k++;
  }
  return(std::tuple<arma::vec, arma::vec, arma::mat, arma::mat> {coef,score, fisher, coefs.rows(0,k-1)});
}

std::tuple<arma::vec, arma::vec, arma::mat, arma::mat> 
  cond_estimation_nonpopularity_pl(
    arma::vec coef,
    arma::uvec &i_vec,
    arma::uvec &j_vec,
    arma::uvec &overlap_vec, // Assuming this is arma::uvec (0 or 1)
    bool directed,
    std::tuple<arma::mat, arma::vec> &pseudo_lh,
    int max_iteration,
    double tol,
    arma::vec &coef_popularity,
    double offset_nonoverlap,
    bool &non_stop,
    std::string attr_x_type,
    std::string attr_y_type,
    double x_scale,
    double y_scale,
    bool fix_x) {
    
    int n_coef = coef.size();
    unsigned int n_actor = coef_popularity.n_elem / 2;
    unsigned int n_net = i_vec.n_elem;
    
    // Extract the full design matrix and response vector
    const arma::mat& X_all = std::get<0>(pseudo_lh);
    const arma::vec& Y_all = std::get<1>(pseudo_lh);
    
    // --- 1. Define subviews for each component ---
    // Network component is always present
    const arma::mat X_net = X_all.rows(0, n_net - 1);
    const arma::vec Y_net = Y_all.subvec(0, n_net - 1);
    arma::mat X_x, X_y;
    arma::vec Y_x, Y_y;
    
    if (fix_x == false) {
      X_x = X_all.rows(n_net, n_net + n_actor - 1);
      Y_x = Y_all.subvec(n_net, n_net + n_actor - 1);
      
      X_y = X_all.rows(n_net + n_actor, X_all.n_rows - 1);
      Y_y = Y_all.subvec(n_net + n_actor, Y_all.n_elem - 1);
    } else {
      X_y = X_all.rows(n_net, X_all.n_rows - 1);
      Y_y = Y_all.subvec(n_net, Y_all.n_elem - 1);
    }
    
    // Pre-calculate network offsets
    const arma::vec net_offsets = (1.0 - arma::conv_to<arma::vec>::from(overlap_vec)) * offset_nonoverlap;
    
    // Pre-calculate popularity indices
    arma::uvec j_pop_indices; 
    if(directed){
      j_pop_indices = j_vec - 1 + n_actor ;
    } else {
      j_pop_indices =j_vec - 1;
    }
    const arma::uvec i_pop_indices = i_vec - 1;
    
    // --- 2. Initialize estimation variables ---
    arma::vec score(n_coef);
    arma::mat fisher(n_coef, n_coef);
    arma::mat coefs(max_iteration + 1, n_coef);
    coefs.row(0) = coef.t();
    
    bool non_converged = true;
    int k = 1;
    
    // --- 3. Iterative Estimation (Newton-Raphson) ---
    while (non_converged) {
      score.zeros();
      fisher.zeros();
      
      // --- Component 1: Network (Logistic Model) ---
      arma::vec eta_net = X_net * coef + net_offsets + 
        coef_popularity.elem(i_pop_indices) + 
        coef_popularity.elem(j_pop_indices);
      
      arma::vec exp_eta_net = arma::exp(eta_net);
      arma::vec prob_net = exp_eta_net / (1.0 + exp_eta_net);
      arma::vec var_net = prob_net % (1.0 - prob_net);
      
      score += X_net.t() * (Y_net - prob_net);
      fisher += X_net.t() * arma::diagmat(var_net) * X_net;
      
      // --- Component 2: Attribute 'x' ---
      if (fix_x == false) {
        if (attr_x_type == "binomial") {
          arma::vec eta_x = X_x * coef;
          arma::vec exp_eta_x = arma::exp(eta_x);
          arma::vec prob_x = exp_eta_x / (1.0 + exp_eta_x);
          arma::vec var_x = prob_x % (1.0 - prob_x);
          
          score += X_x.t() * (Y_x - prob_x);
          fisher += X_x.t() * arma::diagmat(var_x) * X_x;
          
        } else if (attr_x_type == "poisson") {
          arma::vec eta_x = X_x * coef;
          arma::vec mu_x = arma::exp(eta_x);
          
          score += X_x.t() * (Y_x - mu_x);
          fisher += X_x.t() * arma::diagmat(mu_x) * X_x;
          
        } else if (attr_x_type == "normal") {
          arma::vec mu_x = X_x * coef; 
          score += X_x.t() * (Y_x - mu_x) / x_scale;
          fisher += (X_x.t() * X_x) / x_scale;
        }
      }
      
      // --- Component 3: Attribute 'y' ---
      if (attr_y_type == "binomial") {
        arma::vec eta_y = X_y * coef;
        arma::vec exp_eta_y = arma::exp(eta_y);
        arma::vec prob_y = exp_eta_y / (1.0 + exp_eta_y);
        arma::vec var_y = prob_y % (1.0 - prob_y);
        
        score += X_y.t() * (Y_y - prob_y);
        fisher += X_y.t() * arma::diagmat(var_y) * X_y;
        
      } else if (attr_y_type == "poisson") {
        arma::vec eta_y = X_y * coef;
        arma::vec mu_y = arma::exp(eta_y);
        
        score += X_y.t() * (Y_y - mu_y);
        fisher += X_y.t() * arma::diagmat(mu_y) * X_y;
        
      } else if (attr_y_type == "normal") {
        arma::vec mu_y = X_y * coef;
        
        score += X_y.t() * (Y_y - mu_y) / y_scale;
        fisher += (X_y.t() * X_y) / y_scale;
      }
      // Rcout << "Iteration " << k << score << std::endl;
      // Rcout << "Iteration " << k << fisher << std::endl;
      // --- 4. Update and Check Convergence ---
      coef += (arma::solve(fisher, score));
      coefs.row(k) = coef.t();
      
      if (k == max_iteration) {
        non_converged = false;
      } else if ((arma::norm(coefs.row(k) - coefs.row(k - 1)) < tol) && !non_stop) {
        non_converged = false;
      }
      
      k++;
    }
    
    return std::make_tuple(coef, score, fisher, coefs.rows(0, k - 1));
  }

double calculate_llh(
    const arma::vec& coef,
    arma::vec& coef_popularity,
    const arma::uvec& i_vec,
    const arma::uvec& j_vec,
    const arma::uvec& overlap_vec,
    bool directed,
    const std::tuple<arma::mat, arma::vec>& pseudo_lh,
    double offset_nonoverlap,
    const std::string& attr_x_type,
    const std::string& attr_y_type,
    double x_scale,
    double y_scale,
    int n_actor,
    bool fix_x) {
  // Rcout << "Start" << std::endl;
  
  if(coef_popularity.size() ==1){
    if(directed){
      arma::vec tmp = arma::vec(2*n_actor);
      tmp.fill(0.0);
      coef_popularity = tmp;
    } else {
      arma::vec tmp = arma::vec(n_actor);
      tmp.fill(0.0);
      coef_popularity = tmp;
    }
  }
  unsigned int n_net = i_vec.n_elem;
  // Rcout << coef_popularity.n_elem << std::endl;
  const arma::mat& X_all = std::get<0>(pseudo_lh);
  const arma::vec& Y_all = std::get<1>(pseudo_lh);
  
  const arma::mat X_net = X_all.rows(0, n_net - 1);
  const arma::vec Y_net = Y_all.subvec(0, n_net - 1);
  arma::mat X_x, X_y;
  arma::vec Y_x, Y_y;
  // Rcout << "c" << std::endl;
  if (fix_x == false) {
    X_x = X_all.rows(n_net, n_net + n_actor - 1);
    Y_x = Y_all.subvec(n_net, n_net + n_actor - 1);
    X_y = X_all.rows(n_net + n_actor, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net + n_actor, Y_all.n_elem - 1);
  } else {
    X_y = X_all.rows(n_net, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net, Y_all.n_elem - 1);
  }
  // Rcout << "b" << std::endl;
  const arma::vec net_offsets = (1.0 - arma::conv_to<arma::vec>::from(overlap_vec)) * offset_nonoverlap;
  
  arma::uvec j_pop_indices;
  if(directed){
    j_pop_indices = j_vec - 1 + n_actor;
  } else {
    j_pop_indices = j_vec - 1;
  }
  const arma::uvec i_pop_indices = i_vec - 1;
  
  double llh = 0.0;
  // Rcout << "a" << std::endl;
  // --- Component 1: Network (Logistic Model) ---
  // logL = sum( Y*eta - log(1 + exp(eta)) )
  // Rcout << i_pop_indices.min() << std::endl;
  // Rcout << i_pop_indices.max() << std::endl;
  // 
  // Rcout << coef_popularity.size() << std::endl;
  // Rcout << coef_popularity.elem(i_pop_indices) << std::endl;
  
  arma::vec eta_net = X_net * coef + net_offsets +
    coef_popularity.elem(i_pop_indices) +
    coef_popularity.elem(j_pop_indices);
  
  // Use arma::log1p(exp_eta_net) for the numerically stable log(1 + exp(eta))
  arma::vec exp_eta_net = arma::exp(eta_net);
  // Rcout << mean(exp_eta_net) << std::endl;
  
  llh += arma::sum(Y_net % eta_net - arma::log1p(exp_eta_net));
  // Rcout << "Log-likelihood after network component: " << llh << std::endl;
  // --- Component 2: Attribute 'x' ---
  if (fix_x == false) {
    if (attr_x_type == "binomial") {
      // logL = sum( Y*eta - log(1 + exp(eta)) )
      arma::vec eta_x = X_x * coef;
      arma::vec exp_eta_x = arma::exp(eta_x);
      llh += arma::sum(Y_x % eta_x - arma::log1p(exp_eta_x));
      
    } else if (attr_x_type == "poisson") {
      // logL = sum( Y*eta - exp(eta) - lgamma(Y+1) )
      arma::vec eta_x = X_x * coef;
      arma::vec mu_x = arma::exp(eta_x);
      llh += arma::sum(Y_x % eta_x - mu_x - arma::lgamma(Y_x + 1.0));
      
    } else if (attr_x_type == "normal") {
      // logL = sum( -0.5*log(2*pi*scale) - (Y - mu)^2 / (2*scale) )
      arma::vec mu_x = X_x * coef;
      double const_x = -0.5 * std::log(2.0 * M_PI * x_scale);
      llh += arma::sum(const_x - arma::pow(Y_x - mu_x, 2) / (2.0 * x_scale));
    }
  }
  // Rcout << "Log-likelihood after x component: " << llh << std::endl;
  
  // --- Component 3: Attribute 'y' ---
  if (attr_y_type == "binomial") {
    // logL = sum( Y*eta - log(1 + exp(eta)) )
    arma::vec eta_y = X_y * coef;
    arma::vec exp_eta_y = arma::exp(eta_y);
    llh += arma::sum(Y_y % eta_y - arma::log1p(exp_eta_y));
    
  } else if (attr_y_type == "poisson") {
    // logL = sum( Y*eta - exp(eta) - lgamma(Y+1) )
    arma::vec eta_y = X_y * coef;
    arma::vec mu_y = arma::exp(eta_y);
    llh += arma::sum(Y_y % eta_y - mu_y - arma::lgamma(Y_y + 1.0));
    
  } else if (attr_y_type == "normal") {
    // logL = sum( -0.5*log(2*pi*scale) - (Y - mu)^2 / (2*scale) )
    arma::vec mu_y = X_y * coef;
    double const_y = -0.5 * std::log(2.0 * M_PI * y_scale);
    llh += arma::sum(const_y - arma::pow(Y_y - mu_y, 2) / (2.0 * y_scale));
  }
  // Rcout << "Log-likelihood after y component: " << llh << std::endl;
  
  return llh;
}


arma::mat get_C_new(
    const arma::vec& coef,
    const arma::uvec& i_vec,
    const arma::uvec& j_vec,
    const arma::uvec& overlap_vec,
    bool directed,
    std::tuple<arma::mat, arma::vec>& pseudo_lh, 
    const arma::vec& coef_popularity,
    double offset_nonoverlap,
    bool fix_x,
    const std::string& attr_x_type,
    const std::string& attr_y_type,
    double attr_x_scale,
    double attr_y_scale)
{
  unsigned int n_actor;
  if (directed) {
    n_actor = coef_popularity.n_elem / 2;
  } else { 
    n_actor = coef_popularity.n_elem;
  }
  // number of network dyads
  const unsigned int n_net = i_vec.n_elem;

  // full design matrix and response (pseudo_lh)
  const arma::mat& X_all = std::get<0>(pseudo_lh);
  const arma::vec& Y_all = std::get<1>(pseudo_lh);
  const arma::mat X_net = X_all.rows(0, n_net - 1);

  arma::mat X_x; arma::vec Y_x;
  arma::mat X_y; arma::vec Y_y;

  if (!fix_x) {
    X_x = X_all.rows(n_net, n_net + n_actor - 1);
    Y_x = Y_all.subvec(n_net, n_net + n_actor - 1);
  
    X_y = X_all.rows(n_net + n_actor, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net + n_actor, Y_all.n_elem - 1);
  } else { 
    X_y = X_all.rows(n_net, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net, Y_all.n_elem - 1);
  }
  // 1) network block, which is logistic
  arma::vec net_offsets = (1.0 - arma::conv_to<arma::vec>::from(overlap_vec)) * offset_nonoverlap;
  arma::vec eta_net = X_net * coef + net_offsets + coef_popularity.elem(i_vec-1) + coef_popularity.elem(j_vec-1+ n_actor*directed);
  arma::vec exp_eta_net = arma::exp(eta_net);
  arma::vec prob_net(eta_net.n_elem);
  
  // Find indices for positive and non-positive eta
  arma::uvec pos_idx = arma::find(eta_net > 0);
  arma::uvec nonpos_idx = arma::find(eta_net <= 0);
  
  if (pos_idx.n_elem > 0) {
    prob_net.elem(pos_idx) = 1.0 / (1.0 + arma::exp(-eta_net.elem(pos_idx)));
  }
  if (nonpos_idx.n_elem > 0) {
    arma::vec exp_eta_neg = arma::exp(eta_net.elem(nonpos_idx));
    prob_net.elem(nonpos_idx) = exp_eta_neg / (1.0 + exp_eta_neg);
  }
  arma::vec var_net = prob_net % (1.0 - prob_net);
  
  // Prepare w vector of length X_all.n_rows
  arma::vec w(X_all.n_rows, arma::fill::zeros);

  // assign network block
  w.rows(0, n_net - 1) = var_net;

  // Rcout << "Min variance network block: " << min(exp_eta_net) << std::endl;
  // Rcout << "Max variance network block: " << max(exp_eta_net) << std::endl;
  // Rcout << "Mean variance network block: " << mean(var_net) << std::endl;
  // Rcout << "Mean per col: " << arma::mean(X_net, 0) << std::endl;
  
  
  // 2) attribute x (if present)
  if (!fix_x) {
    arma::vec var_x;
    if (attr_x_type == "binomial") {
      arma::vec eta_x = X_x * coef;
      arma::vec ex = arma::exp(eta_x);
      arma::vec p = ex / (1.0 + ex);
      var_x = p % (1.0 - p);
    } else if (attr_x_type == "poisson") { 
      arma::vec eta_x = X_x * coef;
      arma::vec mu = arma::exp(eta_x);
      var_x = mu;               // variance = mu
    } else if (attr_x_type == "normal" || attr_x_type == "gaussian") { 
      // For normal, variance is constant = attr_x_scale
      var_x = arma::vec(X_x.n_rows, arma::fill::value(attr_x_scale));
    } else { 
      Rcpp::stop("Unknown attr_x_type: must be 'binomial', 'poisson' or 'normal'.");
    } 
    // place into w
    w.rows(n_net, n_net + X_x.n_rows - 1) = var_x;
    // Rcout << "Mean variance x block: " << arma::mean(var_x) << std::endl;
  } 
  
  // 3) attribute y
  arma::vec var_y;
  if (attr_y_type == "binomial") {
    arma::vec eta_y = X_y * coef;
    arma::vec ey = arma::exp(eta_y);
    arma::vec p = ey / (1.0 + ey);
    var_y = p % (1.0 - p);
  } else if (attr_y_type == "poisson") { 
    arma::vec eta_y = X_y * coef;
    arma::vec mu = arma::exp(eta_y);
    var_y = mu;
  } else if (attr_y_type == "normal" || attr_y_type == "gaussian") { 
    var_y = arma::vec(X_y.n_rows, arma::fill::value(attr_y_scale));
  } else { 
    Rcpp::stop("Unknown attr_y_type: must be 'binomial', 'poisson' or 'normal'.");
  } 
  arma::uword start_y;
  if (!fix_x) start_y = n_net + X_x.n_rows;
  else start_y = n_net;
  w.rows(start_y, start_y + X_y.n_rows - 1) = var_y;
  // Rcout << "Mean variance y block: " << arma::mean(var_y) << std::endl;
  // --- Form weighted design and compute Fisher C = X' W X ---
  // Build X_weighted = X_all each column % w
  arma::mat X_weighted = X_all.each_col() % w;
  arma::mat fisher = X_all.t() * X_weighted;
  // Rcpp::Rcout << mean(w) << std::endl;
  return fisher; 
}

arma::mat get_C(arma::vec coef, arma::uvec &i_vec,
                 arma::uvec  &j_vec,
                 arma::uvec  &overlap_vec,
                 bool directed,
                 std::tuple<arma::mat,arma::vec> &pseudo_lh,
                 arma::vec &coef_popularity,
                 double offset_nonoverlap,
                 const std::string& attr_x_type,
                 const std::string& attr_y_type,
                 double attr_x_scale,
                 double attr_y_scale) {
  unsigned int n_actor;
  if(directed){
    n_actor = coef_popularity.n_elem/2;
  } else {
    n_actor = coef_popularity.n_elem;
  }
  // The idea is that we first go over the connections  and then the attributes 
  unsigned int number_elements_network = i_vec.size();
  
  arma::vec vec_network = arma::exp(offset_nonoverlap*(1-overlap_vec)+std::get<0>(pseudo_lh).rows(0,number_elements_network-1)*coef +
    coef_popularity.elem(i_vec-1) + coef_popularity.elem(j_vec-1+ n_actor*directed));
  arma::vec vec_attribute = arma::exp(std::get<0>(pseudo_lh).rows(number_elements_network,std::get<0>(pseudo_lh).n_rows-1)*coef); 
  arma::vec result = join_cols(vec_network, vec_attribute); 
  arma::vec prob = result/(1+result);
  arma::vec values(i_vec.size());
  values.fill(1.0);
  
  const arma::mat& X = std::get<0>(pseudo_lh);
  const arma::vec& w = prob % (1 - prob);
  // Rcpp::Rcout << "Here" << std::endl;
  arma::mat X_weighted = X.each_col() % w;
  // arma::mat cov_popularity_weighted = cov_popularity.each_col() % w;
  // Rcpp::Rcout << "Here" << std::endl;
  arma::mat fisher_alt = X.t() * X_weighted;
  // Rcpp::Rcout << mean(w) << std::endl;
  return fisher_alt;
  
}

std::tuple<arma::vec,arma::vec, arma::mat, arma::mat>  cond_estimation_popularity_pl(arma::vec coef, 
                                                                                     arma::uvec &i_vec, 
                                                                                     arma::uvec  &j_vec,
                                                                                     arma::uvec  &overlap_vec,
                                                                                     bool directed,
                                                                                     std::tuple<arma::mat,arma::vec> &pseudo_lh, 
                                                                                     int max_iteration, 
                                                                                     double tol, 
                                                                                     arma::vec coef_nonpopularity, 
                                                                                     double offset_nonoverlap, 
                                                                                     arma::mat A_inv, 
                                                                                     int it, 
                                                                                     bool &non_stop) {
  int n_actor;
  if(directed){
    n_actor = coef.n_elem/2;
  } else { 
    n_actor = coef.n_elem;
  }
  
  arma::vec  score;
  arma::mat  coefs;
  if(directed){
    score.reshape(2*n_actor, 1);
    score.fill(0);
    coefs.reshape(max_iteration+1,2*n_actor);
  } else {
    score.reshape(n_actor, 1);
    score.fill(0);
    coefs.reshape(max_iteration+1,n_actor);
  }
  arma::vec weights, exp_tmp; 
  coefs.row(0)= coef.t();
  bool non_converged = true;
  int k = 1;
  while(non_converged) {
    for(int i: seq(0,i_vec.size() -1)){
      if(directed){
        if(overlap_vec.at(i)){
          exp_tmp = arma::exp((coef.at(i_vec.at(i)-1) +
            coef.at(j_vec.at(i)-1 + n_actor)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
        } else {
          exp_tmp = arma::exp(offset_nonoverlap+(coef.at(i_vec.at(i)-1) +
            coef.at(j_vec.at(i)-1 + n_actor)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);  
        }
        weights = exp_tmp/(1+exp_tmp);
        if(j_vec.at(i) != n_actor){
          score.at(j_vec.at(i)-1+ n_actor) += std::get<1>(pseudo_lh).at(i) - weights.at(0);
        }
      } else {
        if(overlap_vec.at(i)){
          exp_tmp = arma::exp( (coef.at(i_vec.at(i)-1) + 
            coef.at(j_vec.at(i)-1)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
        } else {
          exp_tmp = arma::exp(offset_nonoverlap+(coef.at(i_vec.at(i)-1) + 
            coef.at(j_vec.at(i)-1)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
        }
        weights = exp_tmp/(1+exp_tmp);
        score.at(j_vec.at(i)-1) += std::get<1>(pseudo_lh).at(i) - weights.at(0);
        
      }
      score.at(i_vec.at(i)-1) +=  std::get<1>(pseudo_lh).at(i) - weights.at(0);
    }
    coef += (A_inv*score);
    coefs.row(k) = coef.t();
    // If the maximal value of iterations is met end the estimation also if the convergence criteria is met 
    // (otherwise start another iteration and reset score and info)
    if(k == max_iteration){
      non_converged = false;
    } else if ((sqrt(sum(arma::pow(coefs.row(k)- coefs.row(k-1),2)))<tol) & !non_stop){
      non_converged = false;
    } else { 
      // Reset the score and fisher info
      // score_alt.fill(0);
      score.fill(0);
    } 
    k++;
  } 
  
  arma::vec coef_MM = coef + A_inv*score;
  // arma::vec exp_tmp_MM,  exp_rest; 
  // // TODO update the calculation of the llh 
  // double ll_MM,  ll_attributes; 
  // if(directed){
  //   exp_tmp_MM = arma::exp(offset_nonoverlap+(coef_MM.elem(i_vec  -1) +
  //     coef_MM.elem(j_vec  -1 + n_actor)) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity);
  //   Rcout << mean(exp_tmp_MM) << std::endl;
  //   exp_rest = arma::exp(std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity);
  //   arma::vec log_one_min_pi_MM =  - log(1+exp_tmp_MM);
  //   arma::vec log_one_min_pi_rest =  - log(1+exp_rest);
  //   
  //   ll_MM = sum(log_one_min_pi_MM + std::get<1>(pseudo_lh).head_rows(i_vec.size())%(offset_nonoverlap+coef_MM.elem(i_vec-1) +
  //     coef_MM.elem(j_vec  -1 + n_actor) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity));
  //   ll_attributes =  sum(log_one_min_pi_rest + std::get<1>(pseudo_lh).tail_rows(2*n_actor)%
  //     (std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity));
  // } else {
  //   exp_tmp_MM = arma::exp(offset_nonoverlap+(coef_MM.elem(i_vec  -1) +
  //     coef_MM.elem(j_vec  -1)) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity);
  //   // Rcout << mean(exp_tmp_MM) << std::endl;
  //   
  //   exp_rest = arma::exp(std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity);
  //   arma::vec log_one_min_pi_MM =  - log(1+exp_tmp_MM);
  //   arma::vec log_one_min_pi_rest =  - log(1+exp_rest);
  //   
  //   ll_MM = sum(log_one_min_pi_MM + std::get<1>(pseudo_lh).head_rows(i_vec.size())%(offset_nonoverlap+coef_MM.elem(i_vec-1) +
  //     coef_MM.elem(j_vec  -1) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity));
  //   ll_attributes =  sum(log_one_min_pi_rest + std::get<1>(pseudo_lh).tail_rows(2*n_actor)%
  //     (std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity));
  //   
  // }
  // double llh_alt; 
  return(std::tuple<arma::vec, arma::vec, arma::mat, arma::mat> {coef,score, A_inv, coefs.rows(0,k-1)});
}


std::tuple<arma::vec,arma::vec, arma::mat , arma::mat>  cond_estimation_popularity_pl_accelerated(arma::vec coef, 
                                                                                                  arma::uvec &i_vec, 
                                                                                                  arma::uvec  &j_vec,
                                                                                                  arma::uvec  &overlap_vec,
                                                                                                  bool directed,
                                                                                                  std::tuple<arma::mat,arma::vec> &pseudo_lh, 
                                                                                                  int max_iteration, 
                                                                                                  double tol, 
                                                                                                  arma::vec coef_nonpopularity, 
                                                                                                  double offset_nonoverlap, 
                                                                                                  arma::mat A_inv, 
                                                                                                  bool &non_stop, 
                                                                                                  arma::vec & old_score_pop, 
                                                                                                  arma::vec & old_coef_pop, 
                                                                                                  arma::mat & old_M, 
                                                                                                  int it, 
                                                                                                  bool first_it) {
  int n_actor;
  if(directed){
    n_actor = coef.n_elem/2;
  } else {
    n_actor = coef.n_elem;
  }
  
  // Define the additional stuff needed for the quasi Newton acceleration
  arma::mat M_new;
  
  arma::vec  score, score_old;
  arma::mat  coefs;
  if(directed){
    score.reshape(2*n_actor, 1);
    score.fill(0);
    score_old.reshape(2*n_actor, 1);
    score_old.fill(0);
    coefs.reshape(max_iteration+1,2*n_actor);
  } else {
    score.reshape(n_actor, 1);
    score.fill(0);
    score_old.reshape(n_actor, 1);
    score_old.fill(0);
    coefs.reshape(max_iteration+1,n_actor);
  }
  arma::vec weights,weights_old,  exp_tmp, exp_tmp_old; 
  coefs.row(0)= coef.t();
  if(first_it){
    old_coef_pop = coef;
  }
  // Rcout << "Hi"<< std::endl;
  bool non_converged = true;
  int k = 1;
  while(non_converged) {
    for(int i: seq(0,i_vec.size() -1)){
      if(directed){
        if(overlap_vec.at(i)){
          exp_tmp = arma::exp((coef.at(i_vec.at(i)-1) +
            coef.at(j_vec.at(i)-1 + n_actor)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
          exp_tmp_old = arma::exp((old_coef_pop.at(i_vec.at(i)-1) +
            old_coef_pop.at(j_vec.at(i)-1 + n_actor)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
          
        } else {
          exp_tmp = arma::exp(offset_nonoverlap+(coef.at(i_vec.at(i)-1) +
            coef.at(j_vec.at(i)-1 + n_actor)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);  
          exp_tmp_old =  arma::exp(offset_nonoverlap+(old_coef_pop.at(i_vec.at(i)-1) +
            old_coef_pop.at(j_vec.at(i)-1 + n_actor)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);  
        }
        weights = exp_tmp/(1+exp_tmp);
        weights_old = exp_tmp_old/(1+exp_tmp_old);
        if(j_vec.at(i) != n_actor){
          score.at(j_vec.at(i)-1+ n_actor) += std::get<1>(pseudo_lh).at(i) - weights.at(0);
          score_old.at(j_vec.at(i)-1+ n_actor) += std::get<1>(pseudo_lh).at(i) - weights_old.at(0);  
        }
      } else {
        if(overlap_vec.at(i)){
          exp_tmp = arma::exp( (coef.at(i_vec.at(i)-1) + 
            coef.at(j_vec.at(i)-1)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
          exp_tmp_old = arma::exp( (old_coef_pop.at(i_vec.at(i)-1) + 
            old_coef_pop.at(j_vec.at(i)-1)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
        } else { 
          exp_tmp = arma::exp(offset_nonoverlap+(coef.at(i_vec.at(i)-1) + 
            coef.at(j_vec.at(i)-1)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
          exp_tmp_old =  arma::exp(offset_nonoverlap+(old_coef_pop.at(i_vec.at(i)-1) + 
            old_coef_pop.at(j_vec.at(i)-1)) + std::get<0>(pseudo_lh).row(i)*coef_nonpopularity);
        } 
        
        weights = exp_tmp/(1+exp_tmp);
        weights_old = exp_tmp_old/(1+exp_tmp_old);
        score.at(j_vec.at(i)-1) += std::get<1>(pseudo_lh).at(i) - weights.at(0);
        score_old.at(j_vec.at(i)-1) += std::get<1>(pseudo_lh).at(i) - weights_old.at(0);
      } 
      score.at(i_vec.at(i)-1) +=  std::get<1>(pseudo_lh).at(i) - weights.at(0);
      score_old.at(i_vec.at(i)-1) +=  std::get<1>(pseudo_lh).at(i) - weights_old.at(0);
    } 
    
    // Rcout <<score.at(2*n_actor-1) << std::endl;
    // Rcout <<coef.at(2*n_actor-1) << std::endl;
    
    if(first_it){
      M_new = old_M;
      old_score_pop = score;
      old_coef_pop = coef;
    } else {
      // Rcout <<sum(old_score_pop)<< std::endl;
      // Rcout <<sum(score_old)<< std::endl;
      // Rcout <<sum(score)<< std::endl;
      // 
      arma::vec score_change_pop = score - score_old;
      old_score_pop = score;
      arma::vec coef_change_pop = coef - old_coef_pop;
      old_coef_pop = coef;
      
      arma::vec r_new = coef_change_pop  + A_inv*score_change_pop;
      arma::vec q_new = r_new - old_M*score_change_pop;
      arma::vec c_new_alt = q_new.t()*score_change_pop;
      M_new = old_M + q_new*q_new.t()/c_new_alt.at(0);
      old_M = M_new;
    }
    
    arma::vec coef_MM = coef + A_inv*score;
    arma::vec coef_accel = coef + (A_inv-M_new)*score;
    arma::vec exp_tmp_MM, exp_tmp_accel; 
    double ll_MM, ll_accel; 
    if(directed){
      exp_tmp_MM = arma::exp(offset_nonoverlap+(coef_MM.elem(i_vec  -1) +
        coef_MM.elem(j_vec  -1 + n_actor)) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity);
      exp_tmp_accel = arma::exp(offset_nonoverlap+(coef_accel.elem(i_vec  -1) +
        coef_accel.elem(j_vec  -1 + n_actor)) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity);
      // exp_rest = arma::exp(std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity);
      arma::vec log_one_min_pi_MM =  - log(1+exp_tmp_MM);
      arma::vec log_one_min_pi_accel =  - log(1+exp_tmp_accel);
      // arma::vec log_one_min_pi_rest =  - log(1+exp_rest);
      
      ll_MM = sum(log_one_min_pi_MM + std::get<1>(pseudo_lh).head_rows(i_vec.size())%(offset_nonoverlap+coef_MM.elem(i_vec-1) +
        coef_MM.elem(j_vec  -1 + n_actor) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity));
      ll_accel = sum(log_one_min_pi_accel + std::get<1>(pseudo_lh).head_rows(i_vec.size())%(offset_nonoverlap+ coef_accel.elem(i_vec-1) +
        coef_accel.elem(j_vec  -1 + n_actor) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity));
      // ll_attributes =  sum(log_one_min_pi_rest + std::get<1>(pseudo_lh).tail_rows(2*n_actor)%
      //   (std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity));
    } else {
      exp_tmp_MM = arma::exp(offset_nonoverlap+(coef_MM.elem(i_vec  -1) +
        coef_MM.elem(j_vec  -1)) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity);
      exp_tmp_accel = arma::exp(offset_nonoverlap+(coef_accel.elem(i_vec  -1) +
        coef_accel.elem(j_vec  -1)) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity);
      // exp_rest = arma::exp(std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity);
      arma::vec log_one_min_pi_MM =  - log(1+exp_tmp_MM);
      arma::vec log_one_min_pi_accel =  - log(1+exp_tmp_accel);
      // arma::vec log_one_min_pi_rest =  - log(1+exp_rest);
      
      ll_MM = sum(log_one_min_pi_MM + std::get<1>(pseudo_lh).head_rows(i_vec.size())%(offset_nonoverlap+coef_MM.elem(i_vec-1) +
        coef_MM.elem(j_vec  -1) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity));
      ll_accel = sum(log_one_min_pi_accel + std::get<1>(pseudo_lh).head_rows(i_vec.size())%(offset_nonoverlap+ coef_accel.elem(i_vec-1) +
        coef_accel.elem(j_vec  -1) + std::get<0>(pseudo_lh).head_rows(i_vec.size())*coef_nonpopularity));
      // ll_attributes =  sum(log_one_min_pi_rest + std::get<1>(pseudo_lh).tail_rows(2*n_actor)%
      //   (std::get<0>(pseudo_lh).tail_rows(2*n_actor)*coef_nonpopularity));
      
    }
    
    if(ll_accel>ll_MM){
      coef = coef_accel;
      
    } else {
      coef = coef_MM;
      
      
    }
    // coef += (A_inv-M_new)*score;
    coefs.row(k) = coef.t();
    // If the maximal value of iterations is met end the estimation also if the convergence criteria is met 
    // (otherwise start another iteration and reset score and info)
    // TODO there should be a check for identifiabiliy -> popularity 
    if(k == max_iteration){
      non_converged = false;
    } else if ((sqrt(sum(arma::pow(coefs.row(k)- coefs.row(k-1),2)))<tol) & !non_stop){ 
      non_converged = false;
    } else {  
      // Reset the score and fisher info
      // score_alt.fill(0);
      score.fill(0);
    }  
    k++;
  }  
  
  return(std::tuple<arma::vec, arma::vec, arma::mat, arma::mat> {coef,score, M_new, coefs.rows(0,k-1)});
}
// [[Rcpp::export]]
List pl_estimation(arma::vec coef,
                   arma::mat z_network ,
                   arma::vec x_attribute ,
                   arma::vec y_attribute ,
                   arma::mat neighborhood,
                   arma::mat overlap,
                   bool directed,
                   std::vector<std::string> terms,
                   std::vector<arma::mat> &data_list,
                   std::vector<double> &type_list, 
                   bool display_progress, 
                   int max_iteration, 
                   double tol, 
                   double offset_nonoverlap, 
                   bool non_stop, 
                   bool fix_x, 
                   std::string attr_x_type, 
                   std::string attr_y_type, 
                   double attr_x_scale, 
                   double attr_y_scale) {
  List res; 
  std::tuple<arma::mat,arma::vec> pseudo_lh;
  int n_actor = y_attribute.size();
  arma::uvec  i_vec, j_vec, overlap_vec; 
  if(directed){
    i_vec = arma::uvec(n_actor*(n_actor-1)); 
    j_vec = arma::uvec(n_actor*(n_actor-1)); 
    overlap_vec = arma::uvec(n_actor*(n_actor-1)); 
  } else { 
    i_vec= arma::uvec(n_actor*(n_actor-1)/2); 
    j_vec= arma::uvec(n_actor*(n_actor-1)/2); 
    overlap_vec = arma::uvec(n_actor*(n_actor-1)/2); 
  } 
  // Calculates the data in a suitable format -> a vector or 32 x p 
  // (being the dimension of the sufficient statistics) matrices corresponding to the data of each dyad
  XYZ_class object(n_actor,directed, x_attribute, y_attribute,z_network,neighborhood,overlap, attr_x_type, attr_y_type,attr_x_scale, attr_y_scale);
  
  int k = 1;
  bool non_converged = true;
  
  
  if(display_progress) {
    Rcout << "Starting with the preprocessing" << std::endl;
  }
  
  arma::vec exp_tmp, score_tmp; 
  pseudo_lh = xyz_get_info_pl(object,terms,data_list,type_list, display_progress,
                              i_vec,j_vec, overlap_vec, n_actor, fix_x);
  
  arma::uvec where_wrong = find(arma::var(std::get<0>(pseudo_lh), 0) == 0);
  if(where_wrong.size() >0){
    Rcout << "Some statistics do not change over all paris/actors (they are excluded from the model since their MLE is negative infinity)" << std::endl;
    arma::uvec where_right = find(arma::var(std::get<0>(pseudo_lh), 0) != 0);
    std::get<0>(pseudo_lh) = std::get<0>(pseudo_lh).cols(where_right);
    coef = coef.rows(where_right);
  }
  int n_coef = coef.size();
  unsigned int n_net = i_vec.n_elem;
  
  arma::mat coefs(max_iteration+1, coef.size()), fisher(coef.size(),coef.size(), arma::fill::zeros);
  k++;
  // Extract the full design matrix and response vector
  const arma::mat& X_all = std::get<0>(pseudo_lh);
  const arma::vec& Y_all = std::get<1>(pseudo_lh);
  
  // --- 1. Define subviews for each component ---
  // Network component is always present
  const arma::mat X_net = X_all.rows(0, n_net - 1);
  const arma::vec Y_net = Y_all.subvec(0, n_net - 1);
  arma::mat X_x, X_y;
  arma::vec Y_x, Y_y;
  
  if (fix_x == false) {
    X_x = X_all.rows(n_net, n_net + n_actor - 1);
    Y_x = Y_all.subvec(n_net, n_net + n_actor - 1);
    
    X_y = X_all.rows(n_net + n_actor, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net + n_actor, Y_all.n_elem - 1);
  } else {
    X_y = X_all.rows(n_net, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net, Y_all.n_elem - 1);
  }
  
  // Pre-calculate network offsets
  const arma::vec net_offsets = (1.0 - arma::conv_to<arma::vec>::from(overlap_vec)) * offset_nonoverlap;
  
  // --- 2. Initialize estimation variables ---
  arma::vec score(n_coef);
  coefs.row(0) = coef.t();
  arma::vec llhs(max_iteration);
  arma::vec stand_in(1);
  stand_in.fill(0.0);
  
  if(display_progress) {
    Rcout << "Starting with the estimation" << std::endl;
  }
  
  // --- 3. Iterative Estimation (Newton-Raphson) ---
  while (non_converged) {
    if(display_progress){
      Rcout << "Iteration " << k - 1 << "\r";
      Rcpp::Rcout.flush();
            
    }
    score.zeros();
    fisher.zeros();
    Rcpp::checkUserInterrupt();
    // --- Component 1: Network (Logistic Model) ---
    arma::vec eta_net = X_net * coef + net_offsets;
    
    arma::vec exp_eta_net = arma::exp(eta_net);
    arma::vec prob_net = exp_eta_net / (1.0 + exp_eta_net);
    arma::vec var_net = prob_net % (1.0 - prob_net);
    
    score += X_net.t() * (Y_net - prob_net);
    fisher += X_net.t() * arma::diagmat(var_net) * X_net;
    
    // --- Component 2: Attribute 'x' ---
    if (fix_x == false) {
      if (attr_x_type == "binomial") {
        arma::vec eta_x = X_x * coef;
        arma::vec exp_eta_x = arma::exp(eta_x);
        arma::vec prob_x = exp_eta_x / (1.0 + exp_eta_x);
        arma::vec var_x = prob_x % (1.0 - prob_x);
        
        score += X_x.t() * (Y_x - prob_x);
        fisher += X_x.t() * arma::diagmat(var_x) * X_x;
        
      } else if (attr_x_type == "poisson") {
        arma::vec eta_x = X_x * coef;
        arma::vec mu_x = arma::exp(eta_x);
        
        score += X_x.t() * (Y_x - mu_x);
        fisher += X_x.t() * arma::diagmat(mu_x) * X_x;
        
      } else if (attr_x_type == "normal") {
        arma::vec mu_x = X_x * coef; 
        score += X_x.t() * (Y_x - mu_x) / attr_x_scale;
        fisher += (X_x.t() * X_x) / attr_x_scale;
      }
    }
    
    // --- Component 3: Attribute 'y' ---
    if (attr_y_type == "binomial") {
      arma::vec eta_y = X_y * coef;
      arma::vec exp_eta_y = arma::exp(eta_y);
      arma::vec prob_y = exp_eta_y / (1.0 + exp_eta_y);
      arma::vec var_y = prob_y % (1.0 - prob_y);
      
      score += X_y.t() * (Y_y - prob_y);
      fisher += X_y.t() * arma::diagmat(var_y) * X_y;
      
    } else if (attr_y_type == "poisson") {
      arma::vec eta_y = X_y * coef;
      arma::vec mu_y = arma::exp(eta_y);
      
      score += X_y.t() * (Y_y - mu_y);
      fisher += X_y.t() * arma::diagmat(mu_y) * X_y;
      
    } else if (attr_y_type == "normal") {
      arma::vec mu_y = X_y * coef;
      
      score += X_y.t() * (Y_y - mu_y) / attr_y_scale;
      fisher += (X_y.t() * X_y) / attr_y_scale;
    }
    // --- 4. Update and Check Convergence ---
    coef += (arma::solve(fisher, score));
    // Rcout << coef <<  std::endl;
    // Rcout << coefs.n_rows <<  std::endl;
    // Rcout << k <<  std::endl;
    // Rcout << coefs.row(k) <<  std::endl;
    // Rcout << "Here" <<  std::endl;
    coefs.row(k) = coef.t();
    // Rcout << "Here" <<  std::endl;
    // Rcout << stand_in <<  std::endl;
    
    llhs.at(k-2) = calculate_llh(coef, 
            stand_in,
            i_vec,
            j_vec,
            overlap_vec,
            directed,
            pseudo_lh,
            offset_nonoverlap,
            attr_x_type,
            attr_y_type,
            attr_x_scale,
            attr_y_scale,
            n_actor,
            fix_x);
    // Rcout << "Here A" <<  std::endl;
    if (k == max_iteration) {
      non_converged = false;
    } else if ((arma::norm(coefs.row(k) - coefs.row(k - 1)) < tol) && !non_stop) {
      non_converged = false;
    }
    
    k++;
  }
  if(display_progress) {
    Rcpp::Rcout.flush();  
    Rcout << "Done with the estimation" << std::endl;
  }
  // coef.rows(ind_popularity) = coef_popularity;
  // coef.rows(ind_nonpopularity) = coef_nonpopularity;
  return(List::create(_["coefficients"] =coef,
                      _["coefficients_path"] =coefs.rows(2,k-1),
                      _["score"] =score,
                      _["where_wrong"] = where_wrong,
                      _["fisher"] = fisher,
                      _["var"] = arma::inv(fisher), 
                      _["llh"] = llhs.head(k-2)
  ));
}

// [[Rcpp::export]]
List pl_estimation_old(arma::vec coef,
                       arma::mat z_network ,
                       arma::vec x_attribute ,
                       arma::vec y_attribute ,
                       arma::mat neighborhood,
                       arma::mat overlap,
                       bool directed,
                       std::vector<std::string> terms,
                       std::vector<arma::mat> &data_list,
                       std::vector<double> &type_list, 
                       bool display_progress, 
                       int max_iteration, 
                       double tol, 
                       double offset_nonoverlap, 
                       bool non_stop, 
                       bool fix_x, 
                       std::string attr_x_type, 
                       std::string attr_y_type, 
                       double attr_x_scale, 
                       double attr_y_scale) {
  List res; 
  std::tuple<arma::mat,arma::vec> pseudo_lh;
  int n_actor = y_attribute.size();
  arma::uvec  i_vec, j_vec, overlap_vec; 
  if(directed){
    i_vec = arma::uvec(n_actor*(n_actor-1)); 
    j_vec = arma::uvec(n_actor*(n_actor-1)); 
    overlap_vec = arma::uvec(n_actor*(n_actor-1)); 
  } else { 
    i_vec= arma::uvec(n_actor*(n_actor-1)/2); 
    j_vec= arma::uvec(n_actor*(n_actor-1)/2); 
    overlap_vec = arma::uvec(n_actor*(n_actor-1)/2); 
  } 
  
  // Calculates the data in a suitable format -> a vector or 32 x p 
  // (being the dimension of the sufficient statistics) matrices corresponding to the data of each dyad
  XYZ_class object(n_actor,directed, x_attribute, y_attribute,z_network,neighborhood,overlap, attr_x_type, attr_y_type,attr_x_scale, attr_y_scale);
  
  int k = 1;
  bool non_converged = true;
  
  
  arma::vec exp_tmp, score_tmp; 
  pseudo_lh = xyz_get_info_pl(object,terms,data_list,type_list, display_progress,
                              i_vec,j_vec, overlap_vec, n_actor, fix_x);
  
  arma::uvec where_wrong = find(arma::var(std::get<0>(pseudo_lh), 0) == 0);
  if(where_wrong.size() >0){
    Rcout << "Some statistics do not change over all paris/actors (they are excluded from the model since their MLE is negative infinity)" << std::endl;
    arma::uvec where_right = find(arma::var(std::get<0>(pseudo_lh), 0) != 0);
    std::get<0>(pseudo_lh) = std::get<0>(pseudo_lh).cols(where_right);
    coef = coef.rows(where_right);
  }
  arma::mat coefs(max_iteration+1, coef.size()), fisher(coef.size(),coef.size(), arma::fill::zeros);
  coefs.row(0)= coef.t();
  arma::vec score(coef.size(), arma::fill::zeros);
  coefs.row(0)= coef.t();
  score.fill(0);
  fisher.fill(0);
  unsigned int number_elements_network = i_vec.size();
  
  while(non_converged) {
    if(display_progress) {
      Rcout << "Iteration = " + std::to_string(k)  << std::endl;
    }
    for(unsigned int i = 0; i < std::get<0>(pseudo_lh).n_rows; i++){
      // What to do if we are regarding network information (relating to the first entries)
      if(i < number_elements_network){
        if(overlap_vec.at(i)){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
        } else {
          exp_tmp = arma::exp(offset_nonoverlap+std::get<0>(pseudo_lh).row(i)*coef);
        } 
        score_tmp = std::get<0>(pseudo_lh).row(i).t()*
          std::get<1>(pseudo_lh).at(i) -
          exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
        score += score_tmp;
        fisher += exp_tmp.at(0)/(pow(1+exp_tmp.at(0), 2))*
          std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
        
      } else if((i > number_elements_network) & (i < number_elements_network +n_actor) & !fix_x){ 
        // Rcpp::Rcout << "X" << std::endl;
        
        // What to do if we are regarding attribute information (relating to the last entries)
        if(attr_x_type == "binomial"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score_tmp = std::get<0>(pseudo_lh).row(i).t()*
            std::get<1>(pseudo_lh).at(i) -
            exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          score += score_tmp;
          fisher += exp_tmp.at(0)/(pow(1+exp_tmp.at(0), 2))*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
          
        } 
        if(attr_x_type == "poisson"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          Rcpp::Rcout << "Here" << std::endl;
          
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          Rcpp::Rcout << "Here" << std::endl;
          
          fisher += exp_tmp.at(0)*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
        } else if(attr_x_type == "normal"){ 
          exp_tmp = std::get<0>(pseudo_lh).row(i)*coef;
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp)*std::get<0>(pseudo_lh).row(i).t()/attr_x_scale;
          fisher += std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i)/attr_x_scale;
        } 
        
        
      } else {
        // Rcpp::Rcout << "Y" << std::endl;
        
        if(attr_y_type == "binomial"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score_tmp = std::get<0>(pseudo_lh).row(i).t()*
            std::get<1>(pseudo_lh).at(i) -
            exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          score += score_tmp;
          fisher += exp_tmp.at(0)/(pow(1+exp_tmp.at(0), 2))*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
          
        } 
        if(attr_y_type == "poisson"){
          exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
          fisher += exp_tmp.at(0)*
            std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
        } else if(attr_y_type == "normal"){ 
          exp_tmp = std::get<0>(pseudo_lh).row(i)*coef;
          score += (std::get<1>(pseudo_lh).at(i) - exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t()/attr_y_scale;
          fisher += std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i)/attr_y_scale;
        } 
        
      }
      
    }
    coef += (arma::inv(fisher)*score);
    coefs.row(k) = coef.t();
    if(k == max_iteration){
      non_converged = false;
    } else if ((sqrt(sum(arma::pow(coefs.row(k)- coefs.row(k-1),2)))<tol) & !non_stop){
      non_converged = false;
    } else {
      score.fill(0);
      fisher.fill(0);
    }
    k++;
  }
  
  // coef.rows(ind_popularity) = coef_popularity;
  // coef.rows(ind_nonpopularity) = coef_nonpopularity;
  return(List::create(_["coefficients"] =coef,
                      _["coefficients_path"] =coefs.rows(0,k-1),
                      _["score"] =score,
                      _["fisher"] = fisher,
                      _["var"] = arma::inv(fisher)
  ));
}



// // Here I am planning on an alternative approach to solve the maximization problem via iterative least squares
// // [[Rcpp::export]]
// List pl_estimation_ils(arma::vec coef,
//                        arma::mat z_network ,
//                        arma::vec x_attribute ,
//                        arma::vec y_attribute ,
//                        arma::mat neighborhood,
//                        arma::mat overlap,
//                        bool directed,
//                        std::vector<std::string> terms,
//                        std::vector<arma::mat> &data_list,
//                        std::vector<double> &type_list, 
//                        bool display_progress, 
//                        int max_iteration, 
//                        double tol, 
//                        double offset_nonoverlap, 
//                        bool non_stop) {
//   List res; 
//   std::tuple<arma::mat,arma::vec> pseudo_lh;
//   int n_actor = y_attribute.size();
//   arma::mat coefs(max_iteration+1, coef.size()), fisher(coef.size(),coef.size(), arma::fill::zeros);
//   coefs.row(0)= coef.t();
//   arma::vec score(coef.size(), arma::fill::zeros), i_vec, j_vec,overlap_vec;
//   if(directed){
//     // network_vec = arma::vec(n_actor*(n_actor-1)); 
//     i_vec = arma::vec(n_actor*(n_actor-1)); 
//     j_vec = arma::vec(n_actor*(n_actor-1)); 
//     overlap_vec = arma::vec(n_actor*(n_actor-1)); 
//   } else {  
//     // network_vec= arma::vec(n_actor*(n_actor-1)/2); 
//     i_vec= arma::vec(n_actor*(n_actor-1)/2); 
//     j_vec= arma::vec(n_actor*(n_actor-1)/2);
//     overlap_vec= arma::vec(n_actor*(n_actor-1)/2);
//   }  
//   
//   // Calculates the data in a suitable format -> a vector or 32 x p 
//   // (being the dimension of the sufficient statistics) matrices corresponding to the data of each dyad
//   
//   // Rcout << "Start preparation"<< std::endl;
//   
//   // Rcout << "Here we are"<< std::endl;
//   XYZ_class object(n_actor,directed, x_attribute, y_attribute,z_network,neighborhood,overlap);
//   // Rcout << "Network generated"<< std::endl;
//   pseudo_lh = xyz_get_info_pl(object,terms,data_list,type_list, display_progress,
//                               i_vec,j_vec,overlap_vec, n_actor);
//   int k = 1;
//   bool non_converged = true;
//   
//   coefs.row(0)= coef.t();
//   score.fill(0);
//   fisher.fill(0);
//   arma::vec exp_tmp, score_tmp; 
//   arma::vec weights; 
//   arma::vec target; 
//   while(non_converged) {
//     // Rcout << "Iteration = " + std::to_string(k)  << std::endl;
//     if(display_progress) {
//       Rcout << "Iteration = " + std::to_string(k)  << std::endl;
//     } 
//     for(unsigned int i = 0; i < std::get<0>(pseudo_lh).n_rows; i++){
//       exp_tmp = arma::exp(std::get<0>(pseudo_lh).row(i)*coef);
//       score_tmp = std::get<0>(pseudo_lh).row(i).t()*std::get<1>(pseudo_lh).at(i) - 
//         exp_tmp.at(0)/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t();
//       score += score_tmp;
//       fisher += exp_tmp.at(0)/(1+exp_tmp.at(0))*1/(1+exp_tmp.at(0))*std::get<0>(pseudo_lh).row(i).t()*std::get<0>(pseudo_lh).row(i);
//     }  
//     coef += (arma::inv(fisher)*score);
//     weights = arma::exp( std::get<0>(pseudo_lh)*coef);
//     weights = weights/(1+weights);; 
//     
//     target = std::get<0>(pseudo_lh); 
//     coefs.row(k) = coef.t();
//     if(k == max_iteration){
//       non_converged = false;
//     } else if ((sqrt(sum(arma::pow(coefs.row(k)- coefs.row(k-1),2)))<tol) & !non_stop){ 
//       non_converged = false;
//     } else { 
//       score.fill(0);
//       fisher.fill(0);
//     } 
//     k++;
//     // Rcout << "Check done"<< std::endl;
//   }
//   // coef.rows(ind_popularity) = coef_popularity;
//   // coef.rows(ind_nonpopularity) = coef_nonpopularity;
//   return(List::create(_["coefficients"] =coef,
//                       _["coefficients_path"] =coefs.rows(0,k-1),
//                       _["score"] =score,
//                       _["fisher"] = fisher,
//                       _["var"] = arma::inv(fisher)
//   ));
// }

// [[Rcpp::export]]
arma::mat invert_mat(double diag, double offdiag,int n_actor){
  return(1/(diag-offdiag)*arma::mat(n_actor,n_actor, arma::fill::eye) - arma::mat(n_actor, n_actor, arma::fill::value(1/((1/offdiag + n_actor/(diag-offdiag))*pow(diag-offdiag,2)))));
}

// [[Rcpp::export]]
arma::mat get_A_inv(double n_actor){
  arma::mat A = (n_actor-1)/4*arma::mat(n_actor,n_actor, arma::fill::eye);
  arma::mat B = arma::mat(n_actor,n_actor-1, arma::fill::value(0.25));
  B.diag() = arma::vec(n_actor-1, arma::fill::zeros);
  
  arma::mat D = (n_actor-1)/4*arma::mat(n_actor-1,n_actor-1, arma::fill::eye);
  arma::mat A_inv = 1/A.at(0,0)*arma::mat(n_actor,n_actor, arma::fill::eye);
  arma::mat D_inv = 1/D.at(0,0)*arma::mat(n_actor-1,n_actor-1, arma::fill::eye);
  arma::mat part_2 = arma::mat(n_actor-1,n_actor-1, arma::fill::value((n_actor-2)/(n_actor-1)*0.25));
  part_2.diag() = arma::vec(n_actor-1, arma::fill::value(0.25));
  part_2 = D - part_2;
  arma::mat part_1 = arma::mat(n_actor,n_actor, arma::fill::value((n_actor-3)/(n_actor-1)*0.25));
  part_1.diag() = arma::vec(n_actor, arma::fill::value((n_actor-2)/(n_actor-1)*0.25));
  part_1.row(n_actor-1) = arma::vec(n_actor, arma::fill::value((n_actor-2)/(n_actor-1)*0.25)).as_row();
  part_1.col(n_actor-1) = arma::vec(n_actor, arma::fill::value((n_actor-2)/(n_actor-1)*0.25)).as_col();
  part_1.at(n_actor-1, n_actor-1) = 0.25;
  part_1 = A - part_1;
  
  // arma::mat part_2 = D - B.t()*A_inv*B;
  // arma::mat part_1 = A - B*D_inv*B.t();
  
  arma::mat a = part_1.submat(0,0,n_actor-2,n_actor-2);
  arma::mat b = part_1.submat(n_actor-1,0,n_actor-1,n_actor-2);
  arma::mat c = b.t();
  arma::mat d = arma::mat(1,1, arma::fill::zeros);
  d.at(0,0) = part_1.at(n_actor-1,n_actor-1);
  
  
  arma::mat Part_2 = d-b*invert_mat(a.at(1,1), a.at(0,1), a.n_rows)*b.t();
  arma::mat Part_2_invert = arma::mat(1,1, arma::fill::zeros);
  Part_2_invert.at(0,0) = 1/Part_2.at(0,0);
  arma::mat Part_1 = a-b.t()*b/d.at(0,0);
  arma::mat Part_1_invert = invert_mat(Part_1.at(1,1), Part_1.at(0,1), Part_1.n_rows);
  arma::mat Part_3 = -b*1/d.at(0,0);
  arma::mat Part_4 = -c.t()*invert_mat(a.at(1,1), a.at(0,1), a.n_rows);
  arma::mat inv_part_1 = arma::join_rows(arma::join_cols( Part_1_invert, (Part_1_invert*Part_3.t()).t()), 
                                         arma::join_cols( (Part_2_invert.at(0,0)*Part_4).t(), Part_2_invert));
  arma::mat inv_part_2 = invert_mat(part_2.at(1,1), part_2.at(0,1), part_2.n_rows);
  arma::mat part_3 = -arma::mat(n_actor, n_actor-1, arma::fill::value(1/(n_actor-1)));
  part_3.diag() = arma::vec(n_actor-1, arma::fill::zeros); 
  // arma::mat part_3 = -B*D_inv;
  // arma::mat tmp_mat = (inv_part_1* part_3); 
  // Rcout << arma::accu(part_3)<< std::endl;
  // Rcout << "Diagonal"<< std::endl;
  // Rcout << part_3.at(0,1)*(inv_part_1.at(1,2)*(n_actor-2)+inv_part_1.at(1,n_actor-1))<< std::endl;
  // Rcout << "Off-Diagonal"<< std::endl;
  // Rcout << part_3.at(0,1)*(inv_part_1.at(0,1)*(n_actor-3)+
  //   inv_part_1.at(0,0)+inv_part_1.at(0,n_actor-1))<< std::endl;
  // Rcout << "Last Row"<< std::endl;
  // Rcout << part_3.at(0,1)*(inv_part_1.at(0,n_actor-1)*(n_actor-2)+inv_part_1.at(n_actor-1,n_actor-1))<< std::endl;
  
  
  arma::mat tmp_mat = arma::mat(n_actor, n_actor-1,
                                arma::fill::value(part_3.at(0,1)*(inv_part_1.at(0,1)*(n_actor-3)+
                                  inv_part_1.at(0,0)+inv_part_1.at(0,n_actor-1)))); 
  tmp_mat.diag() = arma::vec(n_actor-1, arma::fill::value(part_3.at(0,1)*(inv_part_1.at(1,2)*(n_actor-2)+inv_part_1.at(1,n_actor-1))));
  tmp_mat.row(n_actor-1) = arma::vec(n_actor-1,arma::fill::value(part_3.at(0,1)*(inv_part_1.at(0,n_actor-1)*(n_actor-2)+inv_part_1.at(n_actor-1,n_actor-1)))).as_row();
  // Rcout << tmp_mat<< std::endl;
  
  arma::mat res = arma::join_rows(arma::join_cols( inv_part_1, tmp_mat.t()), 
                                  arma::join_cols( tmp_mat, inv_part_2));
  // Rcout << arma::accu(tmp_mat)<< std::endl;
  res = arma::join_cols(res, arma::mat(1,2*n_actor-1, arma::fill::zeros));
  res = arma::join_rows(res, arma::mat(2*n_actor,1, arma::fill::zeros));
  return(res);
}

// [[Rcpp::export]]
List outerloop_estimation_pl(arma::vec coef,
                             arma::vec coef_popularity,
                             arma::mat z_network ,
                             arma::vec x_attribute ,
                             arma::vec y_attribute ,
                             arma::mat neighborhood,
                             arma::mat overlap,
                             bool directed,
                             std::vector<std::string> terms,
                             std::vector<arma::mat> &data_list,
                             std::vector<double> &type_list, 
                             bool display_progress, 
                             int max_iteration_outer, 
                             int max_iteration_inner_popularity, 
                             int max_iteration_inner_nonpopularity, 
                             double tol, 
                             double offset_nonoverlap, 
                             bool non_stop, 
                             bool var, 
                             bool accelerated, 
                             bool fix_x, 
                             std::string type_x, 
                             std::string type_y, 
                             double attr_x_scale, 
                             double attr_y_scale, 
                             int start = 0) {
  List res; 
  // iglm_print_registered_functions();
  // auto m = iglm::Registry::instance().info("repetition");
  // Rcpp::Rcout << m.fn << "  " << m.value << "\n";
  // 
  // auto n = iglm::Registry::instance().info("edges");
  // Rcpp::Rcout << n.short_name << "  " << n.value << "\n";
  // 
  std::tuple<arma::mat,arma::vec> pseudo_lh;
  int n_actor = y_attribute.size();
  
  arma::mat id_mat = arma::mat((double) n_actor,n_actor,arma::fill::eye);
  arma::mat ones_mat = arma::mat(n_actor,n_actor,arma::fill::ones);
  arma::mat A_inv = 4/((double)n_actor -2)*(id_mat - 1/(2*(double)n_actor-2)*ones_mat);
  if(directed){
    A_inv = get_A_inv(n_actor);
  }
  
  std::tuple<arma::vec, arma::vec,  arma::mat, arma::mat> res_popularity, res_nonpopularity, res_nonpopularity_alt;
  arma::uvec i_vec, j_vec,overlap_vec;
  arma::mat coefs_popularity; 
  //  coefs_nonpopularity(terms.size()), coefs_popularity(n_actor)
  if(directed){
    // network_vec = arma::vec(n_actor*(n_actor-1)); 
    i_vec = arma::uvec(n_actor*(n_actor-1)); 
    j_vec = arma::uvec(n_actor*(n_actor-1)); 
    overlap_vec = arma::uvec(n_actor*(n_actor-1)); 
    coef_popularity.reshape(n_actor*2,1);
    coefs_popularity.reshape(0, n_actor*2);  
  } else {
    // network_vec= arma::vec(n_actor*(n_actor-1)/2); 
    i_vec= arma::uvec(n_actor*(n_actor-1)/2); 
    j_vec= arma::uvec(n_actor*(n_actor-1)/2);
    overlap_vec= arma::uvec(n_actor*(n_actor-1)/2);
    coef_popularity.reshape(n_actor,1);
    coefs_popularity.reshape(max_iteration_inner_popularity, n_actor);
  }
  XYZ_class object(n_actor,directed, x_attribute, y_attribute,z_network,neighborhood,overlap, type_x, type_y,attr_x_scale, attr_y_scale);
  // Rcout << "Start"<< std::endl;
  // Rcout << object.n_actor<< std::endl;
  
  if(display_progress) {
    Rcout << "Starting with the preprocessing" << std::endl;
  }
  pseudo_lh = xyz_get_info_pl(object,terms,data_list,type_list, display_progress,
                              i_vec,j_vec,overlap_vec, n_actor, fix_x);
  arma::uvec where_wrong = find(arma::sum(std::get<0>(pseudo_lh), 0) == 0);
  if(where_wrong.size() >0){
    Rcout << "Some statistics do not change over all paris/actors (they are excluded from the model since their MLE is negative infinity)" << std::endl;
    arma::uvec where_right = find(arma::var(std::get<0>(pseudo_lh), 0) != 0);
    std::get<0>(pseudo_lh) = std::get<0>(pseudo_lh).cols(where_right);
    coef = coef.rows(where_right);
  }
  arma::mat coefs(max_iteration_outer+1, coef.size() + coef_popularity.size()),
  fisher_popularity(terms.size(),terms.size()), fisher_nonpopularity(terms.size(),terms.size()), 
  coefs_nonpopularity(max_iteration_inner_nonpopularity, terms.size());
  
  coefs.row(0)= arma::join_rows(coef.t(), coef_popularity.t());
  arma::vec score_popularity(n_actor), 
  score_nonpopularity(terms.size()),coef_nonpopularity(terms.size());
  
  int k = 1;
  bool non_converged = true;
  bool first_it = true;
  
  coef_nonpopularity = coef;
  arma::vec old_coef_pop, old_score_pop, llh(max_iteration_outer +1);
  arma::mat old_M;
  llh.at(0) = 0; 
  
  if(accelerated){
    old_M.reshape(coef_popularity.size(), coef_popularity.size());
    old_M.fill(arma::fill::zeros); 
    old_coef_pop.reshape(coef_popularity.size(),1);
    old_coef_pop.fill(arma::fill::zeros);
  }
  if(display_progress) {
    Rcout << "Starting with the estimation" << std::endl;
  }
  
  while(non_converged) {
    Rcpp::checkUserInterrupt();
    
    if(display_progress) {
      Rcpp::Rcout.flush();
      Rcout << "Iteration = " + std::to_string(k+start)  << "\r";
    }
    if(accelerated){
      res_popularity = cond_estimation_popularity_pl_accelerated(coef_popularity,
                                                                 i_vec, 
                                                                 j_vec,
                                                                 overlap_vec,
                                                                 directed,
                                                                 pseudo_lh, 
                                                                 max_iteration_inner_popularity, 
                                                                 tol, 
                                                                 coef_nonpopularity, 
                                                                 offset_nonoverlap, 
                                                                 A_inv, 
                                                                 non_stop, 
                                                                 old_score_pop,
                                                                 old_coef_pop, 
                                                                 old_M,  k, first_it);
      first_it = false;  
    } else {
      res_popularity = cond_estimation_popularity_pl(coef_popularity,
                                                     i_vec, 
                                                     j_vec,
                                                     overlap_vec,
                                                     directed,
                                                     pseudo_lh, 
                                                     max_iteration_inner_popularity, 
                                                     tol, 
                                                     coef_nonpopularity, 
                                                     offset_nonoverlap, 
                                                     A_inv, k, non_stop);
      
    }
    coef_popularity = std::get<0>(res_popularity);
    // res_nonpopularity_alt = cond_estimation_nonpopularity_pl_old(coef_nonpopularity, 
    //                                                      i_vec, 
    //                                                      j_vec,
    //                                                      overlap_vec,
    //                                                      directed,
    //                                                      pseudo_lh, 
    //                                                      max_iteration_inner_nonpopularity, 
    //                                                      tol, 
    //                                                      coef_popularity, 
    //                                                      offset_nonoverlap, 
    //                                                      non_stop, 
    //                                                      type_x, type_y,  
    //                                                      attr_x_scale, attr_y_scale, fix_x);
    // Rcout << std::get<0>(res_nonpopularity_alt)<< std::endl;
    res_nonpopularity = cond_estimation_nonpopularity_pl(coef_nonpopularity, 
                                                         i_vec, 
                                                         j_vec,
                                                         overlap_vec,
                                                         directed,
                                                         pseudo_lh, 
                                                         max_iteration_inner_nonpopularity, 
                                                         tol, 
                                                         coef_popularity, 
                                                         offset_nonoverlap, 
                                                         non_stop, 
                                                         type_x, type_y, 
                                                         attr_x_scale, attr_y_scale, fix_x);
    // Rcout << std::get<0>(res_nonpopularity)<< std::endl;
    // Rcout << "Done 2. Stage"<< std::endl;
    coef_nonpopularity = std::get<0>(res_nonpopularity);
    llh.at(k) = calculate_llh(coef_nonpopularity, 
           coef_popularity,
           i_vec,
           j_vec,
           overlap_vec,
           directed,
           pseudo_lh,
           offset_nonoverlap,
           type_x,
           type_y,
           attr_x_scale,
           attr_y_scale,
           n_actor,
           fix_x);
    // Rcout << coef_nonpopularity<< std::endl;
    coefs.row(k) = join_cols(coef_nonpopularity, coef_popularity).t();
    if(k == max_iteration_outer){
      non_converged = false;
    } else if ((arma::max(arma::vec{arma::norm(coefs.row(k)- coefs.row(k-1), 2), 
                          std::abs((llh.at(k) -llh.at(k-1))/llh.at(k))})<tol) & !non_stop){
      non_converged = false;
    }
    k++;
    // Rcout << "Check done"<< std::endl;
  }
  if(display_progress) {
    Rcpp::Rcout.flush();  
    Rcout << "Done with the estimation" << std::endl;
  }
  
  std::tie(coef_nonpopularity,score_nonpopularity,fisher_nonpopularity,coefs_nonpopularity) = res_nonpopularity;
  std::tie(coef_popularity,score_popularity,fisher_popularity,coefs_popularity) = res_popularity;
  
  if(var){
    std::tuple<arma::vec,  arma::mat> res_mat = get_B_pl(i_vec, 
                                                         j_vec,overlap_vec,
                                                         pseudo_lh, 
                                                         coef_popularity, 
                                                         coef_nonpopularity, offset_nonoverlap, 
                                                         object.z_network.directed,object.n_actor );
    arma::mat exact_A = get_A_exact(i_vec, 
                                    j_vec,overlap_vec,
                                    pseudo_lh, 
                                    coef_popularity, 
                                    coef_nonpopularity, offset_nonoverlap, 
                                    object.z_network.directed,object.n_actor );
    arma::mat B_mat;
    arma::vec A_diag;
    std::tie(A_diag,B_mat) = res_mat;
    // coefs(ind_popularity) = coef_popularity;
    // coef.rows(ind_nonpopularity) = coef_nonpopularity;
    if(accelerated){
      return(List::create(_["coefficients_nonpopularity"] =coef_nonpopularity,
                          _["coefficients_popularity"] =coef_popularity,
                          _["coefficients_path"] =coefs.rows(0,k-1),
                          _["score_popularity"] =score_popularity,
                          _["score_nonpopularity"] =score_nonpopularity,
                          _["fisher_popularity"] = fisher_popularity, 
                          _["fisher_nonpopularity"] = fisher_nonpopularity, 
                          _["A_diag"] = A_diag, 
                          _["M"] = old_M, 
                          _["A_approx"] = A_inv, 
                          _["exact_A"] =exact_A,
                          _["B_mat"] = B_mat, 
                          _["llh"] = llh.rows(1,k-1), 
                          _["where_wrong"] = where_wrong));
    } else {
      return(List::create(_["coefficients_nonpopularity"] =coef_nonpopularity,
                          _["coefficients_popularity"] =coef_popularity,
                          _["coefficients_path"] =coefs.rows(0,k-1),
                          _["score_popularity"] =score_popularity,
                          _["score_nonpopularity"] =score_nonpopularity,
                          _["fisher_popularity"] = fisher_popularity, 
                          _["fisher_nonpopularity"] = fisher_nonpopularity, 
                          _["A_diag"] = A_diag, 
                          _["B_mat"] = B_mat, 
                          _["llh"] = llh.rows(1,k-1), 
                          _["where_wrong"] = where_wrong));
    }
  } else {
    // coefs(ind_popularity) = coef_popularity;
    // coef.rows(ind_nonpopularity) = coef_nonpopularity;
    return(List::create(_["coefficients_nonpopularity"] =coef_nonpopularity,
                        _["coefficients_popularity"] =coef_popularity,
                        _["coefficients_path"] =coefs.rows(0,k-1),
                        _["score_popularity"] =score_popularity,
                        _["score_nonpopularity"] =score_nonpopularity,
                        _["fisher_popularity"] = fisher_popularity, 
                        _["fisher_nonpopularity"] = fisher_nonpopularity, 
                        _["A_inv"] = A_inv,
                        _["llh"] = llh.rows(1,k-1), 
                        _["where_wrong"] = where_wrong
    ));
  }
  
  
  
}



// std::vector<arma::mat> xyz_prepare_composite_estimation_internal_approx(XYZ_class object,
//                                                                         std::vector<std::string> terms,
//                                                                         std::vector<arma::mat> &data_list,
//                                                                         std::vector<double> &type_list, 
//                                                                         bool add_info, double prob, int seed) {
//   // Set up objects
//   int n_actor = object.n_actor;
//   bool is_full_neighborhood = object.check_if_full_neighborhood();
//   
//   std::vector<xyz_ValidateFunction> functions;
//   functions = xyz_change_statistics_generate(terms);
//   std::string z = "z", x = "x", y = "y";
//   arma::vec change_stat_x_i(functions.size()),
//   change_stat_x_j(functions.size()), change_stat_y_i(functions.size()),
//   change_stat_y_j(functions.size()), change_stat_z_ij(functions.size());
//   bool x_i, x_j, z_ij, y_i,y_j;
//   std::vector<arma::mat> res(n_actor*(n_actor-1)/2);
//   // int now = 0;
//   set_seed(seed);
//   arma::vec change_stat;
//   NumericVector random_accept= runif(n_actor*(n_actor-1),0,1);
//   // Number of trials
//   int n_trials = 0; 
//   // Number of success
//   int n_success = 0; 
//   
//   for(int i: seq(1,n_actor-1)){
//     for(int j: seq(i+1,n_actor)){
//       if(random_accept.at(n_trials)>prob){
//         n_trials++;
//         continue;
//       }
//       n_trials++;
//       
//       // Get present values of x_ij, y_i, y_j
//       x_i = object.x_attribute.get_val(i);
//       x_j = object.x_attribute.get_val(j);
//       y_i = object.y_attribute.get_val(i);
//       y_j = object.y_attribute.get_val(j);
//       z_ij = object.z_network.get_val(i,j);
//       
//       if(add_info) {
//         res.at(n_success) = get_all_combinations_xyz(i,j,x_i, x_j,y_i,y_j,z_ij,
//                n_success,  object,functions, data_list, type_list, is_full_neighborhood);
//       } else { 
//         res.at(n_success) = get_all_combinations_xyz(i,j,x_i, x_j,y_i,y_j,z_ij,
//                n_success,  object,functions, data_list, type_list, is_full_neighborhood).submat(0,8,31,7+functions.size());
//       } 
//       // offset_new.at(n_success) = offset.at()
//       n_success ++;
//     } 
//   }
//   return(res);
// } 

// arma::vec calculate_score_approx(XYZ_class object,
//                                  arma::vec coef,
//                                  std::vector<std::string> terms,
//                                  std::vector<arma::mat> &data_list,
//                                  std::vector<double> &type_list, 
//                                  double prob, int seed) {
//   double w_tmp;
//   std::vector<arma::mat> composite_lh;
//   // Rcout << "A"<< std::endl;
//   
//   composite_lh = xyz_prepare_composite_estimation_internal_approx(object,
//                                                                   terms,
//                                                                   data_list,
//                                                                   type_list,
//                                                                   false, prob, seed);
//   int n_coef = terms.size();
//   arma::vec score(n_coef), score_tmp(n_coef), m(n_coef), exp_tmp(n_coef);
//   // // Rcout << "Iteration = " + std::to_string(k)  << std::endl;
//   // Update the score and fisher info
//   int end = 0;
//   for(unsigned int i = 0; i < composite_lh.size(); i++){
//     if(composite_lh.at(i).size()== 0) {
//       continue;
//     }  
//     end ++;
//     // Rcout << "Iteration = " + std::to_string(i - composite_lh.size())  << std::endl;
//     exp_tmp = arma::exp(composite_lh.at(i)*coef);
//     // Rcout << exp_tmp  << std::endl;
//     // Rcout << "B"  << std::endl;
//     w_tmp = sum(exp_tmp);
//     // Rcout << w_tmp  << std::endl;
//     // Rcout << "C"  << std::endl;
//     m = arma::sum(composite_lh.at(i).t()*exp_tmp, 1);
//     // Rcout << m  << std::endl;
//     // Rcout << "D"  << std::endl;
//     score_tmp = m/w_tmp;
//     // Rcout << score_tmp  << std::endl;
//     // Rcout << "E"  << std::endl;
//     score -= score_tmp;
//   }
//   // Rcout << composite_lh.size()  << std::endl;
//   // Rcout << end  << std::endl;
//   score = score*composite_lh.size()/end;
//   return(score);
// }


arma::vec calculate_score_pl(XYZ_class & object,
                             arma::vec coef,
                             std::vector<std::string> terms,
                             std::vector<arma::mat> &data_list,
                             std::vector<double> &type_list, 
                             double &offset_nonoverlap, 
                             bool fix_x, 
                             std::string attr_x_type, 
                             std::string attr_y_type, 
                             double attr_x_scale, 
                             double attr_y_scale) {
  // double w_tmp;
  arma::uvec i_vec, j_vec,overlap_vec;
  if(object.z_network.directed){
    // network_vec = arma::vec(n_actor*(n_actor-1)); 
    i_vec = arma::uvec(object.n_actor*(object.n_actor-1)); 
    j_vec = arma::uvec(object.n_actor*(object.n_actor-1)); 
    overlap_vec = arma::uvec(object.n_actor*(object.n_actor-1)); 
  } else { 
    // network_vec= arma::vec(n_actor*(n_actor-1)/2); 
    i_vec= arma::uvec(object.n_actor*(object.n_actor-1)/2); 
    j_vec= arma::uvec(object.n_actor*(object.n_actor-1)/2); 
    overlap_vec= arma::uvec(object.n_actor*(object.n_actor-1)/2); 
  } 
  
  std::tuple<arma::mat, arma::vec> pseudo_lh = xyz_get_info_pl(object,terms,data_list,
                                                               type_list, false,
                                                               i_vec,j_vec, overlap_vec,
                                                               object.n_actor, fix_x);
  
  int n_actor = object.n_actor, n_coef = coef.size();
  unsigned int n_net = i_vec.n_elem;
  
  const arma::mat& X_all = std::get<0>(pseudo_lh);
  const arma::vec& Y_all = std::get<1>(pseudo_lh);
  
  // --- 1. Define subviews for each component ---
  // Network component is always present
  const arma::mat X_net = X_all.rows(0, n_net - 1);
  const arma::vec Y_net = Y_all.subvec(0, n_net - 1);
  arma::mat X_x, X_y;
  arma::vec Y_x, Y_y;
  
  if (fix_x == false) {
    X_x = X_all.rows(n_net, n_net + n_actor - 1);
    Y_x = Y_all.subvec(n_net, n_net + n_actor - 1);
    
    X_y = X_all.rows(n_net + n_actor, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net + n_actor, Y_all.n_elem - 1);
  } else {
    X_y = X_all.rows(n_net, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net, Y_all.n_elem - 1);
  }
  
  // Pre-calculate network offsets
  const arma::vec net_offsets = (1.0 - arma::conv_to<arma::vec>::from(overlap_vec)) * offset_nonoverlap;
  // --- 2. Initialize estimation variables ---
  arma::vec score(n_coef);
  score.zeros();
  
  // --- Component 1: Network (Logistic Model) ---
  arma::vec eta_net = X_net * coef + net_offsets;
  
  arma::vec exp_eta_net = arma::exp(eta_net);
  arma::vec prob_net = exp_eta_net / (1.0 + exp_eta_net);
  arma::vec var_net = prob_net % (1.0 - prob_net);
  
  score += X_net.t() * (Y_net - prob_net);
  // --- Component 2: Attribute 'x' ---
  if (fix_x == false) {
    if (attr_x_type == "binomial") {
      arma::vec eta_x = X_x * coef;
      arma::vec exp_eta_x = arma::exp(eta_x);
      arma::vec prob_x = exp_eta_x / (1.0 + exp_eta_x);
      arma::vec var_x = prob_x % (1.0 - prob_x);
      score += X_x.t() * (Y_x - prob_x);
    } else if (attr_x_type == "poisson") {
      arma::vec eta_x = X_x * coef;
      arma::vec mu_x = arma::exp(eta_x);
      score += X_x.t() * (Y_x - mu_x);
    } else if (attr_x_type == "normal") {
      arma::vec mu_x = X_x * coef; 
      score += X_x.t() * (Y_x - mu_x) / attr_x_scale;
    }
  }
  
  // --- Component 3: Attribute 'y' ---
  if (attr_y_type == "binomial") {
    arma::vec eta_y = X_y * coef;
    arma::vec exp_eta_y = arma::exp(eta_y);
    arma::vec prob_y = exp_eta_y / (1.0 + exp_eta_y);
    arma::vec var_y = prob_y % (1.0 - prob_y);
    score += X_y.t() * (Y_y - prob_y);
  } else if (attr_y_type == "poisson") {
    arma::vec eta_y = X_y * coef;
    arma::vec mu_y = arma::exp(eta_y);
    score += X_y.t() * (Y_y - mu_y);
  } else if (attr_y_type == "normal") {
    arma::vec mu_y = X_y * coef;
    score += X_y.t() * (Y_y - mu_y) / attr_y_scale;
  }
  return(score);
}

arma::vec calculate_score_pl_popularity(XYZ_class & object,
                                        arma::vec coef_nonpopularity,
                                        arma::vec coef_popularity,
                                        std::vector<std::string> terms,
                                        std::vector<arma::mat> &data_list,
                                        std::vector<double> &type_list, 
                                        double &offset_nonoverlap, 
                                        bool fix_x, 
                                        bool updated_uncertainty,
                                        bool exact, 
                                        std::string attr_x_type, 
                                        std::string attr_y_type, 
                                        double attr_x_scale, 
                                        double attr_y_scale) {
  // double w_tmp;
  arma::uvec i_vec, j_vec,overlap_vec;
  if(object.z_network.directed){
    i_vec = arma::uvec(object.n_actor*(object.n_actor-1)); 
    j_vec = arma::uvec(object.n_actor*(object.n_actor-1)); 
    overlap_vec = arma::uvec(object.n_actor*(object.n_actor-1)); 
  } else {  
    i_vec= arma::uvec(object.n_actor*(object.n_actor-1)/2); 
    j_vec= arma::uvec(object.n_actor*(object.n_actor-1)/2); 
    overlap_vec= arma::uvec(object.n_actor*(object.n_actor-1)/2); 
  }  
  std::tuple<arma::mat, arma::vec> pseudo_lh = xyz_get_info_pl(object,terms,data_list,
                                                               type_list, false,
                                                               i_vec,j_vec, overlap_vec,
                                                               object.n_actor, fix_x);
  
  
  int n_coef = coef_nonpopularity.size();
  unsigned int n_net = i_vec.n_elem;
  
  // Extract the full design matrix and response vector
  const arma::mat& X_all = std::get<0>(pseudo_lh);
  const arma::vec& Y_all = std::get<1>(pseudo_lh);
  
  // --- 1. Define subviews for each component ---
  // Network component is always present
  const arma::mat X_net = X_all.rows(0, n_net - 1);
  const arma::vec Y_net = Y_all.subvec(0, n_net - 1);
  arma::mat X_x, X_y;
  arma::vec Y_x, Y_y;
  
  if (fix_x == false) {
    X_x = X_all.rows(n_net, n_net + object.n_actor - 1);
    Y_x = Y_all.subvec(n_net, n_net + object.n_actor - 1);
    
    X_y = X_all.rows(n_net + object.n_actor, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net + object.n_actor, Y_all.n_elem - 1);
  } else {
    X_y = X_all.rows(n_net, X_all.n_rows - 1);
    Y_y = Y_all.subvec(n_net, Y_all.n_elem - 1);
  }
  
  // Pre-calculate network offsets
  const arma::vec net_offsets = (1.0 - arma::conv_to<arma::vec>::from(overlap_vec)) * offset_nonoverlap;
  
  // Pre-calculate popularity indices
  arma::uvec j_pop_indices; 
  if(object.z_network.directed){
    j_pop_indices = j_vec - 1 + object.n_actor ;
  } else {
    j_pop_indices =j_vec - 1;
  }
  const arma::uvec i_pop_indices = i_vec - 1;
  
  // --- 2. Initialize estimation variables ---
  arma::vec score_nonpopularity(n_coef, arma::fill::zeros);
  
  // --- Component 1: Network (Logistic Model) ---
  arma::vec eta_net = X_net * coef_nonpopularity + net_offsets + 
    coef_popularity.elem(i_pop_indices) + 
    coef_popularity.elem(j_pop_indices);
  
  arma::vec exp_eta_net = arma::exp(eta_net);
  arma::vec prob_net = exp_eta_net / (1.0 + exp_eta_net);
  arma::vec var_net = prob_net % (1.0 - prob_net);
  
  score_nonpopularity += X_net.t() * (Y_net - prob_net);
  // --- Component 2: Attribute 'x' ---
  if (fix_x == false) {
    if (attr_x_type == "binomial") {
      arma::vec eta_x = X_x * coef_nonpopularity;
      arma::vec exp_eta_x = arma::exp(eta_x);
      arma::vec prob_x = exp_eta_x / (1.0 + exp_eta_x);
      arma::vec var_x = prob_x % (1.0 - prob_x);
      
      score_nonpopularity += X_x.t() * (Y_x - prob_x);
    } else if (attr_x_type == "poisson") {
      arma::vec eta_x = X_x * coef_nonpopularity;
      arma::vec mu_x = arma::exp(eta_x);
      score_nonpopularity += X_x.t() * (Y_x - mu_x);
    } else if (attr_x_type == "normal") {
      arma::vec mu_x = X_x * coef_nonpopularity; 
      score_nonpopularity += X_x.t() * (Y_x - mu_x) / attr_x_scale;
    }
  }
  
  // --- Component 3: Attribute 'y' ---
  if (attr_y_type == "binomial") {
    arma::vec eta_y = X_y * coef_nonpopularity;
    arma::vec exp_eta_y = arma::exp(eta_y);
    arma::vec prob_y = exp_eta_y / (1.0 + exp_eta_y);
    arma::vec var_y = prob_y % (1.0 - prob_y);
    
    score_nonpopularity += X_y.t() * (Y_y - prob_y);
  } else if (attr_y_type == "poisson") {
    arma::vec eta_y = X_y * coef_nonpopularity;
    arma::vec mu_y = arma::exp(eta_y);
    
    score_nonpopularity += X_y.t() * (Y_y - mu_y);
  } else if (attr_y_type == "normal") {
    arma::vec mu_y = X_y * coef_nonpopularity;
    score_nonpopularity += X_y.t() * (Y_y - mu_y) / attr_y_scale;
  }
  
  arma::vec score_popularity(coef_popularity.size(), arma::fill::zeros);
  for(unsigned int i = 0; i < i_vec.size(); i++){
    score_popularity.at(i_vec.at(i)-1) +=  std::get<1>(pseudo_lh).at(i) - prob_net.at(i);
    if(object.z_network.directed){
      score_popularity.at(j_vec.at(i)-1+ object.n_actor) +=  std::get<1>(pseudo_lh).at(i) - prob_net.at(i);
    } else {
      score_popularity.at(j_vec.at(i)-1) +=  std::get<1>(pseudo_lh).at(i) - prob_net.at(i);  
    }
  }  
  arma::vec res_vec;
  if(updated_uncertainty){
    // arma::mat C = get_C(coef_nonpopularity,i_vec,j_vec,overlap_vec,
    //                        object.z_network.directed,
    //                        pseudo_lh,coef_popularity,
    //                        offset_nonoverlap, 
    //                        attr_x_type, 
    //                        attr_y_type, 
    //                        attr_x_scale, 
    //                        attr_y_scale);
    arma::mat C = get_C_new(coef_nonpopularity,i_vec,j_vec,overlap_vec,
                        object.z_network.directed,
                        pseudo_lh,coef_popularity, 
                        offset_nonoverlap, 
                        fix_x,
                        attr_x_type, 
                        attr_y_type, 
                        attr_x_scale, 
                        attr_y_scale);
    arma::mat A = get_A_exact(i_vec, j_vec,overlap_vec,pseudo_lh, 
                              coef_popularity, 
                              coef_nonpopularity, 
                              offset_nonoverlap,
                              object.z_network.directed, 
                              object.n_actor);
    // clock.tock("A");
    
    // clock.tick("B");
    arma::mat B = get_B(i_vec, j_vec,overlap_vec,
                        pseudo_lh, coef_popularity,coef_nonpopularity, offset_nonoverlap,
                        object.z_network.directed,
                        object.n_actor).t();
    
    arma::mat X;
    // clock.tick("solve");
    if(exact){
      X = arma::pinv(A) * B;
      // Rcout << "Using exact inversion for A" << std::endl;
      arma::mat residual = A * X - B;
      
    } else {
      // Step 1: Inverse of diagonal A
      arma::vec ainv = 1.0 / A.diag();  // Elementwise inverse
      // Step 2: Compute X = A^{-1} B using column-wise scaling
      X = B.each_col() % ainv;
    }
    // Compute Schur complement
    arma::mat S = C - B.t() * X;
    // Invert S
    // Rcout << "Inversion for S" << std::endl;
    // Rcout << C << std::endl;
    // Rcout << C_new << std::endl;
    // Rcout << X << std::endl;
    // Rcout << B << std::endl;
    arma::mat S_inv = arma::pinv(S);
    
    // Compute inverse blocks
    arma::mat M22_inv = S_inv;
    arma::mat M12_inv = -X * S_inv;
    res_vec = join_cols(M12_inv.t()*score_popularity + M22_inv*score_nonpopularity, score_popularity); 
  } else {
    res_vec = join_cols(score_nonpopularity, score_popularity); 
  }
  return(res_vec);
} 
// 
// // [[Rcpp::export]]
// void update_progress_cpp(int total_iterations = 100, int width = 50) {
//   for (int i = 0; i <= total_iterations; ++i) {
//     // 1. Calculate progress
//     double progress = static_cast<double>(i) / total_iterations;
//     int pos = static_cast<int>(width * progress);
//     
//     // 2. Build the progress bar string
//     std::string bar = "[";
//     for (int j = 0; j < width; ++j) {
//       if (j < pos) bar += "=";
//       else if (j == pos) bar += ">";
//       else bar += " ";
//     } 
//     bar += "]";
//      
//     // 3. Print the progress bar with carriage return
//     //    Use Rcout for R console compatibility
//     //    Add spaces at the end to overwrite longer previous lines if necessary
//     Rcpp::Rcout << "Progress: " << static_cast<int>(progress * 100.0) << "% " << bar << "   \r";
//      
//     // 4. Flush the output stream
//     Rcpp::Rcout.flush();
//      
//     // 5. Pause briefly to simulate work (optional)
//     //    usleep takes microseconds (1 million = 1 second)
//     usleep(50000); // Sleep for 50 milliseconds
//   } 
//   
//   // 6. Print a newline at the end to finish cleanly
//   Rcpp::Rcout << std::endl;
// }

// [[Rcpp::export]]
List xyz_approximate_variability(arma::vec& coef,
                                 arma::vec& coef_popularity,
                                 std::vector<std::string>& terms,
                                 int& n_actor,
                                 arma::mat z_network,
                                 arma::mat neighborhood,
                                 arma::mat overlap,
                                 arma::vec y_attribute,
                                 arma::vec x_attribute,
                                 bool init_empty,
                                 bool directed,
                                 std::vector<arma::mat>& data_list,
                                 std::vector<double>& type_list,
                                 int n_proposals_x,
                                 int seed_x,
                                 int n_proposals_y,
                                 int seed_y,
                                 int n_proposals_z,
                                 int seed_z,
                                 int n_burn_in,
                                 int n_simulation,
                                 bool display_progress,
                                 bool popularity,
                                 double offset_nonoverlap, 
                                 bool return_samples, 
                                 bool fix_x, 
                                 bool updated_uncertainty, 
                                 bool exact, 
                                 std::string type_x, 
                                 std::string type_y, 
                                 double attr_x_scale, 
                                 double attr_y_scale){
  // Generate the class with the provided information
  XYZ_class object(n_actor,directed, neighborhood, overlap, type_x, type_y,attr_x_scale, attr_y_scale);
  if(!init_empty){
    object.set_info_arma(x_attribute,y_attribute, z_network);
  }
  bool is_full_neighborhood = object.check_if_full_neighborhood();
  // object.print();
  // Rcout <<is_full_neighborhood<< std::endl;
  
  // object.set_neighborhood_from_mat(neighborhood);
  // Generate change statistic function from the terms
  std::vector<xyz_ValidateFunction> functions;
  
  functions = xyz_change_statistics_generate_new(terms);
  
  // arma::vec at_zero;
  // at_zero = xyz_eval_at_empty_network_new(terms);
  // // Start for a burn in period with the normal number of proposals
  // // Intialize global statistics and then adapt them peu a peu
  // arma::vec global_stats(functions.size());
  // global_stats.fill(0);
  // global_stats = global_stats + at_zero;
  arma::vec global_stats = xyz_count_global_internal( object,
                                                      terms,
                                                      n_actor,
                                                      data_list,
                                                      type_list, 
                                                      type_x, type_y, 
                                                      attr_x_scale, 
                                                      attr_y_scale);
  // Matrix of statistics that hold the simulated statistics
  // (we should probably return those as well, since they might be useful later on)
  arma::mat stats(n_simulation,functions.size());
  // The statistics are filled with zeroes at first
  stats.fill(0);
  // We also need to evaluate the gradients of the composite lh for each simulation
  // These are saved in the matrix "gradients"
  int size_gradient;
  if(popularity){
    size_gradient = coef.size() + coef_popularity.size();
  } else {
    size_gradient = coef.size();
  }
  arma::mat gradients(n_simulation,size_gradient);
  arma::mat current_hession(size_gradient, size_gradient);
  
  gradients.fill(0);
  
  // arma::mat gradients_approx(n_simulation,functions.size());
  // gradients_approx.fill(0);
  
  // Progress p(n_simulation + n_burn_in, display_progress);
  // We generate vectors for the exact values of all simulates of xyz objects
  // (since x and y are attributes, they are represented as sets of integers,
  //  z is a networks and hence an unordered_map of sets)
  std::vector<arma::vec> res_x(n_simulation);
  std::vector<arma::vec> res_y(n_simulation);
  std::vector<std::unordered_map< int, std::unordered_set<int>>> res_z(n_simulation);
  
  arma::vec z_tmp(32); 
  z_tmp.zeros();
  for(int i: seq(0,31)){
    if (i % 2 == 0)
      z_tmp.at(i) = 1;
  }
  // arma::vec gradient_tmp;
  for(int i = 1; i <=(n_simulation+n_burn_in);i ++) {
    // Rcout << " Mean x \n";
    // Rcout << mean(object.x_attribute.attribute) << std::endl;
    // Rcout << " Mean y \n";
    // Rcout << mean(object.y_attribute.attribute) << std::endl;
    // Rcout << " Number z\n";
    // Rcout << object.z_network.number_edges() << std::endl;
    Rcpp::checkUserInterrupt();
    if(display_progress) {
      Rcpp::Rcout.flush();
      Rcout << "Sample: " + std::to_string(i)  << "\r";
    }
    // Simulate the network
    // Rcout << "Start Network"  << std::endl;
    // p.increment(); // update progress
    
    if(!fix_x){
      // Sample X| Y,Z
      
      
      xyz_simulate_attribute_mh(coef,object,
                                n_proposals_x,seed_x +i,
                                data_list,
                                type_list,
                                is_full_neighborhood,
                                functions,
                                global_stats, "x");  
    }
    // Sample Y| X,Z
    // Rcout << "Sampling Y| X,Z" << std::endl;
    
    xyz_simulate_attribute_mh(coef,object,
                              n_proposals_y,seed_y +i,
                              data_list, type_list,
                              is_full_neighborhood, functions,
                              global_stats, "y");
    
    // Sample Z|X,Y
    // Rcout << "Sampling Z|X,Y" << std::endl;
    if(popularity){
      xyz_simulate_network_mh_popularity(coef,
                                         coef_popularity,
                                         object,
                                         n_proposals_z, seed_z +i,
                                         data_list, type_list,
                                         is_full_neighborhood, functions,
                                         global_stats,offset_nonoverlap); 
    } else {
      // Rcout << "Start Network"  << std::endl;
      xyz_simulate_network_mh(coef,object,
                              n_proposals_z, seed_z +i,
                              data_list, type_list,
                              is_full_neighborhood, functions,
                              global_stats,  offset_nonoverlap);  
    }
    // Rcout << "Sampling Z_neigh|X,Y" << std::endl;
    
    
    if(popularity){
      xyz_simulate_network_consecutive_popularity_mh(coef,
                                                     coef_popularity,object, seed_z*2 +i,
                                                     data_list, type_list,
                                                     is_full_neighborhood, functions,
                                                     global_stats, offset_nonoverlap);
    } else {
      // Rcout << "Start Network"  << std::endl;
      xyz_simulate_network_consecutive_mh(coef,object,seed_z*2 +i,
                                          data_list, type_list,
                                          is_full_neighborhood, functions,
                                          global_stats, offset_nonoverlap);
    }
    if(i>n_burn_in){
      if(popularity){
        // Rcout <<  "Trying to calculate the score" << std::endl;
        // Rcout <<  updated_uncertainty << std::endl;
        // Rcout <<  exact << std::endl;
        gradients.row(i - n_burn_in-1) = calculate_score_pl_popularity(object,
                      coef,
                      coef_popularity,
                      terms,
                      data_list,
                      type_list,
                      offset_nonoverlap, fix_x, updated_uncertainty, exact, 
                      type_x, 
                      type_y, 
                      attr_x_scale, 
                      attr_y_scale).as_row();
      } else {
        gradients.row(i - n_burn_in-1) = calculate_score_pl(object,
                      coef,
                      terms,
                      data_list,
                      type_list,
                      offset_nonoverlap, fix_x, 
                      type_x, 
                      type_y, 
                      attr_x_scale, 
                      attr_y_scale).t();
      }
      if(return_samples){
        res_x.at(i - n_burn_in-1) =object.x_attribute.attribute;
        res_z.at(i - n_burn_in-1) =object.z_network.adj_list;  
        res_y.at(i - n_burn_in-1) =object.y_attribute.attribute;  
      }
      // Count global statistics
      stats.row(i - n_burn_in-1) = global_stats.as_row();
      // n ++;
    }
    
  }
  if(popularity){
    arma::uvec ind_nonpopularity, ind_popularity;
    if(directed) {
      ind_nonpopularity = arma::regspace<arma::uvec>(0, terms.size() -1); 
      ind_popularity = arma::regspace<arma::uvec>(terms.size(), terms.size() +n_actor*2-1); 
    } else {
      ind_nonpopularity = arma::regspace<arma::uvec>(0, terms.size() -1); 
      ind_popularity = arma::regspace<arma::uvec>(terms.size(), terms.size() +n_actor-1); 
    }
    
    arma::mat gradients_popularity,gradients_nonpopularity;
    gradients_popularity = gradients.cols(ind_popularity);
    gradients_nonpopularity = gradients.cols(ind_nonpopularity);
    // Rcout << "Here" << std::endl;
    // Rcout << gradients_nonpopularity << std::endl;
    if(return_samples){
      return(List::create(_["simulation_x_attributes"] =res_x,
                          _["simulation_y_attributes"] =res_y,
                          _["simulation_z_networks"] =res_z,
                          _["stats"] = stats, _["gradients"] = gradients, 
                          _["gradients_nonpopularity"] = gradients_nonpopularity,
                          _["gradients_popularity"] = gradients_popularity));  
    } else {
      return(List::create(_["stats"] = stats, _["gradients"] = gradients, 
                          _["gradients_nonpopularity"] = gradients_nonpopularity,
                          _["gradients_popularity"] = gradients_popularity));
    }
    
  } else{
    if(return_samples){
      return(List::create(_["simulation_x_attributes"] =res_x,
                          _["simulation_y_attributes"] =res_y,
                          _["simulation_z_networks"] =res_z,
                          _["stats"] = stats, _["gradients"] = gradients));
    }else {
      return(List::create(_["stats"] = stats, _["gradients"] = gradients));
    }
    
  }
  
}


// [[Rcpp::export]]
arma::mat xyz_prepare_pseudo_estimation(arma::mat z_network,
                                        arma::vec x_attribute,
                                        arma::vec y_attribute ,
                                        arma::mat neighborhood,
                                        arma::mat overlap,
                                        bool directed,
                                        std::vector<std::string> terms,
                                        std::vector<arma::mat> &data_list,
                                        std::vector<double> &type_list, 
                                        bool display_progress, 
                                        std::string type_x,
                                        std::string type_y,
                                        double attr_x_scale,
                                        double attr_y_scale,
                                        bool return_x = false,
                                        bool return_y = false,
                                        bool return_z = false) {
  // Set up objects
  int n_actor = y_attribute.size();
  // Rcout << "Read Data" << std::endl;
  XYZ_class object(n_actor,directed, x_attribute, y_attribute,z_network,neighborhood,overlap, type_x, type_y,attr_x_scale, attr_y_scale);
  // Rcout << "Done" << std::endl;
  // object.x_attribute.print();
  // Check whether its a fully observed neighbhorhood (this means that everyone knows everyone)
  // This is provided to the sufficient statistics as this might make some calculations unnecessary
  bool is_full_neighborhood = object.check_if_full_neighborhood();
  // Generate vector of functions that calculate the sufficient statistics
  std::vector<xyz_ValidateFunction> functions;
  functions = xyz_change_statistics_generate_new(terms);
  // strings z, x, and y later needed to tell the sufficient statistics
  // what type of change statistic is wanted
  std::string z = "z", x = "x", y = "y";
  // Just a temporary vector of the change statistics for one dyad of the network
  arma::vec change_stat_network(functions.size()),
  // The same thing but for the attributes i and j
  change_stat_attribute_i(functions.size());
  int x_i, y_i, z_ij;
  // int ncores = 5;
  arma::mat res_covs, res_target;
  // Rcout << "Here" << std::endl;
  // Rcout << terms.size() << std::endl;
  // Rcout << n_actor << std::endl;
  if(directed){
    res_covs = arma::mat(n_actor*(n_actor-1) + n_actor*2,terms.size());
    res_target = arma::mat(n_actor*(n_actor-1) +n_actor*2,1);
  } else {
    res_covs = arma::mat(n_actor*(n_actor-1)/2 + n_actor*2,terms.size());
    res_target = arma::mat(n_actor*(n_actor-1)/2 +n_actor*2,1);
  }
  // Rcout << object.z_network.number_edges() << std::endl;
  int now = 0;
  if(return_z){
    for(int i: seq(1,n_actor-1)){
      
      
      for(int j: seq(1,n_actor)){
        if(i == j){
          continue;
        }
        if((i > j) & (!directed)){
          continue;
        }
        // Rcout << i  << std::endl;
        
        // Get present values of x_ij, y_i, y_j
        z_ij = object.z_network.get_val(i,j);
        xyz_calculate_change_stats(change_stat_network, i,
                                                         j,
                                                         object,
                                                         data_list,
                                                         type_list,
                                                         z,
                                                         is_full_neighborhood,
                                                         functions);
        res_covs.row(now)= (change_stat_network).as_row();
        res_target.row(now) = z_ij;
        now += 1;
      }
    }
  }
  if(return_x){
    for(int i: seq(1,n_actor)){
      // Rcout << i  << std::endl;
      x_i = object.x_attribute.get_val(i);
      xyz_calculate_change_stats(change_stat_attribute_i, i,
                                                           i,
                                                           object,
                                                           data_list,
                                                           type_list,
                                                           x,
                                                           is_full_neighborhood,
                                                           functions);
      res_covs.row(now)= (change_stat_attribute_i).as_row();
      res_target.row(now) = x_i;
      now += 1;
    }
  }
  
  if(return_y){
    for(int i: seq(1,n_actor)){
      // Rcout << i  << std::endl;
      y_i = object.y_attribute.get_val(i);
      xyz_calculate_change_stats(change_stat_attribute_i, i,
                                                           i,
                                                           object,
                                                           data_list,
                                                           type_list,
                                                           y,
                                                           is_full_neighborhood,
                                                           functions);
      res_covs.row(now)= (change_stat_attribute_i).as_row();
      res_target.row(now) = y_i;
      now += 1;
    }
  }
  arma::mat res;
  // res_covs.insert_cols(res_target);
  res = arma::join_rows(res_target,res_covs);
  return(res.rows(0, now-1));
}

