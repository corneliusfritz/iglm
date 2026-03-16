#include <RcppArmadillo.h>
#include "iglm/attribute_class.h"
#include "iglm/network_class.h"
#include "iglm/xz_class.h"
#include "iglm/xyz_class.h"
#include "iglm/helper_functions.h"

// Attribute implementations
Attribute::Attribute(int a, std::string type_, double scale_) {
    n_actor = a;
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
}

Attribute::Attribute(int a, arma::vec attribute_tmp, std::string type_, double scale_) {
    n_actor = a;
    scale = scale_;
    attribute = attribute_tmp;
    if(type_ == "binomial" || type_ == "poisson" || type_ == "normal"){
        type = type_;   
    } else {
        Rcpp::Rcout << "Invalid type, we assume that it is binomial (only binomial, poisson, and normal are implemented.)\n";
        type = "binomial";
    }
}

bool Attribute::check() const {
    int mycount = std::count_if(attribute.begin(),
                                attribute.end(),
                                IsMoreThan(n_actor));
    return(mycount==0);
}

void Attribute::print() {
    attribute.print();
}

// Network implementations
Network::Network(int n_actor_, bool directed_) {
    n_actor = n_actor_;
    directed = directed_;
    adj_list.resize(n_actor + 1);
    if (directed) {
        adj_list_in.resize(n_actor + 1);
    }
    adj_mat.assign(n_actor * n_actor, 0);
    number_edges = 0;
}

Network::Network(int n_actor_, bool directed_, arma::mat mat) {
    n_actor = n_actor_;
    directed = directed_;
    adj_list.resize(n_actor + 1);
    if (directed) {
        adj_list_in.resize(n_actor + 1);
    }
    adj_mat.assign(n_actor * n_actor, 0);
    mat_to_map_vec(mat, n_actor, directed, adj_list, adj_list_in, adj_mat);
    number_edges = count_edges(); 
}

void Network::set_network_from_mat(int n_actor_, bool directed_, arma::mat mat){
    n_actor = n_actor_;
    directed = directed_;
    adj_list.clear();
    adj_list.resize(n_actor + 1);
    if (directed) {
        adj_list_in.clear();
        adj_list_in.resize(n_actor + 1);
    }
    adj_mat.assign(n_actor * n_actor, 0);
    mat_to_map_vec(mat, n_actor, directed, adj_list, adj_list_in, adj_mat);
    number_edges = count_edges(); 
}

void Network::change_edge(int from, int to) {
    if(adj_mat[get_mat_idx(from, to)]) {
        delete_edge(from, to);
    } else {
        add_edge(from, to);
    }
}

size_t Network::count_common_partners(unsigned int from, unsigned int to, std::string type) const {
    if (type == "OTP") {
        const std::vector<int>& l1 = adj_list[from];
        if (directed) {
            const std::vector<int>& l2 = adj_list_in[to];
             if (l1.size() < l2.size()) {
                size_t count = 0;
                for (int k : l1) if (adj_mat[get_mat_idx(k, to)]) count++;
                return count;
            } else {
                size_t count = 0;
                for (int k : l2) if (adj_mat[get_mat_idx(from, k)]) count++;
                return count;
            }
        }
    } else if (type == "ISP") {
        if (directed) {
            const std::vector<int>& l1 = adj_list_in[from];
            const std::vector<int>& l2 = adj_list_in[to];
            if (l1.size() < l2.size()) {
                size_t count = 0;
                for (int k : l1) if (adj_mat[get_mat_idx(k, to)]) count++;
                return count;
            } else {
                size_t count = 0;
                for (int k : l2) if (adj_mat[get_mat_idx(k, from)]) count++;
                return count;
            }
        }
    } else if (type == "OSP") {
        const std::vector<int>& l1 = adj_list[from];
        const std::vector<int>& l2 = adj_list[to];
        if (l1.size() < l2.size()) {
            size_t count = 0;
            for (int k : l1) if (adj_mat[get_mat_idx(to, k)]) count++;
            return count;
        } else {
            size_t count = 0;
            for (int k : l2) if (adj_mat[get_mat_idx(from, k)]) count++;
            return count;
        }
    } else if (type == "ITP") {
        if (directed) {
            const std::vector<int>& l1 = adj_list_in[from];
            const std::vector<int>& l2 = adj_list[to];
            if (l1.size() < l2.size()) {
                size_t count = 0;
                for (int k : l1) if (adj_mat[get_mat_idx(to, k)]) count++;
                return count;
            } else {
                size_t count = 0;
                for (int k : l2) if (adj_mat[get_mat_idx(k, from)]) count++;
                return count;
            }
        }
    }
    return 0;
}

std::vector<int> Network::get_common_partners(unsigned int from,unsigned int to, std::string type) const {
    if(type == "OTP"){
        return(get_intersection_vec(adj_list[from], adj_list_in[to])); 
    } else  if(type == "ISP"){
        return(get_intersection_vec(adj_list_in[from], adj_list_in[to])); 
    }else  if(type == "OSP"){
        return(get_intersection_vec(adj_list[from], adj_list[to])); 
    }
    else  if(type == "ITP"){
        return(get_intersection_vec(adj_list_in[from], adj_list[to])); 
    }
    return std::vector<int>();
}

double Network::count_edges() const {
    double count = 0.0;
    for(int i=1; i<=n_actor; i++){
        count += adj_list[i].size();
    }
    return(count);
}

void Network::add_edge(int from, int to) {
    if(directed){
        if (!adj_mat[get_mat_idx(from, to)]) {
            auto& o_list = adj_list[from];
            o_list.insert(std::lower_bound(o_list.begin(), o_list.end(), to), to);
            auto& i_list = adj_list_in[to];
            i_list.insert(std::lower_bound(i_list.begin(), i_list.end(), from), from);
            adj_mat[get_mat_idx(from, to)] = 1;
            number_edges ++;
        }
    } else{
        if (!adj_mat[get_mat_idx(from, to)]) {
            auto& list_f = adj_list[from];
            list_f.insert(std::lower_bound(list_f.begin(), list_f.end(), to), to);
            auto& list_t = adj_list[to];
            list_t.insert(std::lower_bound(list_t.begin(), list_t.end(), from), from);
            adj_mat[get_mat_idx(from, to)] = 1;
            adj_mat[get_mat_idx(to, from)] = 1;
            number_edges ++;
        }
    }
}

void Network::delete_edge(int from, int to) {
    if(directed){
        if (adj_mat[get_mat_idx(from, to)]) {
            adj_list[from].erase(std::remove(adj_list[from].begin(), adj_list[from].end(), to), adj_list[from].end());
            adj_list_in[to].erase(std::remove(adj_list_in[to].begin(), adj_list_in[to].end(), from), adj_list_in[to].end());
            adj_mat[get_mat_idx(from, to)] = 0;
            number_edges --;
        }
    } else{
        if (adj_mat[get_mat_idx(from, to)]) {
            adj_list[from].erase(std::remove(adj_list[from].begin(), adj_list[from].end(), to), adj_list[from].end());
            adj_list[to].erase(std::remove(adj_list[to].begin(), adj_list[to].end(), from), adj_list[to].end());
            adj_mat[get_mat_idx(from, to)] = 0;
            adj_mat[get_mat_idx(to, from)] = 0;
            number_edges --;
        }
    }
}

void Network::add_edges_from_mat(arma::mat mat) {
    mat_to_map_vec(mat, n_actor, directed, adj_list, adj_list_in, adj_mat);
    number_edges = count_edges();
}

// XZ_class implementations
XZ_class::XZ_class(int n_actor_, bool directed_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),              
    x_attribute(n_actor_, type_, scale_)         
{
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
    overlap_mat = arma::zeros<arma::mat>(0, 2);
    
    for (int i = 1; i <= n_actor; ++i) { 
        all_actors.push_back(i);
    } 
    initialize_overlap_counts();
}

XZ_class::XZ_class(int n_actor_, bool directed_, arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),             
    x_attribute(n_actor_, type_, scale_)        
{
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);

    mat_to_map_vec(neighborhood_, n_actor, directed_, neighborhood, neighborhood, neighborhood_bool_mat);
    mat_to_map_vec(overlap_, n_actor, directed_, overlap, overlap, overlap_bool_mat);
    overlap_mat = overlap_;
    
    for (int i = 1; i <= n_actor; i++){
        all_actors.push_back(i);
    } 
    initialize_overlap_counts();
}

XZ_class::XZ_class(int n_actor_, bool directed_, std::vector<std::vector<int>> neighborhood_,
                   std::vector<std::vector<int>> overlap_,
                   arma::mat overlap_mat_, std::string type_, double scale_):
    n_actor(n_actor_),                             
    z_network(n_actor_, directed_),                
    x_attribute(n_actor_, type_, scale_),         
    overlap_mat(overlap_mat_)                     
{
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);

    for (int i = 1; i <= n_actor; i++){ 
        for(int neighbor : neighborhood_[i]) {
            neighborhood[i].push_back(neighbor);
            neighborhood_bool_mat[get_mat_idx(i, neighbor)] = 1;
        }
        for(int over : overlap_[i]) {
            overlap[i].push_back(over);
            overlap_bool_mat[get_mat_idx(i, over)] = 1;
        }
        std::sort(neighborhood[i].begin(), neighborhood[i].end());
        std::sort(overlap[i].begin(), overlap[i].end());

        adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
        if(z_network.directed){
            adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
        }
        all_actors.push_back(i);
    }
    initialize_overlap_counts();
}

XZ_class::XZ_class(int n_actor_, bool directed_, arma::mat z_network_, arma::vec x_attribute_,
                   arma::mat neighborhood_, arma::mat overlap_, std::string type_, double scale_):
    n_actor(n_actor_),                               
    z_network(n_actor_, directed_, z_network_),      
    x_attribute(n_actor_, x_attribute_, type_, scale_)
{
    overlap.resize(n_actor + 1);
    neighborhood.resize(n_actor + 1);
    adj_list_nb.resize(n_actor + 1);
    adj_list_in_nb.resize(n_actor + 1);
    overlap_bool_mat.assign(n_actor * n_actor, 0);
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);

    mat_to_map_vec(neighborhood_, n_actor, directed_, neighborhood, neighborhood, neighborhood_bool_mat);
    mat_to_map_vec(overlap_, n_actor, directed_, overlap, overlap, overlap_bool_mat);
    overlap_mat = overlap_;
    for (int i = 1; i <= n_actor; i++){
        adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
        if(z_network.directed){
            adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
        }
        all_actors.push_back(i);
    }
    initialize_overlap_counts();
}

void XZ_class::set_network_from_mat(int n_actor_, bool directed_, arma::mat mat){
    z_network.set_network_from_mat(n_actor_, directed_, mat); 
    for (int i = 1; i <= n_actor; i++){
        adj_list_nb[i] = get_intersection_vec(z_network.adj_list[i], overlap[i]);
        if(z_network.directed){
            adj_list_in_nb[i] = get_intersection_vec(z_network.adj_list_in[i], overlap[i]);
        } 
    }
    initialize_overlap_counts();
}

void XZ_class::initialize_overlap_counts() {
    if (overlap_mat.is_empty() || overlap_mat.n_cols < 2) {
        N_total_overlap = 0;
        N_1_overlap = 0;
        return;
    }
    int K = z_network.directed ? 1 : 2;
    N_total_overlap = (int)overlap_mat.n_rows / K;
    N_1_overlap = 0;
    
    for (int idx = 0; idx < (int)overlap_mat.n_rows; ++idx) {
        if (z_network.get_val((int)overlap_mat(idx, 0), (int)overlap_mat(idx, 1))) {
            N_1_overlap++;
        }
    }
    if (!z_network.directed) N_1_overlap /= 2;
}

void XZ_class::add_edge(int from, int to) {
    if(z_network.directed){
        if(!z_network.get_val(from, to)){
            z_network.add_edge(from, to);
            if(overlap_bool_mat[get_mat_idx(from, to)]){
                N_1_overlap++;
                auto& l_nb = adj_list_nb[from];
                l_nb.insert(std::lower_bound(l_nb.begin(), l_nb.end(), to), to);
                auto& li_nb = adj_list_in_nb[to];
                li_nb.insert(std::lower_bound(li_nb.begin(), li_nb.end(), from), from);
            }
        }
    } else{
        if(!z_network.get_val(from, to)){
            z_network.add_edge(from, to);
            if(overlap_bool_mat[get_mat_idx(from, to)]){
                N_1_overlap++;
                auto& l_f = adj_list_nb[from];
                l_f.insert(std::lower_bound(l_f.begin(), l_f.end(), to), to);
                auto& l_t = adj_list_nb[to];
                l_t.insert(std::lower_bound(l_t.begin(), l_t.end(), from), from);
            }
        }
    } 
}

void XZ_class::delete_edge(int from, int to) {
    if(z_network.directed){
        if(z_network.get_val(from, to)){
            z_network.delete_edge(from, to);
            if(overlap_bool_mat[get_mat_idx(from, to)]){
                N_1_overlap--;
                adj_list_nb[from].erase(std::remove(adj_list_nb[from].begin(), adj_list_nb[from].end(), to), adj_list_nb[from].end());
                adj_list_in_nb[to].erase(std::remove(adj_list_in_nb[to].begin(), adj_list_in_nb[to].end(), from), adj_list_in_nb[to].end());
            }
        }
    } else{ 
        if(z_network.get_val(from, to)){
            z_network.delete_edge(from, to);
            if(overlap_bool_mat[get_mat_idx(from, to)]){
                N_1_overlap--;
                adj_list_nb[from].erase(std::remove(adj_list_nb[from].begin(), adj_list_nb[from].end(), to), adj_list_nb[from].end());
                adj_list_nb[to].erase(std::remove(adj_list_nb[to].begin(), adj_list_nb[to].end(), from), adj_list_nb[to].end());
            }
        }
    } 
}

double XZ_class::count_edges() const {
    return z_network.count_edges();
}

double XZ_class::count_nb_edges() const {
    double count = 0.0;
    for(int i=1; i<=n_actor; i++){
        count += adj_list_nb[i].size();
    }
    return(count);
}

std::vector<int> XZ_class::get_common_partners_nb(unsigned int from,unsigned int to, std::string type) const {
    if(type == "OTP"){
        return(get_intersection_vec(adj_list_nb[from], adj_list_in_nb[to])); 
    } else  if(type == "ISP"){ 
        return(get_intersection_vec(adj_list_in_nb[from], adj_list_in_nb[to])); 
    }else  if(type == "OSP"){ 
        return(get_intersection_vec(adj_list_nb[from], adj_list_nb[to])); 
    } 
    else  if(type == "ITP"){
        return(get_intersection_vec(adj_list_in_nb[from], adj_list_nb[to])); 
    } 
    return std::vector<int>();
}

size_t XZ_class::count_common_partners_nb(unsigned int from, unsigned int to, std::string type) const {
    const std::vector<int>* l1_ptr = nullptr;
    if (type == "OTP") {
        l1_ptr = &adj_list_nb[from];
    } else if (type == "ISP") {
        l1_ptr = &adj_list_in_nb[from];
    } else if (type == "OSP") {
        l1_ptr = &adj_list_nb[from];
    } else if (type == "ITP") {
        l1_ptr = &adj_list_in_nb[from];
    } else {
        return 0;
    }

    const std::vector<int>& l1 = *l1_ptr;
    size_t count = 0;
    if (type == "OTP") {
        for (int k : l1) if (z_network.adj_mat[get_mat_idx(k, to)] && overlap_bool_mat[get_mat_idx(k, to)]) count++;
    } else if (type == "ISP") {
        for (int k : l1) if (z_network.adj_mat[get_mat_idx(k, to)] && overlap_bool_mat[get_mat_idx(k, to)]) count++;
    } else if (type == "OSP") {
        for (int k : l1) if (z_network.adj_mat[get_mat_idx(to, k)] && overlap_bool_mat[get_mat_idx(to, k)]) count++;
    } else if (type == "ITP") {
        for (int k : l1) if (z_network.adj_mat[get_mat_idx(to, k)] && overlap_bool_mat[get_mat_idx(to, k)]) count++;
    }
    return count;
}

bool XZ_class::check_if_full_neighborhood() const {
    for(int i = 1; i <= n_actor; ++i) {
        if(neighborhood[i].size() != static_cast<size_t>(n_actor)) {
            return false; 
        }
    }
    return true;
}

void XZ_class::print() {
    Rcout << "Network: Implemented as dense flat structures" << std::endl;
}

void XZ_class::copy_from(XZ_class obj) {
    z_network = obj.z_network;
    x_attribute = obj.x_attribute;
    neighborhood = obj.neighborhood;
    overlap = obj.overlap;
    adj_list_nb = obj.adj_list_nb;
    adj_list_in_nb = obj.adj_list_in_nb;
    overlap_bool_mat = obj.overlap_bool_mat;
    neighborhood_bool_mat = obj.neighborhood_bool_mat;
    n_actor = obj.n_actor;
    overlap_mat = obj.overlap_mat;
    all_actors = obj.all_actors;
    N_total_overlap = obj.N_total_overlap;
    N_1_overlap = obj.N_1_overlap;
}

void XZ_class::set_neighborhood_from_mat(arma::mat mat) {
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
    for (int i = 1; i <= n_actor; i++){
        neighborhood[i].clear();
        for(int j = 1; j <= n_actor; j++) {
            if(mat(i-1, j-1) == 1) {
                neighborhood[i].push_back(j);
                neighborhood_bool_mat[get_mat_idx(i, j)] = 1;
            }
        }
    }
}

void XZ_class::neighborhood_initialize() {
    for (int i = 1; i <= n_actor; i++){
        neighborhood[i].clear();
    }
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
}

void XZ_class::assign_neighborhood(const std::unordered_map< int, std::unordered_set<int>>& new_neighborhood) {
    neighborhood_bool_mat.assign(n_actor * n_actor, 0);
    for (int i = 1; i <= n_actor; i++){
        neighborhood[i].clear();
        for(int neighbor : new_neighborhood.at(i)) {
            neighborhood[i].push_back(neighbor);
            neighborhood_bool_mat[get_mat_idx(i, neighbor)] = 1;
        }
        std::sort(neighborhood[i].begin(), neighborhood[i].end());
    }
}

void XZ_class::change_neighborhood(int actor, std::unordered_set<int> new_neighborhood) {
    for(int old_n : neighborhood[actor]) {
        neighborhood_bool_mat[get_mat_idx(actor, old_n)] = 0;
    }
    neighborhood[actor].clear();
    for(int new_n : new_neighborhood) {
        neighborhood[actor].push_back(new_n);
        neighborhood_bool_mat[get_mat_idx(actor, new_n)] = 1;
    }
    std::sort(neighborhood[actor].begin(), neighborhood[actor].end());
}

// XYZ_class implementations
void XYZ_class::print() {
    Rcout << "XYZ Class: Implemented as dense flat structures" << std::endl;
}

void XYZ_class::set_info_arma(arma::vec x_attribute_, arma::vec y_attribute_, arma::mat z_network_) {
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

void XYZ_class::copy_from(const XYZ_class& obj) {
    XZ_class::copy_from(obj);
    y_attribute = obj.y_attribute;
}
