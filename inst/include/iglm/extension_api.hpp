#pragma once
#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <mutex>
#include <stdexcept>

class XYZ_class;

namespace iglm {
// --- Function type ---
using ExtFn = double(*)(const ::XYZ_class&,const int&,const int&,const arma::mat&,
                     const double&,const std::string&,const bool&);


// --- Meta info attached to each function ---
struct FUN {
  ExtFn fn;
  std::string short_name;
  double value;
};

// --- Registry singleton ---
class Registry {
public:
  static Registry& instance() {
    static Registry inst;
    return inst;
  }
  
  bool add(const std::string& name,
           ExtFn fn,
           const std::string& short_name,
           double value)
  {
    std::lock_guard<std::mutex> lock(mu_);
    return map_.emplace(name, FUN{fn, short_name, value}).second;
  }
  
  bool has(const std::string& name) const {
    std::lock_guard<std::mutex> lock(mu_);
    return map_.count(name);
  }
  
  ExtFn get(const std::string& name) const {
    std::lock_guard<std::mutex> lock(mu_);
    auto it = map_.find(name);
    if (it == map_.end())
      throw std::out_of_range("No extension named '" + name + "'");
    return it->second.fn;
  }
  
  FUN info(const std::string& name) const {
    std::lock_guard<std::mutex> lock(mu_);
    auto it = map_.find(name);
    if (it == map_.end())
      throw std::out_of_range("No extension named '" + name + "'");
    return it->second;
  }
  
  std::vector<std::string> names() const {
    std::lock_guard<std::mutex> lock(mu_);
    std::vector<std::string> out;
    out.reserve(map_.size());
    for (auto& kv : map_) out.push_back(kv.first);
    return out;
  }
  
  std::vector<FUN> all_meta() const {
    std::lock_guard<std::mutex> lock(mu_);
    std::vector<FUN> out;
    out.reserve(map_.size());
    for (auto& kv : map_) out.push_back(kv.second);
    return out;
  }
  
private:
  Registry() = default;
  mutable std::mutex mu_;
  std::unordered_map<std::string, FUN> map_;
};

// --- Registrar helper (runs at static initialization) ---
struct Registrar {
  Registrar(const std::string& name,
            ExtFn fn,
            const std::string& short_name,
            double value)
  {
    if (!Registry::instance().add(name, fn, short_name, value)) {
      Rcpp::Rcerr << "Duplicate extension name '" << name << "' ignored.\n";
    }
  }
};

// ---- helper macros to build unique identifiers ----
#define iglm_JOIN_IMPL(a,b) a##b
#define iglm_JOIN(a,b)      iglm_JOIN_IMPL(a,b)

#if defined(__COUNTER__)  // works on GCC/Clang/MSVC
#define iglm_UNIQ(prefix) iglm_JOIN(prefix, __COUNTER__)
#else                      // portable fallback
#define iglm_UNIQ(prefix) iglm_JOIN(prefix, __LINE__)
#endif

// ---- your Registrar stays the same ----
// struct Registrar { Registrar(const std::string&, ExtFn, const std::string&, double); ... };

// ---- fixed registration macro ----
#define EFFECT_REGISTER(NAME, FN, SHORT, VAL) \
static ::iglm::Registrar iglm_UNIQ(_iglm_registrar_){ (NAME), (FN), (SHORT), (VAL) }

} // namespace iglm
