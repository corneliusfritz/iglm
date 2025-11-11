#pragma once

#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <mutex>
#include <stdexcept>
#include <vector>

/*
 * This logic handles symbol visibility.
 * - When your package (iglm) is compiling, it must define
 * IGLM_COMPILING_IGLM. This tells the compiler to EXPORT symbols.
 * - When another package (OtherPackage) includes this header, this
 * define will be absent, and the compiler will IMPORT symbols.
 *
 * You set IGLM_COMPILING_IGLM in your src/Makevars file.
 */
#if defined(_WIN32)
#ifdef IGLM_COMPILING_IGLM
#define IGLM_API __declspec(dllexport)
#else
#define IGLM_API __declspec(dllimport)
#endif
#else
// On Linux/macOS, R's symbol resolution handles this.
#define IGLM_API
#endif
class XYZ_class; 

#include "iglm/xyz_class.h"

namespace iglm {

// --- Function type (Unchanged) ---
using ExtFn = double(*)(const ::XYZ_class&,const int&,const int&,const arma::mat&,
                     const double&,const std::string&,const bool&);


struct FUN {
  ExtFn fn;
  std::string short_name;
  double value;
};


class IGLM_API Registry {
public:
  
  static Registry& instance();
  bool add(const std::string& name,
           ExtFn fn,
           const std::string& short_name,
           double value);
  
  bool has(const std::string& name) const;
  
  ExtFn get(const std::string& name) const;
  
  FUN info(const std::string& name) const;
  
  std::vector<std::string> names() const;
  
  std::vector<FUN> all_meta() const;
  
private:
  Registry() = default;
  Registry(const Registry&) = delete;
  Registry& operator=(const Registry&) = delete;
  
  mutable std::mutex mu_;
  std::unordered_map<std::string, FUN> map_;
}; 

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

#define iglm_JOIN_IMPL(a,b) a##b
#define iglm_JOIN(a,b)      iglm_JOIN_IMPL(a,b)

#if defined(__COUNTER__)  // works on GCC/Clang/MSVC
#define iglm_UNIQ(prefix) iglm_JOIN(prefix, __COUNTER__)
#else                     // portable fallback
#define iglm_UNIQ(prefix) iglm_JOIN(prefix, __LINE__)
#endif

#define EFFECT_REGISTER(NAME, FN, SHORT, VAL) \
static ::iglm::Registrar iglm_UNIQ(_iglm_registrar_){ (NAME), (FN), (SHORT), (VAL) }

} // namespace iglm