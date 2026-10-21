#include <Rcpp.h>
#include <dlfcn.h>  // For dynamic loading of shared libraries
#ifdef _OPENMP
#include <omp.h>    // For OpenMP
#endif

// Define a function pointer type for the tims_ccs_to_oneoverk0_for_mz function
typedef double (*tims_ccs_to_oneoverk0_for_mz_t)(double ccs, int charge, double mz);

// Load the Bruker library once per R session. Repeated dlopen/dlclose (e.g.
// once per batch inside getqlist) is wasteful, and dlclose while converted
// values may still be referenced is risky. The handle is intentionally kept
// open until the R process exits.
static tims_ccs_to_oneoverk0_for_mz_t load_tims_ccs_to_oneoverk0(const std::string& lib_path) {
  static void* handle = nullptr;
  static tims_ccs_to_oneoverk0_for_mz_t fun = nullptr;
  static std::string loaded_path;

  if (fun == nullptr || loaded_path != lib_path) {
    if (handle != nullptr) {
      dlclose(handle);
      handle = nullptr;
      fun = nullptr;
    }
    handle = dlopen(lib_path.c_str(), RTLD_LAZY);
    if (!handle) {
      Rcpp::stop("Cannot open library: " + std::string(dlerror()));
    }
    dlerror(); // Clear any existing errors
    fun = (tims_ccs_to_oneoverk0_for_mz_t) dlsym(handle, "tims_ccs_to_oneoverk0_for_mz");
    const char* dlsym_error = dlerror();
    if (dlsym_error) {
      fun = nullptr;
      dlclose(handle);
      handle = nullptr;
      Rcpp::stop("Cannot load symbol 'tims_ccs_to_oneoverk0_for_mz': " + std::string(dlsym_error));
    }
    loaded_path = lib_path;
  }
  return fun;
}

// [[Rcpp::export]]
Rcpp::NumericVector ccs_to_one_over_k0_parallel(Rcpp::NumericVector ccs, 
                                               Rcpp::IntegerVector charge, 
                                               Rcpp::NumericVector mz, 
                                               std::string lib_path) {
  int n = ccs.size();
  
  // Check that inputs are the same length
  if (charge.size() != n || mz.size() != n) {
    Rcpp::stop("Input vectors must have the same length");
  }
  
  // Resolve the converter once (single threaded) before the parallel loop
  tims_ccs_to_oneoverk0_for_mz_t tims_ccs_to_oneoverk0_for_mz =
      load_tims_ccs_to_oneoverk0(lib_path);
  
  // Create a vector to store the results
  Rcpp::NumericVector result(n);
  
  // Parallel for loop using OpenMP
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < n; ++i) {
    result[i] = tims_ccs_to_oneoverk0_for_mz(ccs[i], charge[i], mz[i]);
  }
  
  return result;
}
