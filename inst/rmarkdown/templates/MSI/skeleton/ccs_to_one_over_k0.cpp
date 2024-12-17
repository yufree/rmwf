#include <Rcpp.h>
#include <dlfcn.h>  // For dynamic loading of shared libraries
#include <omp.h>    // For OpenMP

// Define a function pointer type for the tims_ccs_to_oneoverk0_for_mz function
typedef double (*tims_ccs_to_oneoverk0_for_mz_t)(double ccs, int charge, double mz);

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
  
  // Load the shared library
  void* handle = dlopen(lib_path.c_str(), RTLD_LAZY);
  if (!handle) {
    Rcpp::stop("Cannot open library: " + std::string(dlerror()));
  }
  
  // Load the function from the shared library
  dlerror(); // Clear any existing errors
  tims_ccs_to_oneoverk0_for_mz_t tims_ccs_to_oneoverk0_for_mz = 
    (tims_ccs_to_oneoverk0_for_mz_t) dlsym(handle, "tims_ccs_to_oneoverk0_for_mz");
  const char* dlsym_error = dlerror();
  if (dlsym_error) {
    dlclose(handle);
    Rcpp::stop("Cannot load symbol 'tims_ccs_to_oneoverk0_for_mz': " + std::string(dlsym_error));
  }
  
  // Create a vector to store the results
  Rcpp::NumericVector result(n);
  
  // Parallel for loop using OpenMP
#pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    result[i] = tims_ccs_to_oneoverk0_for_mz(ccs[i], charge[i], mz[i]);
  }
  
  // Close the library
  dlclose(handle);
  
  return result;
}