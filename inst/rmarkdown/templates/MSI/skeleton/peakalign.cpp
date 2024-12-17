#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <string>
#include <unordered_map>
using namespace Rcpp;

struct ThreadLocalResult {
  std::vector<int> frame;
  std::vector<std::string> mz_ccs;
  std::vector<double> intensity;
};
// [[Rcpp::export]]
DataFrame findpeakalign(NumericVector mz, NumericVector ccs, NumericVector intensity,
                        IntegerVector frame, NumericVector ref_mz, NumericVector ref_ccs,
                        double ppm, double ccs_shift) {
  
  std::unordered_map<int, std::vector<int>> frame_to_index;
  for (int i = 0; i < frame.size(); ++i) {
    frame_to_index[frame[i]].push_back(i);
  }
  
  std::vector<int> unique_frames = as<std::vector<int>>(unique(frame));
  std::sort(unique_frames.begin(), unique_frames.end());
  
  std::vector<ThreadLocalResult> thread_results(

#ifdef _OPENMP
      omp_get_max_threads()
#else
    1
#endif
  );
  
#pragma omp parallel for schedule(guided)
  for (int i = 0; i < ref_mz.size(); i++) {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    ThreadLocalResult& local_result = thread_results[thread_id];
    
    double mz_tol = ppm * ref_mz[i] / 1e6;
    double ccs_tol = ccs_shift * ref_ccs[i];
    std::string mz_ccs_str = std::to_string(ref_mz[i]) + "_" + std::to_string(ref_ccs[i]);
    
    for (int current_frame : unique_frames) {
      const std::vector<int>& indices = frame_to_index[current_frame];
      
      double total_intensity = 0.0;
      
      for (int j : indices) {
        if (std::abs(mz[j] - ref_mz[i]) <= mz_tol &&
            std::abs(ccs[j] - ref_ccs[i]) <= ccs_tol) {
          total_intensity += intensity[j];
        }
      }
      
      if (total_intensity > 0.0) {
        local_result.frame.push_back(current_frame);
        local_result.mz_ccs.push_back(mz_ccs_str);
        local_result.intensity.push_back(total_intensity);
      }
    }
  }
  
  size_t total_size = 0;
  for (const auto& result : thread_results) {
    total_size += result.frame.size();
  }
  
  std::vector<int> out_frame;
  std::vector<std::string> out_mz_ccs;
  std::vector<double> out_intensity;
  
  out_frame.reserve(total_size);
  out_mz_ccs.reserve(total_size);
  out_intensity.reserve(total_size);
  
  for (const auto& result : thread_results) {
    out_frame.insert(out_frame.end(), result.frame.begin(), result.frame.end());
    out_mz_ccs.insert(out_mz_ccs.end(), result.mz_ccs.begin(), result.mz_ccs.end());
    out_intensity.insert(out_intensity.end(), result.intensity.begin(), result.intensity.end());
  }
  
  return DataFrame::create(
    Named("frame") = out_frame,
    Named("mz_ccs") = out_mz_ccs,
    Named("intensity") = out_intensity
  );
}
