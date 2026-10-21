#include <Rcpp.h>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>

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

  // Group point indices by frame
  std::unordered_map<int, std::vector<int>> frame_to_index;
  for (int i = 0; i < frame.size(); ++i) {
    frame_to_index[frame[i]].push_back(i);
  }

  std::vector<int> unique_frames;
  unique_frames.reserve(frame_to_index.size());
  for (const auto& kv : frame_to_index) {
    unique_frames.push_back(kv.first);
  }
  std::sort(unique_frames.begin(), unique_frames.end());

  // Sort the indices of each frame by m/z so that the matching below can use
  // binary search: each reference peak only visits the points that fall inside
  // its m/z window instead of scanning every point of the frame
  // (O(N_ref x N_raw) -> O(N_raw log N_raw + N_ref x log N_raw + hits)).
  // Plain pointers are stored so the parallel loop only performs read accesses
  // on the shared map data.
  std::vector<const std::vector<int>*> frame_indices;
  frame_indices.reserve(unique_frames.size());
  for (int f : unique_frames) {
    std::vector<int>& idxv = frame_to_index[f];
    std::sort(idxv.begin(), idxv.end(),
              [&mz](int a, int b) { return mz[a] < mz[b]; });
    frame_indices.push_back(&idxv);
  }

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
    double mz_lo = ref_mz[i] - mz_tol;
    double mz_hi = ref_mz[i] + mz_tol;
    std::string mz_ccs_str = std::to_string(ref_mz[i]) + "_" + std::to_string(ref_ccs[i]);

    for (size_t k = 0; k < frame_indices.size(); ++k) {
      const std::vector<int>& indices = *frame_indices[k];

      // first point index with mz >= mz_lo
      auto begin = std::lower_bound(indices.begin(), indices.end(), mz_lo,
                                    [&mz](int a, double v) { return mz[a] < v; });

      double total_intensity = 0.0;
      for (auto it = begin; it != indices.end() && mz[*it] <= mz_hi; ++it) {
        if (std::abs(ccs[*it] - ref_ccs[i]) <= ccs_tol) {
          total_intensity += intensity[*it];
        }
      }

      if (total_intensity > 0.0) {
        local_result.frame.push_back(unique_frames[k]);
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
