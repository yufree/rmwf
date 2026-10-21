library(opentimsr)
library(data.table)
library(Rcpp)
library(enviGCMS)
library(irlba)
library(ClusterR)
library(reticulate)
library(umap)
library(dbscan)
library(pmd)
library(igraph)
library(ggraph)
library(tidygraph)
library(flexdashboard)
library(DT)
library(shinyWidgets)
library(shinydashboard)
library(visNetwork)

sourceCpp("one_over_k0_to_ccs.cpp")
sourceCpp('peakalign.cpp')
sourceCpp("peak_finder.cpp")
#'@title Find 2D Peaks
#'@description This function finds 2D peaks in mass spectrometry data.
#'@param mz A numeric vector of m/z values.
#'@param ccs A numeric vector of CCS values.
#'@param intensity A numeric vector of intensity values.
#'@param mz_ppm The m/z tolerance in ppm. Default is 20.
#'@param ccs_tolerance The CCS tolerance. Default is 0.03.
#'@param snr The signal-to-noise ratio cutoff. Default is 3.
#'@param mz_bins The number of m/z bins. Default is 1000.
#'@param ccs_bins The number of CCS bins. Default is 50.
#'@return A data.table with the 2D peaks.
find_2d_peaks <- function(mz, ccs, intensity, 
                          mz_ppm = 20, 
                          ccs_tolerance = 0.03,
                          snr = 3.0,
                          mz_bins = 1000,
                          ccs_bins = 50) {
  
  if (!is.numeric(mz) || !is.numeric(ccs) || !is.numeric(intensity)) {
    stop("All inputs must be numeric vectors")
  }
  if (length(mz) != length(ccs) || length(mz) != length(intensity)) {
    stop("All input vectors must have the same length")
  }
  
  find_2d_peaks_spatial_parallel_openmp(mz, ccs, intensity, 
                                        mz_ppm, ccs_tolerance, snr,
                                        mz_bins, ccs_bins)
}

#'@title Get Summary
#'@description This function provides a summary of the normalized peak list.
#'@param peak_file The path to the normalized peak list file.
#'@return A summary of the normalized peak list.
getsummary <- function(peak_file) {
  # Read the normalized peak list
  dt <- data.table::fread(peak_file, header = TRUE)
  
  # Extract mz and ccs from column names
  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  ccs <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))
  
  # Extract location coordinates
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  # Calculate ranges
  mz_range <- range(mz, na.rm = TRUE)
  ccs_range <- range(ccs, na.rm = TRUE)
  x_range <- range(x, na.rm = TRUE)
  y_range <- range(y, na.rm = TRUE)
  
  # Print the summary
  cat("Summary of normalized peak list:\n")
  cat("  m/z range: ", mz_range[1], " - ", mz_range[2], "\n")
  cat("  CCS range: ", ccs_range[1], " - ", ccs_range[2], "\n")
  cat("  Location X range: ", x_range[1], " - ", x_range[2], "\n")
  cat("  Location Y range: ", y_range[1], " - ", y_range[2], "\n")
}

#'@title Safe Query for OpenTIMS with Binary-Split Error Recovery
#'@description Queries frames from an OpenTIMS handle with automatic recursive binary-split isolation and skipping of corrupt/unreadable frames.
#'@param handle OpenTIMS data handle.
#'@param frames Integer vector of frame IDs.
#'@param columns Character vector of column names.
#'@param max_depth Maximum recursion depth. Default is 8.
#'@return Named list of column vectors, or NULL if all frames are unreadable.
safe_query <- function(handle, frames, columns, max_depth = 8) {
  res <- tryCatch(
    query(handle, frames = frames, columns = columns),
    error = function(e) e
  )
  if (!inherits(res, "error")) return(res)
  if (length(frames) <= 1 || max_depth <= 0) {
    warning(sprintf("Skipped unreadable/corrupted frame(s): %s", paste(frames, collapse = ",")))
    return(NULL)
  }
  mid <- floor(length(frames) / 2)
  f1 <- frames[1:mid]
  f2 <- frames[(mid + 1):length(frames)]
  d1 <- safe_query(handle, f1, columns, max_depth - 1)
  d2 <- safe_query(handle, f2, columns, max_depth - 1)
  if (is.null(d1) && is.null(d2)) return(NULL)
  if (is.null(d1)) return(d2)
  if (is.null(d2)) return(d1)
  lapply(names(d1), function(nm) c(d1[[nm]], d2[[nm]]))
}

#'@title Get Reference Peaks
#'@description This function processes raw mass spectrometry data to generate a reference peak list.
#'@param path The path to the raw data (e.g., .d folder).
#'@param accept_Bruker_EULA_and_on_Windows_or_Linux Set to TRUE to accept the Bruker EULA.
#'@param libpath The path to the Bruker library files.
#'@param batch_size The number of frames to process in each batch. Default is 25000 (optimized for 50GB RAM environments).
#'@param outref The output file for the reference peaks.
#'@param outcoord The output file for the coordinates.
#'@param xrange The x-range of the region of interest.
#'@param yrange The y-range of the region of interest.
#'@return A reference peak list.
getrefpeak <- function(path,accept_Bruker_EULA_and_on_Windows_or_Linux,libpath = './', batch_size = 25000, outref, outcoord, xrange,yrange){
  # columns to collect
  all_columns <- c('frame', 'intensity', 'mz', 'inv_ion_mobility')
  # download dll/so file if the EULA is accepted; otherwise fall back to an
  # already downloaded library in libpath
  if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
    folder_to_stode_priopriatary_code = libpath
    if (!dir.exists(folder_to_stode_priopriatary_code)) {
      dir.create(folder_to_stode_priopriatary_code, recursive = TRUE, showWarnings = FALSE)
    }
    so_file <- if (.Platform$OS.type == 'windows') 'timsdata.dll' else 'libtimsdata.so'
    target_so <- file.path(folder_to_stode_priopriatary_code, so_file)
    if (!file.exists(target_so)) {
      try({
        py_lib <- system(sprintf("python3 -c 'import opentims_bruker_bridge, pathlib; print(pathlib.Path(opentims_bruker_bridge.__file__).parent / \"%s\")'", so_file), intern = TRUE)
        if (file.exists(py_lib)) file.copy(py_lib, target_so)
      }, silent = TRUE)
    }
    if (!file.exists(target_so)) {
      path_to_bruker_dll = download_bruker_proprietary_code(folder_to_stode_priopriatary_code)
    } else {
      path_to_bruker_dll = target_so
    }
    setup_bruker_so(path_to_bruker_dll)
    lib_file <- path_to_bruker_dll
  }else{
    # libtimsdata.dll on windows, libtimsdata.so on other platforms
    so_file <- if (.Platform$OS.type == 'windows') 'libtimsdata.dll' else 'libtimsdata.so'
    lib_file <- file.path(libpath, so_file)
  }
  # Build the connection with raw data
  D = OpenTIMS(path)
  # extract frame information or define a range to subset the data
  xy <- table2df(D, "MaldiFrameInfo")
  # subset coords idx for region of interest (ROI): frames INSIDE the range
  roi_idx <- which(xy$MaldiFrameInfo$XIndexPos>xrange[1]&xy$MaldiFrameInfo$XIndexPos<xrange[2]&xy$MaldiFrameInfo$YIndexPos>yrange[1]&xy$MaldiFrameInfo$YIndexPos<yrange[2])
  # set a batch size to avoid memory issues by process the data in batches
  total_frames <- length(roi_idx)
  frame <- xy$MaldiFrameInfo$Frame[roi_idx]

  # Process data in batches
  total_batches <- ceiling(total_frames / batch_size)

  dtx_list <- list()
  for (i in 1:total_batches) {
    batch_idx <- ((i - 1) * batch_size + 1):min(i * batch_size, total_frames)

    data <- safe_query(D, frames = frame[batch_idx], columns = all_columns)
    if (is.null(data) || length(data$frame) == 0) next
    data <- as.data.table(data)
    # round up to make sure binning is working and fit the mass/ion mobility accuracy
    data[, mz := round(mz, 4)]
    data[, inv_ion_mobility := round(inv_ion_mobility, 3)]

    dt <- data[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]
    dtx_list[[length(dtx_list) + 1]] <- dt
    rm(data, dt)
    gc(verbose = FALSE)
  }
  # reshape the data by mz and ion mobility
  dtx <- rbindlist(dtx_list)
  dtx <- dtx[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]

  # add ccs value for super pixel peaks
  result <- one_over_k0_to_ccs_parallel(dtx$inv_ion_mobility, rep(1,length(dtx$inv_ion_mobility)), dtx$mz, lib_file)
  dtx[,ccs:=result]
  # save super pixel peaks
  fwrite(dtx,outref)
  # save the location of the pixels inside the ROI; getqlist() uses this file
  # to map frame numbers back to pixel locations
  coord <- xy$MaldiFrameInfo[roi_idx,c('Frame','XIndexPos','YIndexPos')]
  fwrite(coord,outcoord)
}

#'@title Get Quantitative Peak List
#'@description This function generates a quantitative peak list.
#'@param refpath The path to the reference peaks file.
#'@param lib_path The path to the Bruker library file (full path to libtimsdata.so/libtimsdata.dll).
#'@param path The path to the raw data.
#'@param method The normalization method ('tic' or 'rms').
#'@param zero_proportion_cutoff The cutoff for removing peaks with a high proportion of zeros.
#'@param coordpath The path to the coordinates file.
#'@param normpath The output file for the normalized peak list.
#'@param xrange The x-range of the region of interest.
#'@param yrange The y-range of the region of interest.
#'@param mz_ppm The m/z tolerance in ppm used to match peaks to the reference list. Default is 20.
#'@param ccs_tolerance The relative CCS tolerance used to match peaks to the reference list. Default is 0.03.
#'@param batch_size The number of frames to process in each batch. Default is 25000 (optimized for 50GB RAM environments).
#'@return A quantitative peak list.
getqlist <- function(refpath,lib_path,path, method, zero_proportion_cutoff, coordpath, normpath, xrange, yrange, mz_ppm = 20, ccs_tolerance = 0.03, batch_size = 25000){
  ref <- fread(refpath)
  setup_bruker_so(lib_path)
  all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
  D = OpenTIMS(path)
  xy <- table2df(D, "MaldiFrameInfo")
  idx <- which(xy$MaldiFrameInfo$XIndexPos>xrange[1]&xy$MaldiFrameInfo$XIndexPos<xrange[2]&xy$MaldiFrameInfo$YIndexPos>yrange[1]&xy$MaldiFrameInfo$YIndexPos<yrange[2])
  # set a batch size to avoid memory issues by process the data in batches
  total_frames <- dim(xy$MaldiFrameInfo[idx,])[1]
  frame <- xy$MaldiFrameInfo$Frame[idx]
  # this setting will make sure the data.frame's row number is within the limits.
  # the upper cap avoids loading too many frames in a single query when ref is small
  batch <- max(1, min(batch_size, floor(2e+09/nrow(ref))))
  batchs <- ceiling(total_frames/batch)
  result <- list()
  for(i in 1:batchs){
    upper <- min(total_frames,i*batch)
    idx<- c(((i-1)*batch+1):upper)
    data <- safe_query(D, frames = frame[idx], columns = all_columns)
    if (is.null(data) || length(data$frame) == 0) next
    data <- as.data.table(data)
    data[,ccs:=one_over_k0_to_ccs_parallel(data$inv_ion_mobility, rep(1,length(data$inv_ion_mobility)), data$mz, lib_path)]
    resultt <- findpeakalign(data$mz,data$ccs,data$intensity,data$frame,ref$mz,ref$ccs, ppm = mz_ppm, ccs_shift = ccs_tolerance)
    resultt <- as.data.table(resultt)
    if (nrow(resultt) == 0) {
      rm(data, resultt)
      next
    }
    if(method == 'tic'){
      # perform TIC normalization: divide by the total matched intensity of the
      # pixel (note: this only covers intensities matched to the reference list)
      resultt[, TIC := sum(intensity, na.rm = TRUE), by = frame]
      resultt[, normalized_intensity := as.numeric(intensity / TIC), by = frame]
      resultt[, TIC := NULL]
    }else{
      # perform RMS normalization: root mean square of the matched intensities
      resultt[, rms := sqrt(mean(intensity^2, na.rm = TRUE)), by = frame]
      resultt[, normalized_intensity := as.numeric(intensity / rms), by = frame]
      resultt[, rms := NULL]
    }
    result[[length(result) + 1]] <- dcast(
      resultt, 
      frame ~ mz_ccs, 
      value.var = "normalized_intensity", 
      fun.aggregate = sum,
      fill = 0)
    rm(data, resultt)
    gc(verbose = FALSE)
  }
  result <- rbindlist(result,fill=T)
  if (ncol(result) < 2) {
    stop("No peak in the reference list was matched to the raw data; check the mz_ppm/ccs_tolerance settings.")
  }
  for (col in names(result)) {
    set(result, i = which(is.na(result[[col]])), j = col, value = 0)
  }
  # remove the peaks with zero_proportion larger than 95%
  zero_proportion <- result[, lapply(.SD, function(x) mean(x == 0)), .SDcols = 2:ncol(result)]
  columns_to_keep <- names(zero_proportion)[zero_proportion <= zero_proportion_cutoff]
  cn <- columns_to_keep[!is.na(columns_to_keep)]
  result_filtered <- result[, ..cn, with = FALSE]
  result_filtered[,frame:=result$frame]
  # save the data with pixel location information
  coord <- fread(coordpath)
  coord[,location:=paste0(coord$XIndexPos,'_',coord$YIndexPos)]
  result_filtered[,location:=coord$location[match(result_filtered$frame,coord$Frame)]]
  setcolorder(result_filtered, c("location", setdiff(names(result_filtered), "location")))
  result_filtered[, frame := NULL]
  fwrite(result_filtered,normpath)
}

#'@title Get Annotation
#'@description This function annotates peaks against a database.
#'@param database The path to the database file.
#'@param mode The ionization mode ('pos' or 'neg').
#'@param peakpath The path to the peak list file.
#'@param annofile The output file for the annotations.
#'@param ppm The m/z tolerance in ppm. Default is 20 (suitable for TOF instruments such as the timsTOF fleX).
#'@param deltaccs The absolute CCS tolerance in squared angstrom (Å²). Default is 5.
#'@return An annotated peak list.
getanno <- function(database,mode,peakpath,annofile, ppm = 20, deltaccs = 5){
  lipid <- fread(database)
  # set mode
  if(mode == 'pos'){
    lipid <- lipid[adducts%in%c('[M+H]','[M+Na]','[M+NH4]','[M+H-H2O]')]
  }else{
    lipid <- lipid[adducts%in%c('[M-H]','[M-H+FA]','[M+Cl]')]
  }
  
  # load ref peaks
  ref <- fread(peakpath)
  mz <- sapply(strsplit(colnames(ref)[-1],'\\_'),function(x) round(as.numeric(x[1]),4))
  im <- sapply(strsplit(colnames(ref)[-1],'\\_'),function(x) as.numeric(x[2]))
  # align
  align <- enviGCMS::getalign(mz,lipid$mz,im,lipid$ccs,ppm=ppm,deltart = deltaccs)
  # getalign() returns NULL when no peak is within tolerance; cbind.data.frame()
  # then silently drops the NULL db_mz/db_im columns, and the merge() below
  # fails on the missing keys. Short-circuit to a well-formed empty file.
  if (is.null(align) || is.null(align$xid) || length(align$xid) == 0) {
    message("No peak matched the database ", database, "; writing an empty annotation file.")
    anno <- data.frame(mz = numeric(), im = numeric(),
                       db_mz = numeric(), db_im = numeric())
    lipidanno <- merge(anno, lipid, by.x = c('db_mz','db_im'), by.y = c('mz','ccs'))
    fwrite(lipidanno, annofile)
    return(invisible(lipidanno))
  }
  anno <- cbind.data.frame(mz=mz[align$xid],im=im[align$xid],db_mz=align$mz2,db_im=align$rt2)
  lipidanno <- merge(anno,lipid,by.x = c('db_mz','db_im'),by.y = c('mz','ccs'))
  lipidanno <- lipidanno[!duplicated(lipidanno),]
  # save annotation result as csv file
  fwrite(lipidanno, annofile)
}

#'@title Plot Peak Statistics
#'@description This function generates and displays plots for peak statistics.
#'@param peak_file The path to the input peak file.
#'@return Plots for peak statistics.
plot_peak_stats <- function(peak_file) {
  library(data.table)
  library(ggplot2)
  
  dt <- fread(peak_file, header = TRUE)
  # count non-zero peaks per pixel; exclude the character 'location' column
  # (first column), otherwise 'location != 0' is always TRUE and inflates np by 1
  dt[, np := rowSums(.SD != 0), .SDcols = -1]
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  df <- cbind.data.frame(x = x, y = y, np = dt$np)
  
  p1 <- ggplot(df, aes(np)) +
    ggtitle('Peak number for each pixel') + xlab('Peak number') +
    geom_histogram(binwidth = 1) + theme_bw()
  
  p2 <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = np), size = 0.001) +
    scale_color_gradient(low = "white", high = "red") + theme_void()
  
  print(p1)
  print(p2)
}

#'@title Perform PCA Segmentation
#'@description This function performs PCA segmentation on the peak data.
#'@param peak_file The input peak file.
#'@param output_file The output file for the segmentation results.
#'@param n_components The number of principal components to use. Default is 20.
#'@param n_clusters The number of clusters to create. Default is 5.
#'@return A segmentation plot and a segmentation file.
perform_pca_segmentation <- function(peak_file, output_file, n_components = 20, n_clusters = 5) {
  
  dt <- fread(peak_file, header = TRUE)
  dt_values <- as.matrix(dt[, -1, with = FALSE])   # rows = pixels, cols = ions
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  # center each ion (column) across pixels, then project the PIXELS (rows)
  # onto the principal components; previous versions transposed the matrix
  # several times and ended up clustering ions instead of pixels
  mat_centered <- scale(dt_values, center = TRUE, scale = FALSE)
  svd_result <- irlba(mat_centered, nv = n_components)
  pca_scores <- mat_centered %*% svd_result$v      # rows = pixels, cols = PCs
  
  km <- KMeans_arma(as.matrix(pca_scores), clusters = n_clusters, n_iter = 10, seed_mode = "random_subset", verbose = TRUE, CENTROIDS = NULL)
  pr <- predict_KMeans(as.matrix(pca_scores), km)
  
  plot(x, y, col = pr, cex = 0.1, pch = 19)
  legend('topright', legend = unique(pr), col = unique(pr), pch = 19, cex = 1)
  
  seg <- cbind.data.frame(x = x, y = y, seg = pr)
  fwrite(seg, output_file)
}

#'@title Perform UMAP Segmentation
#'@description This function performs UMAP segmentation on the peak data.
#'@param peak_file The input peak file.
#'@param output_file The output file for the segmentation results.
#'@param n_threads The number of threads to use. Default is 50.
#'@param eps The epsilon parameter for DBSCAN. Default is 0.2.
#'@param minPts The minimum number of points for DBSCAN. Default is 20.
#'@return A segmentation plot and a segmentation file.
perform_umap_segmentation <- function(peak_file, output_file, n_threads = 50, eps = 0.2, minPts = 20) {
  
  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  mat <- t(dt_values)
  Sys.setenv(NUMBA_NUM_THREADS = n_threads)
  py_run_string(paste0("import os; print(os.environ['NUMBA_NUM_THREADS'])"))
  
  viz <- umap::umap(t(mat), method = 'umap-learn', metric = 'cosine')
  dbscan_result <- dbscan(viz$layout, eps = eps, minPts = minPts)
  
  plot(x = x, y = y, col = dbscan_result$cluster+1,xlab='',ylab = '',main='',xaxt = "n", yaxt = "n",bty = "n",cex=0.1)
  legend('bottomright',legend = unique(dbscan_result$cluster+1),col = unique(dbscan_result$cluster+1), pch=19,bty = "n")
  
  seg <- cbind.data.frame(x = x, y = y, seg = dbscan_result$cluster+1)
  fwrite(seg, output_file)
}

#'@title Perform Ion Clustering
#'@description This function clusters ions based on their spatial distribution.
#'@param peak_file The input peak file.
#'@param output_file The output file for the ion clusters.
#'@param locindex The index of the location to be clustered. Default is NULL.
#'@param hclust_cutoff The cutoff for hierarchical clustering. Default is 0.6.
#'@param min_cluster_size The minimum number of ions in a cluster. Default is 10.
#'@param folder_name The name of the folder to save the cluster images. Default is 'cluster'.
#'@return A file with the ion clusters and images for each cluster.
perform_cluster_ions <- function(peak_file, output_file, locindex=NULL, hclust_cutoff = 0.6, min_cluster_size = 10, folder_name='cluster') {
  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]
  
  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  if(!is.null(locindex)){
    subdt <- dt_values[locindex]
    x <- x[locindex]
    y <- y[locindex]
  }else{
    subdt <- dt_values
  }
  
  Matrix <- t(subdt)  # rows = ions, cols = pixels
  row_norms <- sqrt(rowSums(Matrix^2))
  # guard: ions that are all-zero within the ROI (e.g. ions enriched outside the
  # ROI) have norm 0 and would produce NaN similarities that crash hclust()
  row_norms[!is.finite(row_norms) | row_norms == 0] <- 1
  sim <- sweep(Matrix, 1, row_norms, "/")
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  
  hc <- hclust(D_sim)
  s <- cutree(hc, h = hclust_cutoff)
  
  name <- as.numeric(names(table(s)[table(s) > min_cluster_size]))
  
  Matrix <- t(subdt)
  split_matrices <- lapply(name, function(category) {
    rows <- which(s == as.numeric(category))
    subset_matrix <- Matrix[rows, , drop = FALSE]
    return(subset_matrix)
  })
  
  summed_matrices <- lapply(split_matrices, function(subset_matrix) {
    # z-score each ion (row) across pixels, then average over the ions of this
    # cluster for every pixel: the composite image of the cluster.
    # constant ions (sd = 0) yield NaN in scale(); treat them as zeros
    z <- t(scale(t(subset_matrix)))
    z[!is.finite(z)] <- 0
    colMeans(z)
  })
  
  result_matrix <- do.call(cbind, summed_matrices)
  
  clpan <- cbind.data.frame(x = x, y = y, result_matrix)
  
  dir.create(folder_name, showWarnings = FALSE)
  for (i in c(1:length(name))) {
    dfx <- result_matrix[, i]
    rng <- max(dfx) - min(dfx)
    norm <- if (rng > 0) (dfx - min(dfx)) / rng else rep(0, length(dfx))
    color_palette <- colorRamp(c("skyblue", "red"))
    color_sequence <- rgb(color_palette(norm) / 255, alpha = 1)
    xlim <- max(x) - min(x) + 1
    ylim <- max(y) - min(y) + 1
    
    png(paste0(folder_name,'/cluster', name[i], '.png'), width = xlim, height = ylim)
    # set margins BEFORE plot.new(): with the default margins plot.new() fails
    # on small canvases ("figure margins too large")
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, xlim), ylim = c(0, ylim), xaxs = "i", yaxs = "i", asp = NA)
    points(x - min(x) + 1, y - min(y) + 1, pch = 16, col = color_sequence, cex = 0.3)
    dev.off()
  }
  ioncluster <- cbind.data.frame(mz, im, class = s)
  fwrite(ioncluster, output_file)
}

#'@title Perform Reactomics Analysis
#'@description This function performs reactomics analysis.
#'@param peak_file The input peak file.
#'@param output_file The output file for the reactomics analysis.
#'@param pmd The PMD values to be analyzed. Default is c(14.016,15.995,2.016,18.011).
#'@param locindex The index of the location to be analyzed. Default is NULL.
#'@return A file with the reactomics analysis results.
perform_reactomics_analysis <- function(peak_file, output_file, pmd = c(14.016,15.995,2.016,18.011), locindex = NULL) {
  
  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]
  
  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  xy <- cbind.data.frame(x,y)
  
  if(!is.null(locindex)){
    subdt <- dt_values[locindex]
    xy <- xy[locindex,]
  }else{
    subdt <- dt_values
  }
  
  df <- getpmddf(mz, pmd = pmd, digits = 3)
  
  df$diff3 <- round(df$diff, 3)
  
  # Normalize the intensity data for the ROI (row-wise) and transpose it
  dfall <- t(apply(subdt, 1, function(x) {
    max_val <- max(x, na.rm = TRUE)
    if (!is.finite(max_val) || max_val == 0) {
      return(rep(0, length(x)))
    }
    return(x / max_val)
  }))
  dfall[is.na(dfall)] <- 0
  
  for (pmd_val in pmd) {
    
    # Filter the PMD data frame for the current, exact pmd_val
    dff_group <- df[df$diff3 == pmd_val, ]
    
    # Initialize sums as zero vectors
    pmd_h <- numeric(nrow(xy))
    pmd_l <- numeric(nrow(xy))
    
    # Proceed only if PMD pairs were found for this value
    if (nrow(dff_group) > 0) {
      
      # Sum intensities for higher m/z values ('ms1')
      dfms1 <- dfall[, mz %in% unique(dff_group$ms1), drop = FALSE]
      if (ncol(dfms1) > 0) {
        pmd_h <- rowSums(dfms1)
      }
      
      # Sum intensities for lower m/z values ('ms2')
      dfms2 <- dfall[, mz %in% unique(dff_group$ms2), drop = FALSE]
      if (ncol(dfms2) > 0) {
        pmd_l <- rowSums(dfms2)
      }
    }
    
    # --- Store results in the output data frame ---
    # Create a clean suffix for column names by replacing '.' with '_'
    col_suffix <- gsub(".", "_", as.character(pmd_val), fixed = TRUE)
    
    # Create dynamic column names
    h_col_name <- paste0("pmd", col_suffix, "h")
    l_col_name <- paste0("pmd", col_suffix, "l")
    total_col_name <- paste0("pmd", col_suffix)
    
    # Assign the calculated sums to the new columns
    xy[[h_col_name]] <- pmd_h
    xy[[l_col_name]] <- pmd_l
    xy[[total_col_name]] <- pmd_h + pmd_l
  }
  fwrite(xy, output_file)
}

#'@title Save Ion Images
#'@description This function saves ion images for the specified m/z values.
#'@param peak_file The input peak file.
#'@param mz_values A vector of m/z values.
#'@return Ion images for the specified m/z values.
save_ion_images <- function(peak_file, mz_values) {
  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))
  
  dt_filtered <- as.data.frame(dt_values[, .SD, .SDcols = mz %in% mz_values])
  
  width <- max(x) - min(x) + 1
  height <- max(y) - min(y) + 1
  
  for (i in 1:ncol(dt_filtered)) {
    dfx <- dt_filtered[[i]]
    rng <- max(dfx) - min(dfx)
    norm <- if (rng > 0) (dfx - min(dfx)) / rng else rep(0, length(dfx))
    png(paste0(colnames(dt_filtered)[i], '.png'), width = width, height = height)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim = c(0, width), ylim = c(0, height), xaxs = "i", yaxs = "i", asp = NA)
    points(x - min(x) + 1, y - min(y) + 1, pch = 16, col = grDevices::gray(1 - norm), cex = 0.3)
    dev.off()
  }
}

#'@title Save PMD Images
#'@description This function saves PMD images for the specified PMD values.
#'@param peak_file The input peak file.
#'@param pmd_values The PMD values to be visualized.
#'@param filename The name of the output file.
#'@return A PMD image.
save_pmd_images <- function(peak_file, pmd_values,filename){
  col_suffix <- gsub(".", "_", as.character(pmd_values), fixed = TRUE)
  lowname <- paste0('pmd',col_suffix,'l')
  highname <- paste0('pmd',col_suffix,'h')
  xy <- fread(peak_file)
  # +1 so the outermost pixels are not clipped, consistent with the other
  # imaging functions that use max - min + 1
  width = diff(range(xy$x)) + 1
  height = diff(range(xy$y)) + 1
  # guard against "figure margins too large" on small ROIs: the legend panel
  # uses mar = c(1,0,4,1), whose top margin exceeds a short canvas at 300 dpi
  png(filename, max(width * 2, 600), max(height * 2, 600), bg = 'transparent', res = 300)
  layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(0.8, 0.2))
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim = c(0,width), ylim = c(0,height), xaxs = "i", yaxs = "i", asp = NA)
  dfx <- xy[,..lowname][[1]]
  rng <- max(dfx) - min(dfx)
  norm <- (if (rng > 0) (dfx - min(dfx)) / rng else rep(0, length(dfx))) * 0.5
  # point_colors <- colors_with_alpha[col_index]
  points(xy$x-min(xy$x)+1, xy$y-min(xy$y)+1, pch = 16, col = rgb(1, 0, 0, norm),cex=0.3)
  dfx <- xy[,..highname][[1]]
  rng <- max(dfx) - min(dfx)
  norm <- (if (rng > 0) (dfx - min(dfx)) / rng else rep(0, length(dfx))) * 0.25
  points(xy$x-min(xy$x)+1, xy$y-min(xy$y)+1, pch = 16, col = rgb(0, 0, 1, norm),cex=0.3)
  par(mar = c(1, 0, 4, 1)) 
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs="i", yaxs="i")
  
  n_colors <- 100
  red_colors <- rgb(1, 0, 0, alpha = seq(1, 0, length.out = n_colors) * 0.5)
  blue_colors <- rgb(0, 0, 1, alpha = seq(1, 0, length.out = n_colors) * 0.25)
  
  rasterImage(as.raster(matrix(blue_colors, ncol=1)), 
              xleft = 0.2, ybottom = 0.1, xright = 0.5, ytop = 0.45)
  text(x = 0.55, y = 0.45, labels = "High", pos = 4, cex = 0.2)
  text(x = 0.55, y = 0.1, labels = "Low", pos = 4, cex = 0.2)
  text(x = 0.5, y = 0.5, labels = paste('PMD',pmd_values,'High'), cex = 0.2, adj = 0.5)
  
  rasterImage(as.raster(matrix(red_colors, ncol=1)), 
              xleft =.2, ybottom = 0.55, xright = 0.5, ytop = 0.9)
  text(x = 0.55, y = 0.9, labels = "High", pos = 4, cex = 0.2)
  text(x = 0.55, y = 0.55, labels = "Low", pos = 4, cex = 0.2)
  text(x = 0.5, y = 0.95, labels = paste('PMD',pmd_values,'low'), cex = 0.2, adj = 0.5)
  dev.off()
}

#'@title Generate Molecular Network
#'@description This function generates a molecular network.
#'@param peak_file The input peak csv file from quantitative peaks list.
#'@param ion_cluster_file The ion cluster csv file.
#'@param annotation_file The annotation csv file from qualitative peaks list.
#'@param output_file The output file for the molecular network.
#'@return A molecular network plot and a csv file for the network.
generate_molecular_network <- function(peak_file, ion_cluster_file, annotation_file, output_file) {
  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]
  
  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))
  
  ioncluster <- fread(ion_cluster_file)
  anno <- fread(annotation_file)
  
  data("sda")
  hfpmd <- round(sda$PMD,3)
  df <- getpmddf(mz, group = ioncluster$class, pmd = hfpmd, digits = 3)
  
  df$anno1 <- anno$name[match(df$ms1, anno$mz)]
  df$anno2 <- anno$name[match(df$ms2, anno$mz)]
  
  dfanno <- df[complete.cases(df), ]
  
  # in untargeted data it is common that no PMD pair has both ends annotated;
  # writing an empty network file instead of failing on graph construction
  if (nrow(dfanno) == 0) {
    message("No PMD pair with both peaks annotated; the network is empty.")
    fwrite(dfanno, output_file)
    return(invisible(NULL))
  }
  
  net <- graph_from_data_frame(dfanno)
  graph <- as_tbl_graph(net)
  p <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(color = net)) + theme_void()
  print(p)
  fwrite(dfanno, output_file)
}