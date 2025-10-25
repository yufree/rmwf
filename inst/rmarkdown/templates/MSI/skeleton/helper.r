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
library(enviGCMS)
library(shinyWidgets)
library(shinydashboard)
library(visNetwork)

sourceCpp("one_over_k0_to_ccs.cpp")
sourceCpp('peakalign.cpp')
sourceCpp("one_over_k0_to_ccs.cpp")
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

#'@title Get Reference Peaks
#'@description This function processes raw mass spectrometry data to generate a reference peak list.
#'@param path The path to the raw data (e.g., .d folder).
#'@param accept_Bruker_EULA_and_on_Windows_or_Linux Set to TRUE to accept the Bruker EULA.
#'@param libpath The path to the Bruker library files.
#'@param batch_size The number of frames to process in each batch.
#'@param outref The output file for the reference peaks.
#'@param outcoord The output file for the coordinates.
#'@param xrange The x-range of the region of interest.
#'@param yrange The y-range of the region of interest.
#'@return A reference peak list.
getrefpeak <- function(path,accept_Bruker_EULA_and_on_Windows_or_Linux,libpath, batch_size, outref, outcoord, xrange,yrange){
  # download dll/so file and set the column to be collected
  if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
    folder_to_stode_priopriatary_code = libpath
    path_to_bruker_dll = download_bruker_proprietary_code(folder_to_stode_priopriatary_code)
    setup_bruker_so(path_to_bruker_dll)
    all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
  }
  # Build the connection with raw data
  D = OpenTIMS(path)
  # extract frame information or define a range to subset the data
  xy <- table2df(D, "MaldiFrameInfo")
  # subset coords idx for region of interests
  idx <- which(xy$MaldiFrameInfo$XIndexPos>xrange[1]&xy$MaldiFrameInfo$XIndexPos<xrange[2]&xy$MaldiFrameInfo$YIndexPos>yrange[1]&xy$MaldiFrameInfo$YIndexPos<yrange[2])
  # set a batch size to avoid memory issues by process the data in batches
  total_frames <- dim(xy$MaldiFrameInfo[idx,])[1]
  frame <- xy$MaldiFrameInfo$Frame[idx]
  
  # Process data in batches
  
  total_batches <- ceiling(total_frames / batch_size)
  
  dtx_list <- list()
  for (i in 1:total_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, total_frames)
    idx <- start_idx:end_idx
    
    data <- query(D, frames = frame[idx], columns = all_columns)
    data <- as.data.table(data)
    # round up to make sure binning is working and fit the mass/ion mobility accuracy
    data[, mz := round(mz, 4)]
    data[, inv_ion_mobility := round(inv_ion_mobility, 3)]
    
    dt <- data[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]
    dtx_list[[i]] <- dt  
    rm(data)
    gc()
  }
  # reshape the data by mz and ion mobility
  dtx <- rbindlist(dtx_list)
  dtx <- dtx[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]
  
  # add ccs value for super pixel peaks
  # change to libtimsdata.dll for windows
  lib_path <- paste0(libpath,"libtimsdata.so")
  result <- one_over_k0_to_ccs_parallel(dtx$inv_ion_mobility, rep(1,length(dtx$inv_ion_mobility)), dtx$mz, lib_path)
  dtx[,ccs:=result]
  # save super pixel peaks
  fwrite(dtx,outref)
  # save the location of pixels
  coord <- xy$MaldiFrameInfo[-idx,c('Frame','XIndexPos','YIndexPos')]
  fwrite(coord,outcoord)
}

#'@title Get Quantitative Peak List
#'@description This function generates a quantitative peak list.
#'@param refpath The path to the reference peaks file.
#'@param lib_path The path to the Bruker library files.
#'@param path The path to the raw data.
#'@param method The normalization method ('tic' or 'rms').
#'@param zero_proportion_cutoff The cutoff for removing peaks with a high proportion of zeros.
#'@param coordpath The path to the coordinates file.
#'@param normpath The output file for the normalized peak list.
#'@param xrange The x-range of the region of interest.
#'@param yrange The y-range of the region of interest.
#'@return A quantitative peak list.
getqlist <- function(refpath,lib_path,path, method, zero_proportion_cutoff, coordpath, normpath, xrange, yrange){
  ref <- fread(refpath)
  setup_bruker_so(lib_path)
  all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
  D = OpenTIMS(path)
  xy <- table2df(D, "MaldiFrameInfo")
  idx <- which(xy$MaldiFrameInfo$XIndexPos>xrange[1]&xy$MaldiFrameInfo$XIndexPos<xrange[2]&xy$MaldiFrameInfo$YIndexPos>yrange[1]&xy$MaldiFrameInfo$YIndexPos<yrange[2])
  # set a batch size to avoid memory issues by process the data in batches
  total_frames <- dim(xy$MaldiFrameInfo[idx,])[1]
  frame <- xy$MaldiFrameInfo$Frame[idx]
  # this setting will make sure the data.frame's row number is within the limits
  batch <- floor(2e+09/nrow(ref))
  batchs <- ceiling(total_frames/batch)
  result <- list()
  for(i in 1:batchs){
    upper <- min(total_frames,i*batch)
    idx<- c(((i-1)*batch+1):upper)
    data <- query(D, frames = frame[idx], columns = all_columns)
    data <- as.data.table(data)
    data[,ccs:=one_over_k0_to_ccs_parallel(data$inv_ion_mobility, rep(1,length(data$inv_ion_mobility)), data$mz, lib_path)]
    resultt <- findpeakalign(data$mz,data$ccs,data$intensity,data$frame,ref$mz,ref$ccs, ppm = 20, ccs_shift = 0.03)
    resultt <- as.data.table(resultt)
    if(method == 'tic'){
      # perform TIC normalization
      resultt[, TIC := sum(intensity, na.rm = TRUE), by = frame]
      resultt[, normalized_intensity := as.numeric(intensity / TIC * .N), by = frame]
      resultt[, TIC := NULL]
    }else{
      # perform RMS normalization
      resultt[, rms := sqrt(mean(sum(intensity^2, na.rm = TRUE))), by = frame]
      resultt[, normalized_intensity := as.numeric(intensity / rms), by = frame]
      resultt[, rms := NULL]
    }
    result[[i]] <- dcast(
      resultt, 
      frame ~ mz_ccs, 
      value.var = "normalized_intensity", 
      fun.aggregate = sum,
      fill = 0)
  }
  result <- rbindlist(result,fill=T)
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
#'@param ppm The m/z tolerance in ppm. Default is 10.
#'@param deltaccs The CCS tolerance. Default is 5.
#'@return An annotated peak list.
getanno <- function(database,mode,peakpath,annofile, ppm = 10, deltaccs = 5){
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
  dt[, np := rowSums(.SD != 0)]
  
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
  dt_values <- dt[, -1, with = FALSE]
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))
  
  mat <- t(dt_values)
  mat_centered <- scale(t(mat), center = TRUE, scale = FALSE)
  svd_result <- irlba(t(mat_centered), nv = n_components)
  pca_scores <- t(mat_centered) %*% svd_result$v
  
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
  
  Matrix <- t(subdt)
  row_norms <- sqrt(rowSums(Matrix^2))
  sim <- sweep(Matrix, 1, row_norms, "/")
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)
  
  t <- hclust(D_sim)
  s <- cutree(t, h = hclust_cutoff)
  
  name <- as.numeric(names(table(s)[table(s) > min_cluster_size]))
  
  Matrix <- t(subdt)
  split_matrices <- lapply(name, function(category) {
    rows <- which(s == as.numeric(category))
    subset_matrix <- Matrix[rows, , drop = FALSE]
    return(subset_matrix)
  })
  
  summed_matrices <- lapply(split_matrices, function(subset_matrix) {
    xxx <- apply(subset_matrix, 1, scale)
    x <- rowSums(xxx) / ncol(xxx)
    return(x)
  })
  
  result_matrix <- do.call(cbind, summed_matrices)
  
  clpan <- cbind.data.frame(x = x, y = y, result_matrix)
  
  dir.create(folder_name)
  for (i in c(1:length(name))) {
    dfx <- result_matrix[, i]
    norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx))
    color_palette <- colorRamp(c("skyblue", "red"))
    color_sequence <- rgb(color_palette(norm) / 255, alpha = 1)
    xlim <- max(x) - min(x) + 1
    ylim <- max(y) - min(y) + 1
    
    png(paste0(folder_name,'/cluster', name[i], '.png'), width = xlim, height = ylim)
    plot.new()
    par(mar = c(0, 0, 0, 0))
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
    dfx <- dt_filtered[, i]
    norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx))
    png(paste0(colnames(dt_filtered)[i], '.png'), width = width, height = height)
    plot.new()
    par(mar = c(0, 0, 0, 0))
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
  width = diff(range(xy$x))
  height = diff(range(xy$y))
  png(filename,width*2, height*2,bg = 'transparent',res = 300)
  layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(0.8, 0.2))
  plot.new()
  par(mar=c(0,0,0,0))
  plot.window(xlim = c(0,width), ylim = c(0,height), xaxs = "i", yaxs = "i", asp = NA)
  dfx <- xy[,..lowname][[1]]
  norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx)) *0.5
  # point_colors <- colors_with_alpha[col_index]
  points(xy$x-min(xy$x)+1, xy$y-min(xy$y)+1, pch = 16, col = rgb(1, 0, 0, norm),cex=0.3)
  dfx <- xy[,..highname][[1]]
  norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx)) *0.25
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
  
  net <- graph_from_data_frame(dfanno)
  graph <- as_tbl_graph(net)
  p <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(color = net)) + theme_void()
  print(p)
  fwrite(dfanno, output_file)
}