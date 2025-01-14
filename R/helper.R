#' Get the metabolomics workbench data
#' @param ID study ID of metabolomics workbench
#' @return list with details of checked study
#' @export
getmw <- function(ID){
    url <- paste0('http://www.metabolomicsworkbench.org/rest/study/study_id/',ID,'/mwtab')
    list <- jsonlite::fromJSON(url)
    return(list)
}

#' Convert metabolomics workbench data into mzrt list object
#' @param ID study ID of metabolomics workbench
#' @return list
#' @export
getmwlist <- function(ID){
    list <- getmw(ID)
    group <- list$SUBJECT_SAMPLE_FACTORS
    data <- list$MS_METABOLITE_DATA$Data[,-1]
    rownames(data) <- list$MS_METABOLITE_DATA$Data[,1]
    mz <- list$MS_METABOLITE_DATA$Metabolites$moverz_quant
    rt <- list$MS_METABOLITE_DATA$Metabolites$ri
    re <- list(data = data, group = group, mz=mz, rt=rt)
    return(re)
}

#' Convert metabolomics workbench factor into dataframe.
#' @param ID study ID of metabolomics workbench
#' @return dataframe
#' @export
getmwfactor <- function(ID){
    url <- paste0('http://www.metabolomicsworkbench.org/rest/study/study_id/',ID,'/factors')
    list <- jsonlite::fromJSON(url)
    df <- data.frame(matrix(unlist(list), nrow=length(list), byrow=T))
    return(df)
}
#' Get pmd details for specific reaction after the removal of isotopouge.
#' @param dt data.table of tims-MSI with first column as pixel ID and other columns for peaks with mz_ccs.
#' @param group mass to charge ratio group from either retention time or mass spectrometry imaging segmentation.
#' @param pmd a specific paired mass distance or a vector of pmds
#' @param digits mass or mass to charge ratio accuracy for pmd, default 2.
#' @param mdrange mass defect range to ignore. Default c(0.25,0.9) to retain the possible reaction related paired mass.
#' @return dataframe with paired peaks for specific pmd or pmds. When group is provided, a column named net will be generated to show if certain pmd will be local(within the same group) or global(across the groups)
#' @export
getpmddf2 <- function (dt, group = NULL, pmd = NULL, digits = 2, mdrange = c(0.25,
                                                                             0.9))
{
    mzccs <- colnames(dt)[-1]
    dt_values <- dt[, -1, with = FALSE]
    cor <-stats::cor(dt_values)

    mz <- sapply(strsplit(mzccs,'\\_'),function(x) round(as.numeric(x[1]),4))
    im <- sapply(strsplit(mzccs,'\\_'),function(x) as.numeric(x[2]))
    name <- paste0(mz,'_',round(im))
    dis <- stats::dist(mz, method = "manhattan")
    df <- cbind.data.frame(from = name[which(lower.tri(dis),
                                             arr.ind = TRUE)[, 1]], to = name[which(lower.tri(dis),
                                                                                    arr.ind = TRUE)[, 2]], ms1=mz[which(lower.tri(dis),
                                                                                                                        arr.ind = TRUE)[, 1]], ms2=mz[which(lower.tri(dis),
                                                                                                                                                            arr.ind = TRUE)[, 2]],
                           diff = as.numeric(dis), diff2 = round(as.numeric(dis),
                                                                 digits = digits), md = as.numeric(dis)%%1,cor = cor[lower.tri(cor)])
    if (!is.null(group)) {
        df$group1 <- group[match(df$ms1, mz)]
        df$group2 <- group[match(df$ms2, mz)]
    }
    else {
        df$group1 <- df$group2 <- rep(1, length(mz))
    }
    idx <- df$md < mdrange[1] | df$md > mdrange[2]
    df <- df[idx, ]
    isoindex <- (round(df$diff, digits) != 0) & ((df$diff%%1 <
                                                      0.01 & df$diff >= 1 & df$diff < 2) | (df$diff%%2 < 0.01 &
                                                                                                df$diff >= 2 & df$diff < 3) | (df$diff%%1 > 0.99 & df$diff >=
                                                                                                                                   1 & df$diff < 2) | (df$diff%%1 > 0.99 & df$diff >= 0 &
                                                                                                                                                           df$diff < 1))
    dfiso <- df[isoindex, ]
    dfdeiso <- df[!isoindex, ]
    if (!is.null(pmd)) {
        dfdeisopmd <- dfdeiso[dfdeiso$diff2 %in% pmd, ]
        dfdeisopmd$net <- ifelse(abs(dfdeisopmd$group1 - dfdeisopmd$group2) ==
                                     0, "local", "global")
        return(dfdeisopmd)
    }
    else {
        dfdeiso$net <- ifelse(abs(dfdeiso$group1 - dfdeiso$group2) ==
                                  0, "local", "global")
        return(dfdeiso)
    }
}
#' Convert an XCMSnExp object to an mzrt S3 object.
#'
#' @noRd
.XCMSnExp2mzrt <-
    function(XCMSnExp,
             method = "medret",
             value = "into",
             mzdigit = 4,
             rtdigit = 1)
    {
        data <- xcms::featureValues(XCMSnExp, value = value, missing = 0)
        group <-
            data.frame(apply(XCMSnExp@phenoData@data, 2, as.character),
                       stringsAsFactors = FALSE)

        # peaks info
        peaks <- xcms::featureDefinitions(XCMSnExp)
        mz <- peaks$mzmed
        rt <- peaks$rtmed
        mzrange <- peaks[, c("mzmin", "mzmax")]
        rtrange <- peaks[, c("rtmin", "rtmax")]
        rownames(data) <-
            paste0("M", round(mz, mzdigit), "T", round(rt, rtdigit))
        mzrt <-
            list(
                data = data,
                group = group,
                mz = mz,
                rt = rt,
                mzrange = mzrange,
                rtrange = rtrange
            )
        class(mzrt) <- "mzrt"
        return(mzrt)
    }
#' Convert an xcmsSet object to an mzrt S3 object.
#'
#' @noRd
.xcmsSet2mzrt <-
    function(xcmsSet,
             method = "medret",
             value = "into",
             mzdigit = 4,
             rtdigit = 1)
    {
        data <- xcms::groupval(xcmsSet, method = method,
                               value = value)
        rownames(data) <- xcms::groupnames(xcmsSet,
                                           template = paste0('M.',
                                                             10 ^ (mzdigit - 1),
                                                             'T.',
                                                             10 ^ (rtdigit - 1)))
        # peaks info
        peaks <- as.data.frame(xcms::groups(xcmsSet))
        mz <- peaks$mzmed
        rt <- peaks$rtmed
        if (dim(xcms::phenoData(xcmsSet))[2] > 1) {
            if (sum(grepl('sample_name', colnames(
                xcms::phenoData(xcmsSet)
            ))) > 0) {
                group <- xcms::phenoData(xcmsSet)
                colnames(data) <- group$sample_name
            } else{
                sample_name <- rownames(xcms::phenoData(xcmsSet))
                sample_group <-
                    xcms::phenoData(xcmsSet)
                group <-
                    cbind.data.frame(sample_name, sample_group)
            }
        } else{
            sample_name <- rownames(xcms::phenoData(xcmsSet))
            sample_group <-
                as.character(xcms::phenoData(xcmsSet)$class)
            group <-
                cbind.data.frame(sample_name, sample_group)
        }

        mzrange <- peaks[, c("mzmin", "mzmax")]
        rtrange <- peaks[, c("rtmin", "rtmax")]
        mzrt <-
            list(
                data = data,
                group = group,
                mz = mz,
                rt = rt,
                mzrange = mzrange,
                rtrange = rtrange
            )
        class(mzrt) <- "mzrt"
        return(mzrt)
    }
#' Get the mzrt profile and group information as a mzrt list and/or save them as csv or rds for further analysis.
#' @param xset xcmsSet/XCMSnExp objects
#' @param name file name for csv and/or eic file, default NULL
#' @param mzdigit m/z digits of row names of data frame, default 4
#' @param rtdigit retention time digits of row names of data frame, default 1
#' @param method parameter for groupval or featureDefinitions function, default medret
#' @param value parameter for groupval or featureDefinitions function, default into
#' @param eic logical, save xcmsSet and xcmsEIC objects for further investigation with the same name of files, you will need raw files in the same directory as defined in xcmsSet to extract the EIC based on the binned data. You could use `plot` to plot EIC for specific peaks. For example, `plot(xcmsEIC,xcmsSet,groupidx = 'M123.4567T278.9')` could show the EIC for certain peaks with m/z 206 and retention time 2789. default F
#' @param type csv format for further analysis, m means  Metaboanalyst, a means xMSannotator, p means Mummichog(NA values are imputed by `getimputation`, and F test is used here to generate stats and p value), o means full information csv (for `pmd` package), default o. mapo could output all those format files.
#' @return mzrt object, a list with mzrt profile and group information
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' getmzrt(xset, name = 'demo', type = 'mapo')
#' }
#' @seealso \code{\link{getdata}},\code{\link{getdata2}}, \code{\link{getdoe}}, \code{\link{getcsv}}, \code{\link{getfilter}}
#' @references
#' Smith, C.A., Want, E.J., O’Maille, G., Abagyan, R., Siuzdak, G., 2006. XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. Anal. Chem. 78, 779–787.
#' @export

getmzrt <-
    function(xset,
             name = NULL,
             mzdigit = 4,
             rtdigit = 1,
             method = "medret",
             value = "into",
             eic = FALSE,
             type = 'o') {
        if (inherits(xset, 'xcmsSet')) {
            if (eic) {
                eic <-
                    xcms::getEIC(xset,
                                 rt = "corrected",
                                 groupidx = seq_len(nrow(xset@groups)))
                eic@groupnames <-
                    xcms::groupnames(xset,
                                     template = paste0(
                                         'M.',
                                         10 ^ (mzdigit - 1),
                                         'T.',
                                         10 ^ (rtdigit - 1)
                                     ))
                saveRDS(eic, file = paste0(name, 'eic.rds'))
            } else {
                if (!is.null(name)) {
                    saveRDS(xset,
                            file = paste0(name, 'xset.rds'))
                }

            }
            result <-
                .xcmsSet2mzrt(
                    xset,
                    mzdigit = mzdigit,
                    rtdigit = rtdigit,
                    method = method,
                    value = value
                )
        }
        else if (inherits(xset, 'XCMSnExp')) {
            if (eic) {
                feature_chroms <-
                    xcms::featureChromatograms(xset, features = rep(T, length(
                        xcms::quantify(xset)@NAMES
                    )))
                saveRDS(feature_chroms,
                        file = paste0(name, 'eic.rds'))
            } else {
                if (!is.null(name)) {
                    xset2 <- methods::as(xset, "xcmsSet")
                    saveRDS(xset2,
                            file = paste0(name, 'xset.rds'))
                }
            }
            result <-
                .XCMSnExp2mzrt(
                    xset,
                    mzdigit = mzdigit,
                    rtdigit = rtdigit,
                    method = method,
                    value = value
                )
        }
        getcsv(
            list = result,
            name = name,
            mzdigit = mzdigit,
            rtdigit = rtdigit,
            type = type
        )
        return(result)
    }

#' Get xcmsset object in one step with optimized methods.
#' @param path the path to your data
#' @param index the index of the files
#' @param BPPARAM used for BiocParallel package
#' @param pmethod parameters used for different instrumentals such as 'hplcorbitrap', 'uplcorbitrap', 'hplcqtof', 'hplchqtof', 'uplcqtof', 'uplchqtof'. The parameters were from the reference
#' @param minfrac minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group, default 0.67
#' @param ... arguments for xcmsSet function
#' @details the parameters are extracted from the papers. If you use name other than the name above, you will use the default setting of XCMS. Also I suggest IPO packages or apLCMS packages to get reasonable data for your own instrumental. If you want to summit the results to a paper, remember to include those parameters.
#' @return a xcmsset object for that path or selected samples
#' @references Patti, G. J.; Tautenhahn, R.; Siuzdak, G. Nat. Protocols 2012, 7 (3), 508–516.
#' @examples
#' \dontrun{
#' library(faahKO)
#' cdfpath <- system.file('cdf', package = 'faahKO')
#' xset <- getdata(cdfpath, pmethod = ' ')
#' }
#' @seealso \code{\link{getdata2}}, \code{\link{getmzrt}}
#' @export
getdata <-
    function(path,
             index = FALSE,
             BPPARAM = BiocParallel::SnowParam(),
             pmethod = "hplcorbitrap",
             minfrac = 0.67,
             ...) {
        cdffiles <- list.files(path, recursive = TRUE, full.names = TRUE)
        if (index) {
            cdffiles <- cdffiles[index]
        }
        if (pmethod == "hplcorbitrap") {
            xset <- xcms::xcmsSet(
                cdffiles,
                BPPARAM = BPPARAM,
                method = "centWave",
                ppm = 2.5,
                peakwidth = c(10,
                              60),
                prefilter = c(3, 5000),
                ...
            )
            if (index & length(index) == 1) {
                xset3 <- xset
            } else {
                xset <- xcms::group(
                    xset,
                    bw = 5,
                    mzwid = 0.015,
                    minfrac = min
                )
                xset2 <- xcms::retcor(xset, "obiwarp")
                # you need group the peaks again for this corrected
                # data
                xset2 <-
                    xcms::group(
                        xset2,
                        bw = 5,
                        mzwid = 0.015,
                        minfrac = minfrac
                    )
                xset3 <-
                    xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
            }
        } else if (pmethod == "uplcorbitrap") {
            xset <- xcms::xcmsSet(
                cdffiles,
                BPPARAM = BPPARAM,
                method = "centWave",
                ppm = 2.5,
                peakwidth = c(5,
                              20),
                prefilter = c(3, 5000),
                ...
            )
            xset <- xcms::group(
                xset,
                bw = 2,
                mzwid = 0.015,
                minfrac = minfrac
            )
            xset2 <- xcms::retcor(xset, "obiwarp")
            # you need group the peaks again for this corrected
            # data
            xset2 <- xcms::group(
                xset2,
                bw = 2,
                mzwid = 0.015,
                minfrac = minfrac
            )
            xset3 <- xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
        } else if (pmethod == "hplcqtof") {
            xset <- xcms::xcmsSet(
                cdffiles,
                BPPARAM = BPPARAM,
                method = "centWave",
                ppm = 30,
                peakwidth = c(10,
                              60),
                prefilter = c(0, 0),
                ...
            )
            if (index & length(index) == 1) {
                xset3 <- xset
            } else {
                xset <- xcms::group(
                    xset,
                    bw = 5,
                    mzwid = 0.025,
                    minfrac = minfrac
                )
                xset2 <- xcms::retcor(xset, "obiwarp")
                # you need group the peaks again for this corrected
                # data
                xset2 <-
                    xcms::group(
                        xset2,
                        bw = 5,
                        mzwid = 0.025,
                        minfrac = minfrac
                    )
                xset3 <-
                    xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
            }
        } else if (pmethod == "hplchqtof") {
            xset <- xcms::xcmsSet(
                cdffiles,
                BPPARAM = BPPARAM,
                method = "centWave",
                ppm = 15,
                peakwidth = c(10,
                              60),
                prefilter = c(0, 0),
                ...
            )
            if (index & length(index) == 1) {
                xset3 <- xset
            } else {
                xset <- xcms::group(
                    xset,
                    bw = 5,
                    mzwid = 0.015,
                    minfrac = minfrac
                )
                xset2 <- xcms::retcor(xset, "obiwarp")
                # you need group the peaks again for this corrected
                # data
                xset2 <-
                    xcms::group(
                        xset2,
                        bw = 5,
                        mzwid = 0.015,
                        minfrac = minfrac
                    )
                xset3 <-
                    xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
            }
        } else if (pmethod == "uplcqtof") {
            xset <- xcms::xcmsSet(
                cdffiles,
                BPPARAM = BPPARAM,
                method = "centWave",
                ppm = 30,
                peakwidth = c(5,
                              20),
                prefilter = c(0, 0),
                ...
            )
            if (index & length(index) == 1) {
                xset3 <- xset
            } else {
                xset <- xcms::group(
                    xset,
                    bw = 2,
                    mzwid = 0.025,
                    minfrac = minfrac
                )
                xset2 <- xcms::retcor(xset, "obiwarp")
                # you need group the peaks again for this corrected
                # data
                xset2 <-
                    xcms::group(
                        xset2,
                        bw = 2,
                        mzwid = 0.025,
                        minfrac = minfrac
                    )
                xset3 <-
                    xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
            }
        } else if (pmethod == "uplchqtof") {
            xset <- xcms::xcmsSet(
                cdffiles,
                BPPARAM = BPPARAM,
                method = "centWave",
                ppm = 15,
                peakwidth = c(5,
                              20),
                prefilter = c(0, 0),
                ...
            )
            if (index & length(index) == 1) {
                xset3 <- xset
            } else {
                xset <- xcms::group(
                    xset,
                    bw = 2,
                    mzwid = 0.015,
                    minfrac = minfrac
                )
                xset2 <- xcms::retcor(xset, "obiwarp")
                # you need group the peaks again for this corrected
                # data
                xset2 <-
                    xcms::group(
                        xset2,
                        bw = 2,
                        mzwid = 0.015,
                        minfrac = minfrac
                    )
                xset3 <-
                    xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
            }
        } else {
            xset <- xcms::xcmsSet(cdffiles, BPPARAM = BPPARAM,
                                  ...)
            if (index & length(index) == 1) {
                xset3 <- xset
            } else {
                xset <- xcms::group(xset, minfrac = minfrac)
                xset2 <- xcms::retcor(xset, "obiwarp")
                # you need group the peaks again for this corrected
                # data
                xset2 <-
                    xcms::group(xset2, minfrac = minfrac)
                xset3 <-
                    xcms::fillPeaks(xset2, BPPARAM = BPPARAM)
            }
        }
        return(xset3)
    }
#' Get XCMSnExp object in one step from structured folder path for xcms 3.
#' @param path the path to your data
#' @param index the index of the files
#' @param snames sample names. By default the file name without extension is used
#' @param sclass sample classes.
#' @param phenoData data.frame or NAnnotatedDataFrame defining the sample names and classes and other sample related properties. If not provided, the argument sclass or the subdirectories in which the samples are stored will be used to specify sample grouping.
#' @param BPPARAM used for BiocParallel package
#' @param mode 'inMemory' or 'onDisk' see `?MSnbase::readMSData` for details, default 'onDisk'
#' @param ppp parameters for peaks picking, e.g. xcms::CentWaveParam()
#' @param rtp parameters for retention time correction, e.g. xcms::ObiwarpParam()
#' @param gpp parameters for peaks grouping, e.g. xcms::PeakDensityParam()
#' @param fpp parameters for peaks filling, e.g. xcms::FillChromPeaksParam(), PeakGroupsParam()
#' @details This is a wrap function for metabolomics data process for xcms 3.
#' @return a XCMSnExp object with processed data
#' @seealso \code{\link{getdata}},\code{\link{getmzrt}}
#' @export
getdata2 <- function(path,
                     index = FALSE,
                     snames = NULL,
                     sclass = NULL,
                     phenoData = NULL,
                     BPPARAM = BiocParallel::SnowParam(),
                     mode = "onDisk",
                     ppp = xcms::CentWaveParam(
                         ppm = 5,
                         peakwidth = c(5, 25),
                         prefilter = c(3, 5000)
                     ),
                     rtp = xcms::ObiwarpParam(binSize = 1),
                     gpp = xcms::PeakDensityParam(
                         sampleGroups = 1,
                         minFraction = 0.67,
                         bw = 2,
                         binSize = 0.025
                     ),
                     fpp = xcms::FillChromPeaksParam()) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE)
    if (index) {
        files <- files[index]
    }

    fromPaths <- xcms::phenoDataFromPaths(files)
    n <- dim(fromPaths)[2]
    sample_group <- NULL
    if (n > 1) {
        sample_group <- fromPaths[, 1]
        for (i in 2:n) {
            sample_group <- paste(sample_group, fromPaths[,
                                                          i], sep = "_")
        }
    } else {
        sample_group <- fromPaths[, 1]
    }
    sample_group <- data.frame(sample_group)

    if (is.null(snames)) {
        snames <- rownames(fromPaths)
    } else {
        rownames(sample_group) <- snames
    }

    pdata <- phenoData
    if (is.null(pdata)) {
        pdata <- sclass
        if (is.null(pdata))
            pdata <- methods::new("NAnnotatedDataFrame",
                                  sample_group)
    } else {
        if (inherits(pdata,"data.frame"))
            pdata <- methods::new("NAnnotatedDataFrame",
                                  sample_group)
        if (!inherits(pdata, "NAnnotatedDataFrame"))
            stop("phenoData has to be a data.frame or NAnnotatedDataFrame!")
    }
    raw_data <- MSnbase::readMSData(files, pdata = pdata,
                                    mode = mode)
    gpp@sampleGroups <- pdata$sample_group
    xod <- xcms::findChromPeaks(raw_data, param = ppp,
                                BPPARAM = BPPARAM)
    xod <- xcms::adjustRtime(xod, param = rtp)
    xod <- xcms::groupChromPeaks(xod, param = gpp)
    xod <-
        xcms::fillChromPeaks(xod, param = fpp, BPPARAM = BPPARAM)
    return(xod)
}
