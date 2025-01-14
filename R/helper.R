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
#' Perform MS/MS dot product annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be list object from `getMSP`
#' @param ppm mass accuracy, default 10
#' @param prems precursor mass range, default 1.1 to include M+H or M-H
#' @param binstep bin step for consin similarity
#' @param consinc consin similarity cutoff for annotation. Default 0.6.
#' @return list with MSMS annotation results
#' @export
dotpanno <- function(file,
                     db = NULL,
                     ppm = 10,
                     prems = 1.1,
                     binstep = 1,
                     consinc = 0.6) {
    namemgf <- basename(file)
    sample <- MSnbase::readMgfData(file)
    prec <- MSnbase::precursorMz(sample)
    mz <- MSnbase::mz(sample)
    ins <- MSnbase::intensity(sample)
    idx <- unlist(ins) / max(unlist(ins)) > 0
    if (sum(idx) > 0) {
        mz <- unlist(mz)[idx]
        ins <- unlist(ins)[idx]
        ins <- ins / max(ins) * 100
        mz2bin <-
            cut(
                mz,
                breaks = seq(min(mz), max(mz) + binstep, binstep),
                include.lowest = TRUE
            )
        mz2bint <- stats::aggregate(ins, list(mz2bin), sum)
        colnames(mz2bint) <- c('factor', 'ins')
        if (inherits(db, 'list')) {
            range <- cbind(db$mz - prems, db$mz + prems)
            idxmz <-
                enviGCMS::getalign(range, matrix(
                    c(
                        prec - ppm / prec * 1e-06,
                        prec + ppm / prec * 1e-06
                    ),
                    nrow = 1
                ))
            idxmz <- idxmz[, 1]
            if (sum(idxmz) > 0) {
                name <- db$name[idxmz]
                mz2 <- db$mz[idxmz]
                msmsraw <- db$msmsraw[idxmz]

                ms2bindb <- function(x) {
                    ms2mz <-
                        cut(
                            x[, 1],
                            breaks = seq(
                                min(x[, 1]),
                                max(x[, 1]) + binstep,
                                binstep
                            ),
                            include.lowest = TRUE
                        )
                    ms2ins <-
                        stats::aggregate(x[, 2], list(ms2mz), sum)
                    colnames(ms2ins) <-
                        c('factor', 'insdb')

                    db <- merge(mz2bint, ms2ins)
                    cosscore <-
                        crossprod(db$ins, db$insdb) / sqrt(crossprod(db$ins) * crossprod(db$insdb))
                    return(cosscore)
                }
                ms2bin <-
                    vapply(msmsraw, ms2bindb, 1)

                t <-
                    list(
                        name = name[ms2bin > consinc & !is.na(ms2bin > consinc)],
                        mz = mz2[ms2bin > consinc &
                                     !is.na(ms2bin > consinc)],
                        msmsraw = msmsraw[ms2bin >
                                              consinc &
                                              !is.na(ms2bin > consinc)],
                        coss = ms2bin[ms2bin >
                                          consinc &
                                          !is.na(ms2bin > consinc)],
                        dmz = mz,
                        dprc = prec,
                        dins = ins,
                        file = namemgf
                    )
            } else{
                return(NULL)
            }
        } else{
            return(NULL)
        }
    } else{
        return(NULL)
    }
}

#' Perform MS/MS X rank annotation for mgf file
#' @param file mgf file generated from MS/MS data
#' @param db database could be list object from `getms2pmd`
#' @param ppm mass accuracy, default 10
#' @param prems precursor mass range, default 1.1 to include M+H or M-H
#' @param intc intensity cutoff for peaks. Default 0.1
#' @param quantile X rank quantiles cutoff for annotation. Default 0.75.
#' @return list with MSMS annotation results
#' @export
xrankanno <- function(file,
                      db = NULL,
                      ppm = 10,
                      prems = 1.1,
                      intc = 0.1,
                      quantile = 0.75) {
    # the score table is copied from MatchWeiz package https://github.com/AharoniLab/MatchWeiz
    xrankm <- data.frame(
        V2 = c(
            5.1141927,
            4.6402166,
            4.0347778,
            3.2659687,
            2.8868788,
            2.5778594,
            2.4243286,
            2.3812965,
            2.3287661,
            2.2854422,
            2.1004184,
            1.8994774,
            1.5654162,
            0.8235242,
            0.4100572,
            0.7430581,
            1.0656059,
            1.0032377,
            0.8241324,
            0.7564311,
            0.917284,
            1.0002008,
            1.0621707,
            1.1566258,
            1.1392231,
            0.9148275,
            1.1658015,
            1.3023251,
            1.3606814,
            1.3072061,
            -1.0277345
        ),
        V3 = c(
            4.4662252,
            4.2625276,
            4.0266953,
            3.7418072,
            3.3847624,
            3.1361971,
            2.9080739,
            2.713817,
            2.4905474,
            2.304285,
            2.049122,
            1.8832838,
            1.868606,
            1.9907963,
            2.1194155,
            1.9654197,
            1.5399935,
            0.9886019,
            0.5569938,
            0.8646485,
            1.2852703,
            1.2970452,
            1.1906578,
            1.0054851,
            0.4643638,
            0.3149796,
            0.6530265,
            0.9365217,
            1.276735,
            1.6331704,
            -0.8794154
        ),
        V4 = c(
            3.6795088,
            3.7669894,
            3.7782841,
            3.7152657,
            3.5313373,
            3.2815733,
            3.0549055,
            2.8764497,
            2.6851887,
            2.4204358,
            2.2089195,
            2.0982346,
            1.9607808,
            1.7542598,
            1.3818804,
            1.198171,
            1.246158,
            1.2886556,
            1.2656819,
            1.2415007,
            1.1999841,
            1.0194242,
            0.9042444,
            0.8420256,
            0.7771686,
            0.4544566,
            0.5260332,
            0.6261827,
            0.7254316,
            0.8462921,
            -0.7569588
        ),
        V5 = c(
            3.2839027,
            3.5150629,
            3.6093468,
            3.6011231,
            3.4368337,
            3.2726368,
            3.1239369,
            2.9851635,
            2.8044835,
            2.4863963,
            2.2902994,
            2.2353244,
            2.2524986,
            2.2618101,
            2.0346991,
            1.7207373,
            1.3948365,
            0.7286686,
            0.1800626,
            0.551818,
            0.9307692,
            1.2736242,
            1.3806847,
            1.2709348,
            1.0789728,
            1.0715978,
            1.3090899,
            1.3045615,
            1.1549414,
            0.5900862,
            -0.654933
        ),
        V6 = c(
            2.5984291,
            3.0046244,
            3.203698,
            3.2877079,
            3.2571621,
            3.1834096,
            3.085427,
            2.9219801,
            2.6760367,
            2.4617591,
            2.3534256,
            2.2377258,
            2.0590923,
            1.9655526,
            1.9060469,
            1.9008196,
            1.7993952,
            1.6434069,
            1.71971,
            1.64201,
            1.4450944,
            1.2230123,
            1.0367481,
            1.1504708,
            1.1010623,
            0.9585488,
            0.8015804,
            0.7383734,
            0.6980518,
            0.7457697,-0.5689874
        ),
        V7 = c(
            2.2602371,
            2.7464055,
            3.0054499,
            3.1477893,
            3.2077435,
            3.165352,
            3.0703202,
            2.9495569,
            2.8437625,
            2.7068907,
            2.3780462,
            2.0283327,
            1.9713898,
            2.0096547,
            1.9692381,
            1.7104057,
            1.3232825,
            1.2789279,
            1.3324883,
            1.3299554,
            1.3211272,
            1.2046809,
            1.1817414,
            0.8955127,
            0.7820397,
            1.0082414,
            0.8953015,
            0.8201049,
            0.7231856,
            0.6023505,
            -0.4953545
        ),
        V8 = c(
            1.780844,
            2.4251461,
            2.7362897,
            2.9017494,
            2.9748077,
            2.9993037,
            2.9341194,
            2.8141705,
            2.6106826,
            2.437433,
            2.3719587,
            2.2569385,
            2.0227376,
            1.7871484,
            1.7484059,
            1.7687429,
            1.6906459,
            1.3787828,
            1.1912703,
            1.1235638,
            0.8294797,
            0.40528,
            0.2093585,
            0.4069855,
            0.6013594,
            0.6819047,
            0.7574698,
            0.7689249,
            0.7047324,
            0.5285583,
            -0.4316488
        ),
        V9 = c(
            1.3428288,
            2.1247415,
            2.4586731,
            2.6210534,
            2.6626053,
            2.6148896,
            2.55033,
            2.4725361,
            2.387539,
            2.2699826,
            2.0674599,
            1.8867581,
            1.8696961,
            1.9600889,
            1.9703743,
            1.9428283,
            1.9346449,
            1.9797194,
            1.866236,
            1.5232692,
            1.3395643,
            1.3270395,
            1.2546527,
            1.2453368,
            0.9798871,
            0.8051503,
            0.8642886,
            0.7171807,
            0.3767922,-0.6804648,
            -0.376095
        ),
        V10 = c(
            1.9469214,
            2.0541428,
            2.2259246,
            2.4281479,
            2.6737424,
            2.8176345,
            2.7858621,
            2.6781424,
            2.6216533,
            2.6928871,
            2.6242811,
            2.2949485,
            1.8479744,
            1.3822823,
            1.4812091,
            1.861727,
            1.9892546,
            1.9760403,
            1.7386677,
            1.4067558,
            1.2361354,
            1.0888057,
            1.2080851,
            1.3518449,
            1.4029794,
            1.3922017,
            1.1391871,
            0.8259857,
            0.3274847,-0.8479935,
            -0.3285679
        ),
        V11 = c(
            1.032439,
            1.868163,
            2.227866,
            2.420178,
            2.49452,
            2.452937,
            2.511973,
            2.608769,
            2.682327,
            2.715736,
            2.574212,
            2.300455,
            2.043854,
            1.940104,
            2.079843,
            2.059508,
            1.800115,
            1.632611,
            1.504172,
            1.570284,
            1.713437,
            1.810089,
            1.908065,
            1.709064,
            1.385227,
            1.284118,
            1.257085,
            1.261826,
            1.320244,
            1.418463,
            -0.288381
        ),
        V12 = c(
            1.7018212,
            1.6708563,
            1.7280936,
            1.8510958,
            2.0413064,
            2.219501,
            2.300718,
            2.410956,
            2.5309432,
            2.642857,
            2.6473313,
            2.533732,
            2.3803898,
            2.1034979,
            1.8401677,
            1.7326038,
            1.6869479,
            1.6361263,
            1.422603,
            1.2904044,
            1.519895,
            1.781235,
            1.7759146,
            1.4436766,
            0.754851,
            0.2016948,
            0.8460149,
            1.0949137,
            1.2622762,
            1.3558464,-0.2535718
        ),
        V13 = c(
            1.0238163,
            1.3920828,
            1.6508524,
            1.8516509,
            2.0140125,
            2.0816779,
            2.1657678,
            2.3173347,
            2.5305612,
            2.646542,
            2.6022692,
            2.4213696,
            2.1419737,
            1.9741682,
            2.019168,
            2.0968997,
            2.1768472,
            2.1030137,
            1.8888521,
            1.6483642,
            1.5851545,
            1.5515542,
            1.3544017,
            1.1527845,
            1.012491,
            1.2094269,
            1.2503864,
            1.2568069,
            1.2503217,
            1.2125666,
            -0.2224962
        ),
        V14 = c(
            0.3468005,
            1.2964461,
            1.684207,
            1.8864016,
            1.9847797,
            1.9246173,
            1.7570264,
            1.5098705,
            1.423666,
            1.6670175,
            1.8653912,
            1.8688202,
            1.9000031,
            2.1181577,
            2.2916102,
            2.2332598,
            2.0661117,
            2.0224624,
            2.1251125,
            2.1607475,
            2.0161522,
            1.778285,
            1.6032933,
            1.4739095,
            1.257173,
            1.153034,
            1.0964124,
            1.0498028,
            1.1531411,
            1.3371823,
            -0.1989078
        ),
        V15 = c(
            0.0507587,
            0.6336722,
            0.9255884,
            1.0862229,
            1.1237131,
            1.3705725,
            1.6614447,
            1.7259731,
            1.7236453,
            1.7218323,
            1.7646869,
            2.0101464,
            2.2578926,
            2.3815446,
            2.2970573,
            2.0461818,
            1.8252976,
            1.7460836,
            1.445195,
            0.826319,
            0.8080851,
            1.4811161,
            1.8048276,
            1.7409854,
            1.299695,
            0.4519152,
            0.4793704,
            0.466986,
            0.5728771,
            0.744037,
            -0.1824328
        ),
        V16 = c(
            -0.4738215,
            0.4727348,
            0.9274847,
            1.2048411,
            1.421322,
            1.5537476,
            1.5433123,
            1.5296399,
            1.639142,
            1.6306439,
            1.6520935,
            1.9257464,
            2.0885409,
            2.0734827,
            1.8262528,
            1.5423633,
            1.5733343,
            1.7136354,
            1.7929254,
            1.8988632,
            1.8928769,
            1.5880022,
            1.0706523,
            0.6158382,
            0.5398145,
            0.7050322,
            0.6450264,
            0.795963,
            1.0421368,
            1.3426467,
            -0.1659159
        ),
        V17 = c(
            0.415412306,
            0.669041161,
            0.887952888,
            1.065570319,
            1.252917184,
            1.437232432,
            1.42398399,
            1.470966126,
            1.480313211,
            1.446783162,
            1.653113718,
            1.796406092,
            1.834872806,
            1.819048923,
            1.769825465,
            1.899545685,
            2.039529007,
            1.9965402,
            1.863095503,
            1.681937696,
            1.521084486,
            1.360576445,
            0.89712193,
            0.453765169,
            0.739650412,
            1.171870428,
            1.286366885,
            1.226049366,
            0.940129628,
            0.006144024,
            -0.148891361
        ),
        V18 = c(
            0.4029431,
            0.578631,
            0.6844446,
            0.7526526,
            0.6999427,
            0.8538586,
            1.3808023,
            1.7124098,
            1.8394157,
            1.8609021,
            1.7973889,
            1.7942522,
            1.7493591,
            1.6624452,
            1.8384434,
            1.9720332,
            1.9264542,
            1.7442466,
            1.6220334,
            1.7065845,
            1.698036,
            1.5573873,
            1.4871053,
            1.5647806,
            1.616974,
            1.4225419,
            1.1540925,
            1.0463909,
            1.1018232,
            1.3238546,
            -0.1346162
        ),
        V19 = c(
            -0.2086738,
            0.104658,
            0.4236391,
            0.7020921,
            0.9581844,
            1.1563411,
            1.4202728,
            1.5447069,
            1.5203407,
            1.3263225,
            1.4763387,
            1.9044777,
            2.1088163,
            2.1955835,
            2.1880875,
            2.0593745,
            1.868448,
            1.6182065,
            1.5225032,
            1.4718988,
            1.3159038,
            1.0773942,
            0.7954514,
            1.2300444,
            1.6228633,
            1.7196485,
            1.50843,
            1.2931858,
            0.9776897,
            0.5346483,
            -0.1219615
        ),
        V20 = c(
            -1.17466298,
            -0.57517985,
            -0.09469124,
            0.29341358,
            0.6275681,
            0.85657918,
            0.97005964,
            1.08544676,
            1.24737813,
            1.28432685,
            1.34063508,
            1.44902602,
            1.47152529,
            1.46396956,
            1.4450822,
            1.50800797,
            1.51509685,
            1.4236821,
            1.36516257,
            1.39389767,
            1.50727537,
            1.61656326,
            1.72162222,
            1.6186774,
            1.26994512,
            1.0449613,
            1.18170636,
            1.20705172,
            1.18974513,
            1.10240654,-0.10986603
        ),
        V21 = c(
            -0.64155598,
            0.05612773,
            0.46926615,
            0.74609668,
            0.9417987,
            1.10472549,
            1.21313115,
            1.32259961,
            1.3824699,
            1.16137716,
            1.14685416,
            1.28688178,
            1.48932069,
            1.65493148,
            1.64094851,
            1.64199637,
            1.61422595,
            1.52110189,
            1.44505869,
            1.29578841,
            1.28058613,
            1.4234758,
            1.61148215,
            1.60230881,
            1.37761036,
            1.07941075,
            1.23067113,
            1.39335228,
            1.59333376,
            1.81099522,-0.09853742
        ),
        V22 = c(
            -0.25237878,
            -0.14055822,
            -0.04412339,
            0.05458255,
            0.15600038,
            0.16282601,
            0.32542217,
            0.48810514,
            0.62898278,
            0.80684616,
            1.19067836,
            1.36729266,
            1.33640048,
            1.39873367,
            1.60794424,
            1.82212914,
            1.77806477,
            1.36875503,
            1.19615472,
            1.62992665,
            1.94805433,
            1.91215583,
            1.48434178,
            1.02365593,
            1.10686547,
            1.42563978,
            1.62560531,
            1.55871876,
            1.26887519,
            0.28447474,-0.08871415
        ),
        V23 = c(
            -0.16115903,
            -0.08391832,
            -0.01647319,
            0.02485493,
            0.05622159,
            0.08114541,
            0.08074582,
            0.286037,
            0.8045066,
            1.11176416,
            1.19684153,
            1.1567433,
            1.16076905,
            1.27902749,
            1.22616615,
            0.91111231,
            0.60854019,
            0.75664073,
            1.21895461,
            1.36165573,
            1.11617657,
            0.91607543,
            0.95890797,
            1.03957485,
            1.16212089,
            1.04901458,
            0.84476625,
            0.67451634,
            0.47356615,
            0.23627478,-0.08018362
        ),
        V24 = c(
            -0.71580503,
            -0.27217826,
            -0.0436726,
            0.05140871,
            0.06902796,
            0.04674192,
            0.1338606,
            0.29404185,
            0.67017027,
            0.97898454,
            1.07231188,
            1.15895578,
            1.07801359,
            0.91601943,
            0.94290241,
            0.87139595,
            1.0393812,
            1.35791832,
            1.41817197,
            1.40620617,
            1.23069528,
            1.12921573,
            1.26412864,
            1.27032549,
            1.23857962,
            1.20233683,
            1.12570276,
            1.14053571,
            1.23413259,
            1.3865318,-0.07274541
        ),
        V25 = c(
            -0.80447615,
            -0.36978346,
            -0.1149565,
            0.06442209,
            0.15423803,
            0.11491788,
            0.24894237,
            0.46719968,
            0.75581492,
            0.79729455,
            0.56349886,
            0.75343982,
            1.19469527,
            1.33759283,
            1.18629656,
            0.69901844,
            0.49800694,
            1.0295121,
            1.21341662,
            1.31957509,
            1.50155538,
            1.44061521,
            1.50586986,
            1.59421589,
            1.49538858,
            1.42624738,
            1.29586545,
            1.11294976,
            0.85942509,
            0.4968517,-0.06641275
        ),
        V26 = c(
            -0.88028928,
            -0.87742406,
            -0.81625187,
            -0.71987089,-0.56049376,
            -0.55713467,
            -0.351607,
            0.1028415,
            0.32873133,
            0.21449839,
            0.21092383,
            0.43334989,
            0.63600237,
            0.52970601,
            0.27485446,
            0.66284492,
            1.08874861,
            1.27593303,
            1.23844583,
            0.9993165,
            0.68768581,
            0.7926368,
            1.23239054,
            1.44573497,
            1.42903186,
            1.19000917,
            1.02234449,
            0.80780815,
            0.52961934,
            0.06918921,-0.06128711
        ),
        V27 = c(
            0.75839995,
            -1.3518103,
            -0.66685695,
            -0.35336782,-0.18595038,
            -0.2319309,
            -0.08208004,
            0.17125204,
            0.24066322,
            0.02703125,
            -0.3487187,-0.35423592,
            -0.1599313,
            0.13835646,
            0.38145214,
            0.79939512,
            1.1332018,
            1.08063781,
            0.94091869,
            0.89119766,
            1.00726349,
            1.16833776,
            1.02817672,
            0.78434039,
            0.5377862,
            0.32062451,
            0.30234168,
            0.51861785,
            0.87403661,
            1.25392001,-0.05731214
        ),
        V28 = c(
            -1.64425473,
            -1.46995102,
            -1.05072124,
            -0.60697974,-0.14356611,
            -0.04617828,
            -0.22441688,-0.23753935,
            0.0992326,
            0.49057842,
            0.54421444,
            0.184126,
            -0.46020092,
            -0.58709485,
            0.09349281,
            0.49746967,
            0.61976976,
            0.5339048,
            0.30912358,
            0.41219825,
            0.72645343,
            0.87199364,
            0.89251851,
            0.64601449,
            0.52924344,
            0.7819728,
            0.89644932,
            0.84901067,
            0.6284136,
            0.03421563,-0.05443636
        ),
        V29 = c(
            -0.91133014,
            -0.50186916,
            -0.36481741,
            -0.37947046,-0.58111115,
            -0.70460767,
            -0.18327545,
            0.16473325,
            0.33814655,
            0.49096046,
            0.89235074,
            1.23157558,
            1.18777596,
            0.97225027,
            0.67527971,
            0.50100392,
            0.66630452,
            0.77393421,
            1.18952284,
            1.46158781,
            1.37466942,
            1.1170623,
            0.83685894,
            0.90365669,
            1.08840309,
            0.98018994,
            0.88652635,
            0.96012355,
            1.12100016,
            1.34711566,-0.05258754
        ),
        V30 = c(
            -0.20610999,
            -0.27010751,
            -0.33343182,
            -0.39049313,-0.47370422,
            -0.45503989,
            -0.71344963,-1.15403172,
            -0.32208072,
            0.32604399,
            0.41668408,
            0.14504194,
            -0.03456606,
            0.35753346,
            0.75340612,
            0.67381566,
            0.446751,
            0.75319635,
            1.1025185,
            1.21055796,
            1.03444437,
            0.69610997,
            0.58632956,
            0.57997311,
            0.62885149,
            0.91179354,
            0.89000259,
            0.88535445,
            0.90948994,
            0.93495735,-0.05178215
        ),
        V31 = c(
            -3.10534277,
            -2.08074716,
            -1.58449918,
            -1.28784782,-1.0321364,
            -0.9420936,
            -1.38186672,-1.4291478,
            -0.7798967,
            -0.33427203,
            0.09787527,
            0.20863284,
            0.22152833,
            0.17018238,
            0.11081253,
            0.36453404,
            0.39386712,
            0.47080299,
            0.55998291,
            0.48994446,
            0.64023255,
            0.6359818,
            0.45718184,
            0.39041628,
            0.29089896,
            0.45555962,
            0.67865541,
            0.82889259,
            0.94950862,
            1.05494627,-0.05210723
        )
    )
    namemgf <- basename(file)
    sample <- MSnbase::readMgfData(file)
    prec <- MSnbase::precursorMz(sample)
    mz <- MSnbase::mz(sample)
    ins <- MSnbase::intensity(sample)
    idx <- unlist(ins) / max(unlist(ins)) > intc
    if (sum(idx) > 0) {
        mz <- unlist(mz)[idx]
        ins <- unlist(ins)[idx]
        ins <- ins / max(ins) * 100
        mzr <- mz[order(ins, decreasing = TRUE)]
        if (length(mzr) > 30) {
            mzr <- mzr[1:30]
        }
        if (inherits(db, 'list')) {
            range <- cbind(db$mz - prems, db$mz + prems)
            idxmz <-
                enviGCMS::getalign(range, matrix(
                    c(
                        prec - ppm / prec * 1e-06,
                        prec + ppm / prec * 1e-06
                    ),
                    nrow = 1
                ))
            idxmz <- idxmz[, 1]
            if (sum(idxmz) > 0) {
                name <- db$name[idxmz]
                mz2 <- db$mz[idxmz]
                msmsraw <- db$msmsraw[idxmz]

                ms2rank <- function(x) {
                    mz <- x[, 1]
                    ins <- x[, 2]
                    idx <-
                        (ins / max(ins)) > intc
                    mz <- mz[idx]
                    ins <- ins[idx]
                    mzrdb <-
                        mz[order(ins, decreasing = TRUE)]
                    if (length(mzrdb) > 30) {
                        mzrdb <- mzrdb[1:30]
                    }
                    score <- rep(NA, length(mzr))
                    for (i in seq_along(mzr)) {
                        match <- which.min(abs(mzrdb - mzr[i]) / mzr[i] * 1e06 < ppm)
                        if (length(match) > 0) {
                            match <-  match[which.min(
                                abs(
                                    mzrdb - mzr[i]
                                ) / mzr[i] * 1e06 < ppm
                            )]
                            score[i] <-
                                xrankm[match, i]
                        } else{
                            score[i] <- xrankm[31, i]
                        }
                    }
                    # reverse score
                    score2 <-
                        rep(NA, length(mzrdb))
                    for (i in seq_along(mzrdb)) {
                        match <- which.min(
                            abs(mzr - mzrdb[i]) / mzrdb[i] * 1e06 < ppm
                        )
                        if (length(match) > 0) {
                            match <-  match[which.min(
                                abs(
                                    mzr - mzrdb[i]
                                ) / mzrdb[i] * 1e06 < ppm
                            )]
                            score2[i] <-
                                xrankm[match, i]
                        } else{
                            score2[i] <- xrankm[31, i]
                        }
                    }

                    return(mean(c(
                        sum(score),
                        sum(score2)
                    )))
                }
                ms2xr <- vapply(msmsraw, ms2rank, 1)
                idx <-
                    ms2xr / max(ms2xr, na.rm = TRUE) > quantile

                t <-
                    list(
                        name = name[idx & !is.na(ms2xr)],
                        mz = mz2[idx &
                                     !is.na(ms2xr)],
                        msmsraw = msmsraw[idx &
                                              !is.na(ms2xr)],
                        score = ms2xr[idx &
                                          !is.na(ms2xr)],
                        dmz = mz,
                        dprc = prec,
                        dins = ins,
                        file = namemgf
                    )
            } else{
                return(NULL)
            }
        } else{
            return(NULL)
        }
    } else{
        return(NULL)
    }
}

#' Show MS/MS pmd annotation result
#' @param anno list from MSMS anno function
#' @param ... other parameter for plot function
#' @return NULL
#' @export
plotanno <- function(anno, ...) {
    for (i in seq_along(anno$name)) {
        graphics::plot(anno$msmsraw[[i]],
                       type = 'h',
                       main = anno$name[[i]],
                       ...)
        graphics::points(anno$dins ~ anno$dmz,
                         type = 'h',
                         col = 'red')
    }
}
