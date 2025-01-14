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
