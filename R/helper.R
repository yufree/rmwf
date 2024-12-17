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
