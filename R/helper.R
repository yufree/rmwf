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
