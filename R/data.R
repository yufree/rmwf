#' Demo mzrt list object
#' @docType data
#' @format mzrt list object
"mzrt"

#' Demo mzrt list object for negative mode
#' @docType data
#' @format mzrt list object
"mzrtn"

#' Demo ipo parameters object
#' @docType data
#' @format ipo parameters
"para"

#' Demo ipo parameters object for negative mode
#' @docType data
#' @format ipo parameters
"paran"

#' Demo xcmseic object
#' @docType data
#' @format xcmseic object
"srmeic"

#' Demo xcmseic object for negative mode
#' @docType data
#' @format xcmseic object
"srmneic"

#' Demo xcmsset object
#' @docType data
#' @format xcmsset object
"srmxset"

#' Demo xcmsset object for negative mode
#' @docType data
#' @format xcmsset object
"srmnxset"

#' A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annotation
#' @docType data
#' @format A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annotation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"qqq"

#' A list containing HMDB qtof MS/MS data with peaks larger than 10 percentage for PMD annotation
#' @docType data
#' @format A list containing HMDB qtof MS/MS data with peaks larger than 10 percentage for PMD annotation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"qtof"

#' A list containing HMDB orbitrap MS/MS data with peaks larger than 10 percentage for PMD annotation
#' @docType data
#' @format A list containing HMDB orbitrap MS/MS data with peaks larger than 10 percentage for PMD annotation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"orb"

#' A data frame with compounds from hmdb and refmet
#' @docType data
#' @format A data frame with compounds from hmdb and refmet including name, InChIKey, chemical_formula, data source and exact mass
#' \describe{
#'   \item{name}{compounds name}
#'   \item{InChIKey}{InChIKey}
#'   \item{chemical_formula}{chemical formula}
#'   \item{db}{database sources}
#'   \item{mass}{exact mass}
#'   }
"hr"

#' A list containing HMDB GC-EI-MS spectra database
#' @docType data
#' @format A list with compounds from hmdb for GC-MS simulation
"hmdbcms"

#' A list containing MoNA LC-MS spectra database
#' @docType data
#' @format A list with compounds from MoNA LC-MS spectra database for GC-MS simulation
"monams1"
