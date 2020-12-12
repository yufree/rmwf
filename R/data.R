#' Demo mzrt list object
#' @format mzrt list object
#' @usage data(mzrt)
"mzrt"

#' Demo mzrt list object for negative mode
#' @format mzrt list object
#' @usage data(mzrtn)
"mzrtn"

#' Demo ipo parameters object
#' @format ipo parameters
#' @usage data(para)
"para"

#' Demo ipo parameters object for negative mode
#' @format ipo parameters
#' @usage data(paran)
"paran"

#' Demo xcmseic object
#' @format xcmseic object
#' @usage data(srmeic)
"srmeic"

#' Demo xcmseic object for negative mode
#' @format xcmseic object
#' @usage data(srmneic)
"srmneic"

#' Demo xcmsset object
#' @format xcmsset object
#' @usage data(srmxset)
"srmxset"

#' Demo xcmsset object for negative mode
#' @format xcmsset object
#' @usage data(srmnxset)
"srmnxset"

#' A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annatation
#' @docType data
#' @usage data(qqq)
#' @format A list containing HMDB qqq MS/MS data with peaks larger than 10 percentage for PMD annatation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"qqq"

#' A list containing HMDB qtof MS/MS data with peaks larger than 10 percentage for PMD annatation
#' @docType data
#' @usage data(qtof)
#' @format A list containing HMDB qtof MS/MS data with peaks larger than 10 percentage for PMD annatation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"qtof"

#' A list containing HMDB orbitrap MS/MS data with peaks larger than 10 percentage for PMD annatation
#' @docType data
#' @usage data(orb)
#' @format A list containing HMDB orbitrap MS/MS data with peaks larger than 10 percentage for PMD annatation
#' \describe{
#'   \item{name}{HMDB ID}
#'   \item{mz}{mass to charge ratio}
#'   \item{msms}{msms pmd}
#'   \item{msmsraw}{raw msms data}
#'   }
"orb"

#' A data frame with compounds from hmdb and refmet
#' @docType data
#' @usage data(hr)
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
#' @usage data(hmdbcms)
#' @format A list with compounds from hmdb for GC-MS simulation
"hmdbcms"

#' A list containing MoNA LC-MS spectra database
#' @docType data
#' @usage data(monams1)
#' @format A list with compounds from MoNA LC-MS spectra database for GC-MS simulation
"monams1"
