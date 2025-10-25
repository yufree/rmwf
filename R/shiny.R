#' Shiny application for Mass Spectrometry Imaging Network Analysis
#' @export
runmsinet <- function() {
    file <- system.file("msinet",
                        package = "rmwf")
    if (file == "") {
        stop("Could not find directory. Try re-installing `rmwf`.",
             call. = FALSE)
    }
    shiny::runApp(file)
}
#' Shiny application for Mass Spectrometry Imaging annotation with mz and ccs
#' @export
runmzccsanno <- function() {
    file <- system.file("mzccsanno",
                        package = "rmwf")
    if (file == "") {
        stop("Could not find directory. Try re-installing `rmwf`.",
             call. = FALSE)
    }
    shiny::runApp(file)
}
