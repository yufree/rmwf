#' Shiny application for Mass Spectrometry Imaging Network Analysis
#' @export
runmsinet <- function() {
    file <- system.file("msinet",
                        package = "rmwf")
    if (file == "") {
        stop("Could not find directory. Try re-installing `enviGCMS`.",
             call. = FALSE)
    }
    shiny::runApp(file)
}
