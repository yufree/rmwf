#' Download files from a Figshare article
#'
#' This function downloads all files attached to a public Figshare article
#' using the Figshare API v2. It replaces the deprecated rfigshare package.
#'
#' @param article_id numeric, Figshare article ID
#' @param path character, destination directory (created if non-existent), default '.'
#' @return invisible character vector of downloaded file destination paths
#' @export
#' @examples
#' \dontrun{
#' download_figshare(13252196, "data/NIST1950")
#' }
download_figshare <- function(article_id, path = ".") {
  url <- paste0("http://api.figshare.com/v2/articles/", article_id)
  req <- curl::curl_fetch_memory(url)
  if (req$status_code != 200) {
    stop("Failed to fetch Figshare metadata for article ", article_id, " (HTTP ", req$status_code, ")")
  }
  files <- jsonlite::fromJSON(rawToChar(req$content))$files
  if (is.null(files) || !nrow(files)) {
    warning("No files found for Figshare article ", article_id)
    return(invisible(character(0)))
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  destfiles <- file.path(path, files$name)
  for (i in seq_len(nrow(files))) {
    if (!file.exists(destfiles[i])) {
      message("Downloading ", files$name[i], "...")
      download.file(files$download_url[i], destfiles[i], method = "libcurl", quiet = TRUE)
    }
  }
  invisible(destfiles)
}
