#' Load the nhanes data
#'
#' @return A data frame containing the nhanes dataset.
#' @export
#' @examples
#' mod <- load_nhanes()
load_nhanes <- function() {
  path <- system.file("extdata", "mlnhanes.rds", package = "BigDataAndAI")
  if (path == "") stop("Model file not found. Try reinstalling the package.")
  readRDS(path)
}