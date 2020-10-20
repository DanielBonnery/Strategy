#' Shiny app 1
#' @description launches shiny app 1
#' @export
shinyapp1 <- function() {
  appDir <- system.file("shinyapp1", "shinyapp1", package = "Strategy")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
