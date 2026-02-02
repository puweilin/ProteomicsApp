#' Launch the Proteomics Analysis App
#'
#' @description
#' Starts the Proteomics Analysis Platform Shiny application.
#'
#' @param ... Additional arguments passed to \code{\link[shiny]{runApp}}
#'
#' @return This function does not return a value; it launches the Shiny app.
#'
#' @examples
#' \dontrun{
#' # Launch the app
#' run_app()
#'
#' # Launch with custom port
#' run_app(port = 8080)
#'
#' # Launch in browser
#' run_app(launch.browser = TRUE)
#' }
#'
#' @export
run_app <- function(...) {
  app_dir <- system.file("app", package = "ProteomicsApp")
  if (app_dir == "") {
    stop("Could not find app directory. Try re-installing `ProteomicsApp`.", call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
