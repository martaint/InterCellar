#' Run the Shiny Application
#'
#' @param reproducible boolean for setting a seed, making plots reproducible
#' @return a running instance of InterCellar
#'
#' @export
#' @examples 
#' \dontrun{
#' run_app()
#' }
#' 
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  reproducible = TRUE
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui, 
      server = app_server,
      options = list("launch.browser" = TRUE)
    ), 
    golem_opts = list(reproducible = reproducible)
  )
}
