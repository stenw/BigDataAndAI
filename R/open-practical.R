#' Open a practical
#' @export
open_practical <- function(name) {
  # Locate the file in the installed package
  template <- system.file("practicals", paste0(name, ".qmd"), package = "BigDataAndAI")
  
  if (template == ""){
    message("Package not installed. Looking in inst/practicals/ ...")
    template <- file.path("inst", "practicals", paste0(name, ".qmd"))
  } 
  
  # Copy it to the student's current working directory so they can edit it
  new_file <- file.path(getwd(), paste0(name, ".qmd"))
  file.copy(template, new_file)
  
  # Open it in RStudio
  if (rstudioapi::isAvailable()) file.edit(new_file)
}