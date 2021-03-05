pd <- NULL
np <- NULL
sts <- NULL
plt <- NULL
sns <- NULL
comhub_main <- NULL
genie_main <- NULL
genie_link_list <- NULL

.onLoad <- function(libname, pkgname){
  pd <<- reticulate::import(module = "pandas", delay_load = T)
  np <<- reticulate::import(module = "numpy", delay_load = T)

  path <- system.file("python", package = "ComhubbeR")
  comhub_module <- reticulate::import_from_path(module = "comhub", path = path)

  path <- system.file("python", package = "ComhubbeR")
  genie_module <- reticulate::import_from_path(module = "GENIE3", path = path)

  comhub_main <<- comhub_module$comhub
  genie_main <<- genie_module$GENIE3
  genie_link_list <<- genie_module$get_link_list
}

