# Common header: activate the project renv library and load shared utilities.
# Stage scripts start with: source(file.path(Sys.getenv("PM_PROJECT_ROOT", "."), "03_scripts", "R", "header.R"))
local({
  p <- Sys.getenv("PM_PROJECT_ROOT", "")
  if (!nzchar(p)) p <- normalizePath(".")
  Sys.setenv(PM_PROJECT_ROOT = p, RENV_PROJECT = p)
  if (!"renv" %in% loadedNamespaces()) source(file.path(p, "renv", "activate.R"), local = FALSE)
  source(file.path(p, "03_scripts", "R", "utils.R"), local = FALSE)
})
