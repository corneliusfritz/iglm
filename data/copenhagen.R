# Dynamically generate the copenhagen R6 data object on load
path <- system.file("extdata", "copenhagen_raw.rds", package = "iglm")
if (path == "") {
    if (dir.exists("inst")) {
        path <- file.path("inst", "extdata", "copenhagen_raw.rds")
    } else {
        path <- file.path("..", "inst", "extdata", "copenhagen_raw.rds")
    }
}
copen_raw <- readRDS(path)
if (exists("iglm.data", mode = "function")) {
    copenhagen <- get("iglm.data")(
        x_attribute = copen_raw$x_attribute,
        y_attribute = copen_raw$y_attribute,
        z_network = copen_raw$z_network,
        n_actor = copen_raw$n_actor,
        neighborhood = copen_raw$neighborhood,
        directed = copen_raw$directed,
        type_x = copen_raw$type_x,
        type_y = copen_raw$type_y,
        scale_x = copen_raw$scale_x,
        scale_y = copen_raw$scale_y,
        fix_x = copen_raw$fix_x,
        fix_z = copen_raw$fix_z,
        fix_z_alocal = copen_raw$fix_z_alocal
    )
} else {
    copenhagen <- list() # Dummy object for roxygen2 parsing
}
rm(copen_raw)
