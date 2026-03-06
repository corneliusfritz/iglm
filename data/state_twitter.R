# Dynamically generate the state_twitter R6 data object on load
path <- system.file("extdata", "state_twitter_raw.rds", package = "iglm")
if (path == "") {
  if (dir.exists("inst")) {
    path <- file.path("inst", "extdata", "state_twitter_raw.rds")
  } else {
    path <- file.path("..", "inst", "extdata", "state_twitter_raw.rds")
  }
}
twit_raw <- readRDS(path)
if (exists("iglm.data", mode = "function")) {
  state_twitter <- list(
    iglm.data = get("iglm.data")(
      x_attribute = twit_raw$iglm.data$x_attribute,
      y_attribute = twit_raw$iglm.data$y_attribute,
      z_network = twit_raw$iglm.data$z_network,
      n_actor = twit_raw$iglm.data$n_actor,
      neighborhood = twit_raw$iglm.data$neighborhood,
      directed = twit_raw$iglm.data$directed,
      type_x = twit_raw$iglm.data$type_x,
      type_y = twit_raw$iglm.data$type_y,
      scale_x = twit_raw$iglm.data$scale_x,
      scale_y = twit_raw$iglm.data$scale_y,
      fix_x = twit_raw$iglm.data$fix_x,
      fix_z = twit_raw$iglm.data$fix_z,
      fix_z_alocal = twit_raw$iglm.data$fix_z_alocal
    ),
    match_gender = twit_raw$match_gender,
    match_race = twit_raw$match_race,
    match_state = twit_raw$match_state,
    white_attribute = twit_raw$white_attribute,
    gender_attribute = twit_raw$gender_attribute
  )
} else {
  state_twitter <- list() # Dummy object for roxygen2 parsing
}
rm(twit_raw, path)
