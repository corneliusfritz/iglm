# Generate the state_twitter R6 data object for data/
twit_raw <- readRDS("inst/extdata/state_twitter_raw.rds")

state_twitter <- list(
  iglm.data = do.call(iglm.data, twit_raw$iglm.data),
  match_gender = twit_raw$match_gender,
  match_race = twit_raw$match_race,
  match_state = twit_raw$match_state,
  white_attribute = twit_raw$white_attribute,
  gender_attribute = twit_raw$gender_attribute
)

save(state_twitter, file = "data/state_twitter.rda", compress = "xz")
rm(twit_raw)
