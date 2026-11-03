# Generate the copenhagen R6 data object for data/
library(iglm)
copen_raw <- readRDS("inst/extdata/copenhagen_raw.rds")

copenhagen <- do.call(iglm.data, copen_raw)

save(copenhagen, file = "data/copenhagen.rda", compress = "xz")
rm(copen_raw)
