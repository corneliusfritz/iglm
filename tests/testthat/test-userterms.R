# 
# library(iglm)
# 
# test_that("Test some sufficient statistics for undirected networks", {
#   pkg_path <- file.path(tempdir(), "test_pkg")
#   on.exit({
#     if (dir.exists(pkg_path)) unlink(pkg_path, recursive = TRUE)
#   }, add = TRUE)
#   
#   create_userterms_skeleton(path = pkg_path, pkg_name = "iglm.custom")
#   old_wd <- getwd()
#   on.exit(setwd(old_wd), add = TRUE)
#   setwd(file.path(pkg_path, "iglm.custom"))
#   system("R CMD INSTALL .")
#   
#   library(iglm.custom)
#   data(state_twitter)
#   res <- statistics(state_twitter[[1]] ~ mutual + my_mutual)
#   expect_equal(as.numeric(res[1]), as.numeric(res[2]))
# })
