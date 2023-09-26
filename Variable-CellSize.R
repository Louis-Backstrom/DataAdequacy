for (cell_size in c(2.0, 1.0, 0.5, 0.2)) {
  rm(list = setdiff(ls(), "cell_size"))
  gc()
  
  source("Base-Script.R")
}
