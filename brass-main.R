library(fagin)
library(rmonad)
library(knitr)
library(magrittr)

source("run-fagin.R")
source("common.R")

out <- make_out("brass-output")
con <- get_brassicaceae_config()
m   <- .run_fagin(con)

common_stuff(m=m, con=con, out=out)
