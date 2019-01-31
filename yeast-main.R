library(fagin)
library(rmonad)
library(knitr)
library(magrittr)
library(readr)

source("run-fagin.R")
source("yeast-venn.R")
source("common.R")

out <- make_out("yeast-output")
con <- get_yeast_config()
m   <- .run_fagin(con)

strata <- readr::read_tsv("arendsee/fagin-yeast/archive/phylostratr-strata.tab")
common_stuff(m=m, con=con, out=out, strata=strata)

dplyr::group_by(strata, mrca) %>%
  dplyr::summarize(count = n()) %>%
  textab("phylostratr-counts", caption="Number of genes in each phylostratum relative to the focal species *S. cerevisiae*.")

pdf(out("venn.pdf"))
stuff <- get_yeast_uno_comp(m, strata)
par(mfrow=c(2,3))
for(n in names(stuff$comparisons)){
    limma::vennDiagram(stuff$venn(n), cex=1, main=n, asp=1)
}
dev.off()
