library(fagin)
library(rmonad)
library(knitr)
library(magrittr)

source("run-fagin.R")
source("common.R")
source("common-venn.R")
source("brass-venn.R")

out <- make_out("brass-output")
con <- get_brassicaceae_config()
m   <- .run_fagin(con)

strata <- readr::read_tsv("brassicaceae-strata.tab")
species_order <- rev(get_species_phylogenetic_order(con))
common_stuff(m=m, con=con, out=out, species_order=species_order, strata=strata)

pdf(out("brass-venn.pdf"))
stuff <- get_brass_uno_comp(m, strata)
par(mfrow=c(2,2))
for(n in names(stuff$comparisons)){
    limma::vennDiagram(stuff$venn(n), cex=1, main=n, asp=1)
}
dev.off()
