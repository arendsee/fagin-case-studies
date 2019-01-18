library(fagin)
library(rmonad)
library(knitr)
library(magrittr)

source("io.R")
source("run-fagin.R")
source("venn.R")

dir.create("output", showWarnings=FALSE)

con <- get_yeast_config()
m <- .run_fagin(con)

# a focal to deeper ordering for species (used in plots and tables)
species_order <- get_species_phylogenetic_order(con)

textab <- function(x, basename, digits=5, align="r", ...){
  knitr::kable(x, format="latex", digits=digits, align=strsplit(align, "")[[1]], ...) %>%
    gsub(pattern="Saccharomyces\\\\_([a-z]+)", replacement="\\\\textit{S. \\1}", perl=TRUE) %>%
    write(file=out(paste0(basename, ".tex")))
}

missues(m) %>% kable(format='html') %>% write(out('issues.html'))
mtabulate(m, code=TRUE) %>% kable(format='html') %>% write(out('table.html'))

pdf(out('rmonad.pdf'))
plot(m, vertex.size=2, edge.arrow.size=0.2, vertex.label=NA)
dev.off()

initial_protein_residue_counts(m) %>% textab("initial-aa")

final_protein_residue_counts(m) %>% textab("final-aa")
genomic_composition(m) %>% textab("genome-comp")
make_synder_flag_table(m) %>% textab("synder-flags")

## Input data summary

# Get a list of input summaries
ss <- get_summaries(m)

pdf(out("yeast-tree.pdf"))
plot(ape::read.tree(con@input@tree), show.node.label=TRUE)
dev.off()

strata <- readr::read_tsv("arendsee/fagin-yeast/archive/phylostratr-strata.tab")
dplyr::group_by(strata, mrca) %>%
  dplyr::summarize(count = n()) %>%
  textab("phylostratr-counts", caption="Number of genes in each phylostratum relative to the focal species *S. cerevisiae*.")

make_genome_table(m, species_order) %>% textab(
  basename="genome-table",
  align="lrrrrrrrr",
  caption="**Genomic statistics**. *scafs* - number of scaffolds in the genome. *bases* - total number of bases in the genome. *prots* - total number of protein coding models. *GC* - the percent of G and C nucleotides in the genome. *oddstart* - number of coding sequences that do not begin with M. *oddstop* - number of coding sequences that do not end with a stop (\\*). *p_N* - number of unknown nucleotides (N) in the genome. *n_stop* - number of coding sequences with an internal stop."
)

textab(
  make_synmap_table(m),
  basename="synmap-table",
  caption="Lengths of homologous links in each synteny map."
)

## `synder` output

textab(
    make_synder_table(m),
    basename="synder-table",
    caption="Search interval widths"
)

textab(
  get_synder_count_table(m),
  basename="synder-count-table",
  caption="Number of search intervals found for each focal gene against each target genome"
)

genes <- fagin:::.extract(m, 'query_genes')[[1]]
textab(
  get_synder_count_table(m, genes),
  basename="synder-gene-count-table",
  caption="Number of search intervals found for each query focal gene against each target genome"
)

## Fagin output

pdf(out("secondary-labels-1.pdf"))
plot_secondary_labels(m, species_order, fill='secondary')
dev.off()

pdf(out("secondary-labels-2.pdf"))
plot_secondary_labels(m, species_order, fill='species')
dev.off()

## Origin classifications
textab(make_origin_table(m), "origin-table")

pdf(out("venn.pdf"))
stuff <- get_yeast_uno_comp(m, strata)
par(mfrow=c(2,3))
for(n in names(stuff$comparisons)){
    limma::vennDiagram(stuff$venn(n), cex=1, main=n, asp=1)
}
dev.off()
