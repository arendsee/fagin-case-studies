# Count something as ORFic iff it matches an annotated protein
by_strict <- function(node_labels){
  cls <- node_labels$primary
  snd <- node_labels$secondary
  cls[snd == "O2" | snd == "O3"] <- "Non-ORFic"
  cls
}

make_venn_function <- function(strata, x, y){
  # set youngest strata name to "orphan"
  strata$mrca[strata$ps ==  max(strata$ps)] <- "orphan"
  phylostratr_ages <- split(strata, strata$mrca) %>%
      lapply(function(x) x$seqid)

  fagin_ages_x <- lapply(x$strOrigins, names)
  fagin_ages_y <- lapply(y$strOrigins, names)

  stopifnot(setequal(names(phylostratr_ages), names(fagin_ages_x)))
  stopifnot(setequal(names(phylostratr_ages), names(fagin_ages_y)))

  # venn maker function
  function(group){
    xs <- list(
        'standard' = phylostratr_ages[[group]]
      , 'default' = fagin_ages_x[[group]]
      , 'strict' = fagin_ages_y[[group]]
    )
    ids <- unique(unlist(xs))
    lapply(xs, function(x) ids %in% x) %>% do.call(what=rbind) %>% t %>%
    limma::vennCounts()
  }
}

