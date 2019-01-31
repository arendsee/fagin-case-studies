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

do_more <- function(classStr, feats, profiler=profiler){
    # need to generalize this ...
    strOrigins <- list(
        orphan        = classStr[grepl("^[^O]+$", classStr)]
      , s5            = classStr[grepl("O[^O]{4}", classStr, perl=TRUE)]
      , s4            = classStr[grepl(".{1}O[^O]{3}", classStr, perl=TRUE)]
      , s3            = classStr[grepl(".{2}O[^O]{2}", classStr, perl=TRUE)]
      , s2            = classStr[grepl(".{3}O[^O]{1}", classStr, perl=TRUE)]
      , Saccharomyces = classStr[grepl(".{4}O", classStr, perl=TRUE)]
    )
    feats <- lapply(feats, function(d){
        d[, c("seqid",
              "gen", "trn", "orf", # ORFic
              "nuc", "cds", "exo", # non-ORFic
              "rna", "tec", "una", "ind", "nst", "scr" # Unknown
              )]
    })
    strSums <- lapply(strOrigins, function(x) summary(factor(x)))

    group_profile <- function(group, origins, feats, profiler){
        mats <- lapply(names(origins[[group]]), profiler, feats)
        names(mats) <- names(origins[[group]])
        mats
    }

    profiles <- lapply(names(strOrigins)[-1], group_profile, origins=strOrigins, feats=feats, profiler=profiler)
    names(profiles) <- names(strOrigins)[-1]
    list(
        strOrigins = strOrigins
      , strSums    = strSums
      , profiles   = profiles
    )
}

get_yeast_uno_comp <- function(m, strata){

  gene_profiler <- function(seqid, feats){
      y <- lapply(feats, function(f) f[f$seqid == seqid, -1]) %>% do.call(what=rbind)
      rownames(y) <- sub("feature_table/query/Saccharomyces_", "", rownames(y))
      y
  }

  classStr <- rmonad::get_value(m, tag='query_origins')[[1]]$classStr
  feats <- rmonad::get_value(m, tag='feature_table/query')
  x <- do_more(classStr, feats, profiler=gene_profiler)

  .labels <- rmonad::get_value(m, tag="query_labels")[[1]]
  strict_origins <- fagin:::.determine_origins(.labels, con, get_classifier=by_strict)

  y <- do_more(strict_origins$classStr, feats, profiler=gene_profiler)

  list(
    comparisons = list(
       orphan        = merge_named_vectors(x=x$strSums[[1]], y=y$strSums[[1]], .fill=0)
     , s5            = merge_named_vectors(x=x$strSums[[2]], y=y$strSums[[2]], .fill=0)
     , s4            = merge_named_vectors(x=x$strSums[[3]], y=y$strSums[[3]], .fill=0)
     , s3            = merge_named_vectors(x=x$strSums[[4]], y=y$strSums[[4]], .fill=0)
     , s2            = merge_named_vectors(x=x$strSums[[5]], y=y$strSums[[5]], .fill=0)
     , Saccharomyces = merge_named_vectors(x=x$strSums[[6]], y=y$strSums[[6]], .fill=0)
    ),
    venn = make_venn_function(strata, x, y)
  )
}
