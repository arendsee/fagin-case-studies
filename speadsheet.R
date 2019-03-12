make_spreadsheet <- function(m, strata, out){

  seqid <- get_value(m, tag='query_genes')[[1]]

  d <- subset(strata, seqid %in% seqid) %>%
    dplyr::select(seqid, std_strata_name = mrca, std_strata_level = ps)

  vs <- get_value(m, tag='query_labels')[[1]]

  name_conversion <-  c(O1="A_gen", O2="A_trn", O3="A_orf",
                        N1="N_cds", N2="N_exo", N3="N_rna", N4="N_dna",
                        U2="U_ind", U5="U_scr", U6="U_unk", U1="U_una", U3="U_nst", U7="U_tec")

  labels <- rbind_with_name(vs$labels, "target_species") %>%
    dplyr::select(seqid, homology_class=secondary, target_species) %>%
    dplyr::mutate(homology_class = name_conversion[homology_class])

  d <- merge(d, labels, by='seqid')

  clean_name <- function(x, suffix){
    x$target_species <- sub(suffix, "", x$target_species)
    x
  }

  feats <- get_value(m, tag='feature_table') %>%
    rbind_with_name("target_species") %>%
    clean_name("feature_table/query/")

  d <- merge(d, feats, by=c("seqid", "target_species"))

  ori <- get_value(m, tag="query_origins")[[1]]$backbone
  ori <- as.matrix(ori)
  ori[ori == "O"] <- "A"
  ori <- as.data.frame(ori)
  ori$seqid <- rownames(ori)

  d <- merge(d, ori, by='seqid')

  aa2aa <- get_value(m, tag="aa2aa.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, Agen_target=target, Agen_pval=pval.adj) %>%
    dplyr::filter(Agen_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(Agen_target = paste0(Agen_target, "(", Agen_pval, ")", collapse=","), Agen_pval_min=min(Agen_pval))

  aa2orf <- get_value(m, tag="aa2orf.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, Aorf_target=target, Aorf_pval=pval.adj) %>%
    dplyr::filter(Aorf_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(Aorf_target = paste0(Aorf_target, "(", Aorf_pval, ")", collapse=","), Aorf_pval_min=min(Aorf_pval))

  aa2transorf <- get_value(m, tag="aa2transorf.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, Acds_target=target, Acds_pval=pval.adj) %>%
    dplyr::filter(Acds_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(Acds_target = paste0(Acds_target, "(", Acds_pval, ")", collapse=","), Acds_pval_min=min(Acds_pval))

  gene2genome <- get_value(m, tag="gene2genome.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, genomic_target=target, genomic_pval=pval.adj) %>%
    dplyr::filter(genomic_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(genomic_target = paste0(genomic_target, "(", genomic_pval, ")", collapse=","), genomic_pval_min=min(genomic_pval))

  d <- merge(d, aa2aa, by=c('seqid', 'target_species'), all=TRUE)
  d <- merge(d, aa2orf, by=c('seqid', 'target_species'), all=TRUE)
  d <- merge(d, aa2transorf, by=c('seqid', 'target_species'), all=TRUE)
  d <- merge(d, gene2genome, by=c('seqid', 'target_species'), all=TRUE)

  get_alignment_table <- function(xs, prefix=""){
    data.frame(
      seqid = xs@metadata$query,
      target = xs@metadata$target,
      qstart = xs@pattern@range@start,
      qwidth = xs@pattern@range@width,
      qseq = as.character(xs@pattern),
      tstart = xs@subject@range@start,
      twidth = xs@subject@range@width,
      tseq = as.character(xs@subject),
      score = xs@score
    )
  }

  get_filtered <- function(dat, tag){
    merge(
      dat,
      get_value(m, tag=sprintf("%s.pval/query", tag))[[1]] %>%
        dplyr::filter(pval.adj <= 0.05) %>%
        dplyr::select(seqid=query, target_species=species, target=target, pval.adj),
      by=c("seqid", "target_species", "target")
    )
  }

  aatab <- lapply(get_value(m, tag="aa2aa/query"), function(x) get_alignment_table(x$aln)) %>%
    rbind_with_name("target_species") %>% clean_name("aa2aa/query/") %>%
    get_filtered("aa2aa")

  orftab <- lapply(get_value(m, tag="aa2orf/query"), function(x) get_alignment_table(x$aln)) %>%
    rbind_with_name("target_species") %>% clean_name("aa2orf/query/") %>%
    dplyr::group_by(seqid, target_species) %>%
    get_filtered("aa2orf")

  transtab <- lapply(get_value(m, tag="aa2transorf/query"), function(x) get_alignment_table(x$aln)) %>%
    rbind_with_name("target_species") %>% clean_name("aa2transorf/query/") %>%
    get_filtered("aa2transorf")

  gentab <- lapply(get_value(m, tag="gene2genome/query"), function(x){
      map <- x$map
      hits <- as.data.frame(x$hits[map$target]) %>%
        dplyr::select(
          seqid  = query,
          target = second.seqnames,
          tstart = second.start,
          twidth = second.width,
          qstart = first.start,
          qwidth = first.width,
          strand = strand
        )
      hits$id <- map$target
      hits
    }) %>%
    rbind_with_name("target_species") %>% clean_name("gene2genome/query/") %>%
    {
      merge(., get_value(m, tag='gene2genome.pval')[[1]] %>%
          dplyr::select(seqid=query, id=target, target_species=species, pval, pval.adj),
          by=c("seqid", "id", "target_species")
      )
    } %>%
    dplyr::filter(pval.adj < 0.05)

  options(java.parameters = "-Xmx1024m")
  out_file <- out("fagin-data.xlsx")
  if(file.exists(out_file)){
    file.remove(out_file)
  }

  wb <- XLConnect::loadWorkbook(out_file, create=TRUE)
  XLConnect::createSheet(wb, "main")
  XLConnect::createSheet(wb, "aa2aa")
  XLConnect::createSheet(wb, "aa2orf")
  XLConnect::createSheet(wb, "aa2transorf")
  XLConnect::createSheet(wb, "gene2genome")
  XLConnect::writeWorksheet(wb, data=d, sheet="main")
  XLConnect::writeWorksheet(wb, data=aatab, sheet="aa2aa")
  XLConnect::writeWorksheet(wb, data=orftab, sheet="aa2orf")
  XLConnect::writeWorksheet(wb, data=transtab, sheet="aa2transorf")
  XLConnect::writeWorksheet(wb, data=gentab, sheet="gene2genome")
  XLConnect::saveWorkbook(wb)
}
