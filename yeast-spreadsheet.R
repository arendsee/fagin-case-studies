# [x] gene_name
# [ ] locus_id
# [x] standard_ps_name
# [x] standard_ps_level
# [x] target_species
# [x] homology_class
# [ ] homology_class p-value
# [ ] feature_alignment_length
# [ ] feature_seqid
# [ ] feature_start
# [ ] feature_end
# [ ] clade
# [ ] clade_class

make_yeast_spreadsheet <- function(m, strata){

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
    dplyr::summarize(Agen_target = paste0(Agen_target, "(", Agen_pval, ")", collapse=","))

  aa2orf <- get_value(m, tag="aa2orf.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, Aorf_target=target, Aorf_pval=pval.adj) %>%
    dplyr::filter(Aorf_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(Aorf_target = paste0(Aorf_target, "(", Aorf_pval, ")", collapse=","))

  aa2transorf <- get_value(m, tag="aa2transorf.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, Acds_target=target, Acds_pval=pval.adj) %>%
    dplyr::filter(Acds_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(Acds_target = paste0(Acds_target, "(", Acds_pval, ")", collapse=","))

  gene2genome <- get_value(m, tag="gene2genome.pval/query")[[1]] %>%
    dplyr::select(seqid=query, target_species=species, genomic_target=target, genomic_pval=pval.adj) %>%
    dplyr::filter(genomic_pval <= 0.05) %>%
    dplyr::group_by(seqid, target_species) %>%
    dplyr::summarize(genomic_target = paste0(genomic_target, "(", genomic_pval, ")", collapse=","))

  d <- merge(d, aa2aa, by=c('seqid', 'target_species'), all=TRUE)
  d <- merge(d, aa2orf, by=c('seqid', 'target_species'), all=TRUE)
  d <- merge(d, aa2transorf, by=c('seqid', 'target_species'), all=TRUE)
  d <- merge(d, gene2genome, by=c('seqid', 'target_species'), all=TRUE)

  # get coordinates of each match

}

# Another is the xml file that we discussed before several times, for athal and
# yeast. with column info: locus ID (or similar for yeast), gene name,
# phylostratr-assigned PS.  then for EACH clade: specific assignment (i.e.,Agen,
# Atrn, Aorf, Ncda, Nexo.....Uuna,Unst,Utec) and support(P value) for assignment;
# then the length and coordinates of the CDS and (as applicable) matching ORF or
# NT-sequence; then your Fig4-based overall assignment (you can call it
# "potential overall assignment or some-such)
