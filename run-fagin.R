get_brassicaceae_config <- function(){
  con <- fagin::config() 
  con@archive = "brassicaceae-archive"
  con@synder@offsets = c(1L,1L) # offsets for satsuma?
  con@synder@trans = "d" # proportion transform satsuma
  con@alignment@dna2dna_maxspace = 1e8L
  con@input@focal_species = "Arabidopsis_thaliana"
  con@input@gff <- list(
      Arabidopsis_lyrata   = "brassicaceae-input/gff/Arabidopsis_lyrata.gff"
    , Arabidopsis_thaliana = "brassicaceae-input/gff/Arabidopsis_thaliana.gff"
    , Brassica_rapa        = "brassicaceae-input/gff/Brassica_rapa.gff"
    , Capsella_rubella     = "brassicaceae-input/gff/Capsella_rubella.gff"
    , Eutrema_salsugineum  = "brassicaceae-input/gff/Eutrema_salsugineum.gff"
  )
  con@input@fna <- list(
      Arabidopsis_lyrata   = "brassicaceae-input/fna/Arabidopsis_lyrata.fna"
    , Arabidopsis_thaliana = "brassicaceae-input/fna/Arabidopsis_thaliana.fna"
    , Brassica_rapa        = "brassicaceae-input/fna/Brassica_rapa.fna"
    , Capsella_rubella     = "brassicaceae-input/fna/Capsella_rubella.fna"
    , Eutrema_salsugineum  = "brassicaceae-input/fna/Eutrema_salsugineum.fna"
  )
  con@input@syn <- list(
      Arabidopsis_lyrata  = "brassicaceae-input/syn/Arabidopsis_thaliana.vs.Arabidopsis_lyrata.syn"
    , Brassica_rapa       = "brassicaceae-input/syn/Arabidopsis_thaliana.vs.Brassica_rapa.syn"
    , Capsella_rubella    = "brassicaceae-input/syn/Arabidopsis_thaliana.vs.Capsella_rubella.syn"
    , Eutrema_salsugineum = "brassicaceae-input/syn/Arabidopsis_thaliana.vs.Eutrema_salsugineum.syn"
  )
  con@input@tree <- "brassicaceae-input/tree"
  con@input@query_gene_list <- "brassicaceae-input/brassicaceae-specific-models.txt"
  con@input@control_gene_list <- "brassicaceae-input/non-brassicaceae-specific-models.txt"
  fagin::validate_config(con)
  con
}

get_yeast_config <- function(){
  con <- fagin::config() 
  con@archive = "yeast-archive"
  con@synder@offsets = c(1L,1L) # offsets for mummer
  con@synder@trans = "p" # percent identity transform mummer 
  con@alignment@dna2dna_maxspace = 1e8L
  con@input@focal_species = "Saccharomyces_cerevisiae"
  con@input@gff <- list(
    "Saccharomyces_arboricola"   = "arendsee/fagin-yeast/archive/saccharomyces_arboricola_annotation.gff"
  , "Saccharomyces_cerevisiae"   = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae_annotation.gff"
  , "Saccharomyces_eubayanus"    = "arendsee/fagin-yeast/archive/saccharomyces_eubayanus_annotation.gff"
  , "Saccharomyces_kudriavzevii" = "arendsee/fagin-yeast/archive/saccharomyces_kudriavzevii_annotation.gff"
  , "Saccharomyces_mikatae"      = "arendsee/fagin-yeast/archive/saccharomyces_mikatae_annotation.gff"
  , "Saccharomyces_paradoxus"    = "arendsee/fagin-yeast/archive/saccharomyces_paradoxus_annotation.gff"
  , "Saccharomyces_uvarum"       = "arendsee/fagin-yeast/archive/saccharomyces_uvarum_annotation.gff"
  )
  con@input@fna <- list(
    "Saccharomyces_arboricola"   = "arendsee/fagin-yeast/archive/saccharomyces_arboricola_genome.fna"
  , "Saccharomyces_cerevisiae"   = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae_genome.fna"
  , "Saccharomyces_eubayanus"    = "arendsee/fagin-yeast/archive/saccharomyces_eubayanus_genome.fna"
  , "Saccharomyces_kudriavzevii" = "arendsee/fagin-yeast/archive/saccharomyces_kudriavzevii_genome.fna"
  , "Saccharomyces_mikatae"      = "arendsee/fagin-yeast/archive/saccharomyces_mikatae_genome.fna"
  , "Saccharomyces_paradoxus"    = "arendsee/fagin-yeast/archive/saccharomyces_paradoxus_genome.fna"
  , "Saccharomyces_uvarum"       = "arendsee/fagin-yeast/archive/saccharomyces_uvarum_genome.fna"
  )
  con@input@syn <- list(
    "Saccharomyces_arboricola"   = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae.vs.saccharomyces_arboricola.syn"
  , "Saccharomyces_eubayanus"    = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae.vs.saccharomyces_eubayanus.syn"
  , "Saccharomyces_kudriavzevii" = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae.vs.saccharomyces_kudriavzevii.syn"
  , "Saccharomyces_mikatae"      = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae.vs.saccharomyces_mikatae.syn"
  , "Saccharomyces_paradoxus"    = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae.vs.saccharomyces_paradoxus.syn"
  , "Saccharomyces_uvarum"       = "arendsee/fagin-yeast/archive/saccharomyces_cerevisiae.vs.saccharomyces_uvarum.syn"
  )
  con@input@tree <- "arendsee/fagin-yeast/archive/tree."
  con@input@query_gene_list <- "arendsee/fagin-yeast/archive/orphan-list.txt"
  con@input@control_gene_list <- "arendsee/fagin-yeast/archive/control-list.txt"
  fagin::validate_config(con)
  con
}

.run_fagin <- function(con){ 
  # if the last result is present in the archive, use it
  m_sav <- file.path(con@archive, 'm.Rda')
  if(file.exists(m_sav)){
      m <- readRDS(m_sav)
  # otherwise run the full analysis
  } else {
      m <- fagin::run_fagin(con)
      saveRDS(m, m_sav)
  }
  m
}
