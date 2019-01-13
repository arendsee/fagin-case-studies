get_config <- function(){
  # The default config uses the sample data for yeast
  con <- config()
  con@synder@offsets = c(1L,1L) # offsets for mummer
  con@synder@trans = "p" # percent identity transform mummer 
  con@alignment@dna2dna_maxspace = 1e8L
  con@input@focal_species     = "Saccharomyces_cerevisiae"
  con@input@gff_dir           = file.path("INPUT", "gff")
  con@input@syn_dir           = file.path("INPUT", "syn")
  con@input@fna_dir           = file.path("INPUT", "fna")
  con@input@tree              = file.path("INPUT", "tree")
  con@input@query_gene_list   = file.path("INPUT", "orphan-list.txt")
  con@input@control_gene_list = file.path("INPUT", "control-list.txt")
  validate_config(con)
  con
}

run_fagin <- function(con){ 
  # if the last result is present in the archive, use it
  m_sav <- file.path(con@archive, 'm.Rda')
  if(file.exists(m_sav)){
      m <- readRDS(m_sav)
  # otherwise run the full analysis
  } else {
      m <- run_fagin(con)
      saveRDS(m, m_sav)
  }
  m
}
