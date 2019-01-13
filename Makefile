ARC = ${PWD}/arendsee/fagin-yeast/archive
INP = ${PWD}/INPUT

init:
	# data get https://datahub.io/arendsee/fagin-yeast
	mkdir ${INP}
	cp -s ${ARC}/phylostratr-strata.tab ${INP}
	cp -s ${ARC}/phylostratr-species.tab ${INP}
	# symlink all fagin required inputs
	mkdir ${INP}/syn
	mkdir ${INP}/gff
	mkdir ${INP}/fna
	cp -s ${ARC}/*.syn ${INP}/syn
	cp -s ${ARC}/control-list.txt ${INP}
	cp -s ${ARC}/orphan-list.txt ${INP}
	cp -s ${ARC}/tree. ${INP}/tree
	# Explanation of the following kludge:
	# Basically, it seems that DataHub will not allow files with the same
	# basename, so I have to add _annotation and _genome to distinguish
	# betweeen the <species>.gff and <species>.fna files. Then I need to
	# convert back when preparing data for fagin. A better way to do this, on
	# the fagin side, would be to allow fagin to accept a list of species name
	# to filenames, rather than generating expected filenames from directory
	# names.
	cp -s ${ARC}/saccharomyces_uvarum_annotation.gff         ${INP}/gff/saccharomyces_uvarum.gff
	cp -s ${ARC}/saccharomyces_paradoxus_annotation.gff      ${INP}/gff/saccharomyces_paradoxus.gff
	cp -s ${ARC}/saccharomyces_mikatae_annotation.gff        ${INP}/gff/saccharomyces_mikatae.gff
	cp -s ${ARC}/saccharomyces_kudriavzevii_annotation.gff   ${INP}/gff/saccharomyces_kudriavzevii.gff
	cp -s ${ARC}/saccharomyces_eubayanus_annotation.gff      ${INP}/gff/saccharomyces_eubayanus.gff
	cp -s ${ARC}/saccharomyces_cerevisiae_annotation.gff     ${INP}/gff/saccharomyces_cerevisiae.gff
	cp -s ${ARC}/saccharomyces_arboricola_annotation.gff     ${INP}/gff/saccharomyces_arboricola.gff
	cp -s ${ARC}/saccharomyces_uvarum_genome.fna             ${INP}/fna/saccharomyces_uvarum.fna
	cp -s ${ARC}/saccharomyces_paradoxus_genome.fna          ${INP}/fna/saccharomyces_paradoxus.fna
	cp -s ${ARC}/saccharomyces_mikatae_genome.fna            ${INP}/fna/saccharomyces_mikatae.fna
	cp -s ${ARC}/saccharomyces_kudriavzevii_genome.fna       ${INP}/fna/saccharomyces_kudriavzevii.fna
	cp -s ${ARC}/saccharomyces_eubayanus_genome.fna          ${INP}/fna/saccharomyces_eubayanus.fna
	cp -s ${ARC}/saccharomyces_cerevisiae_genome.fna         ${INP}/fna/saccharomyces_cerevisiae.fna
	cp -s ${ARC}/saccharomyces_arboricola_genome.fna         ${INP}/fna/saccharomyces_arboricola.fna

run:
	Rscript main.R

clean:
	rm -rf output
