init:
	data get https://datahub.io/arendsee/fagin-yeast

# run yeast case study 
run-yeast:
	Rscript yeast-main.R

# run brassicaceae case study 
run-brass:
	Rscript brass-main.R

clean:
	rm -rf brass-output yeast-output
