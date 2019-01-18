init:
	data get https://datahub.io/arendsee/fagin-yeast

run:
	Rscript main.R

clean:
	rm -rf output
