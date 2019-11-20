README.md: README.Rmd figure/heat.png SupplementaryTable1.csv
	R -e 'knitr::knit("README.Rmd")'
	sed '/^---$$/,/^---$$/d' README.md --in-place

figure: 
	mkdir figure

figure/heat.png: figure phylogenetic.R SupplementaryTable1.csv
	Rscript phylogenetic.R

.PHONY: clean

clean: 
	rm -f figure/heat.png README.md
