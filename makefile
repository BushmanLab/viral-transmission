README.md: README.Rmd figure/heat.png
	R -e 'knitr::knit("README.Rmd")'
	sed '/^---$$/,/^---$$/d' README.md --in-place

figure: 
	mkdir figure

figure/heat.png: figure phylogenetic.R
	Rscript phylogenetic.R

.PHONY: clean

clean: 
	rm -f figure/heat.png README.md
