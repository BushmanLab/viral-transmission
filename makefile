README.md: README.Rmd
	R -e 'knitr::knit("README.Rmd")'
	sed '/^---$$/,/^---$$/d' README.md --in-place

