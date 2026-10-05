.PHONY: install deps test build

deps:
	Rscript scripts/install_dependencies.R

install: deps
	R CMD INSTALL .

test:
	Rscript -e 'testthat::test_local(".")'

build:
	R CMD build .
