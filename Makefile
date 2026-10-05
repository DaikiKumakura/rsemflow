.PHONY: deps document install test check build example

deps:
	Rscript scripts/install_dependencies.R

document:
	Rscript -e 'devtools::document()'

install: deps
	R CMD INSTALL .

test:
	Rscript -e 'testthat::test_local(".")'

build:
	R CMD build .

check: build
	R CMD check --no-manual rsemflow_*.tar.gz

example:
	bash inst/extdata/example/run_example.sh rsemflow-example-output
