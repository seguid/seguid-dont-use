SHELL=bash

BROWSER=google-chrome

all: install check check-cli

requirements:
	Rscript -e "install.packages(c('knitr', 'rmarkdown'))"

install: doc
	Rscript -e "install.packages('.', repos = NULL)"

doc:
	Rscript -e "devtools::document()"

build: doc
	mkdir -p ".local"
	rm .local/seguid_*.tar.gz || true
	cd ".local" && R CMD build ..

check: build
	cd ".local" && R CMD check --as-cran seguid_*.tar.gz

check-cli:
	module load CBI bats-core bats-assert bats-file; \
	(cd tests/; bats *.bats)

coverage-html:
	tf=$$(mktemp --suffix="-report.html"); \
	Rscript -e "c <- covr::package_coverage(quiet = FALSE); print(c); r <- covr::report(c, file='$${tf}'); utils::browseURL(r, browser = '$(BROWSER)')"

incl/OVERVIEW.md: vignettes/seguid-overview.Rmd
	Rscript -e "rmarkdown::render('vignettes/seguid-overview.Rmd', rmarkdown::md_document(), output_dir = '$(@D)', output_file = '$(@F)')"

README.md: incl/README.md.rsp incl/OVERVIEW.md
	Rscript -e "R.rsp::rfile('$<', postprocess=FALSE)"

spelling:
	Rscript -e "spelling::spell_check_package()"
	Rscript -e "spelling::spell_check_files(c('NEWS.md', dir('vignettes', pattern='[.]Rmd$$', full.names=TRUE)), ignore=readLines('inst/WORDLIST', warn=FALSE))"

WIN_BUILDER = win-builder.r-project.org
win-builder-devel: .local/seguid_*.tar.gz
	curl -v -T "$?" ftp://anonymous@$(WIN_BUILDER)/R-devel/

win-builder-release: .local/seguid_*.tar.gz
	curl -v -T "$?" ftp://anonymous@$(WIN_BUILDER)/R-release/

win-builder: win-builder-devel win-builder-release
