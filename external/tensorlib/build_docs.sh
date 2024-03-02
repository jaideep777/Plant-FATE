R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); pkgdown::clean_site(); pkgdown::build_site()"
doxygen doxygen/Doxyfile
