R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc'); pkgdown::clean_site(); pkgdown::init_site(); pkgdown::build_home(); pkgdown::build_articles(); pkgdown::build_tutorials(); pkgdown::build_news()"
#doxygen doxygen/Doxyfile
