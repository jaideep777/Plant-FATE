pf_read_outputs = function(input_dir, output_dir, expt_dir){
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))

  l = list(
    # seeds1 = read.delim("seeds.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50)),
    # BA1 = read.csv("basal_area.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    co = read.csv("canopy_openness.csv", header=F, col.names = paste0("V", 1:50)),
    lai_v = read.csv("lai_profile.csv", header=F, col.names = paste0("V", 1:27)),
    traits = read.csv("traits.csv"),
    dat_d = readr::read_csv("D_PFATE.csv"),
    # dat$YEAR = decimal_date(as_date(dat$YEAR, format = "%Y-%m-%d %H:%M:%S GMT (doy = %j)"))
    dat2 = read.csv("Y_PFATE.csv"),
    dat3 = read.csv("Y_mean_PFATE.csv"),
    dist = readr::read_csv("size_distributions.csv", col_names = F),
    x = exp(seq(log(0.01), log(10), length.out=100))
  )

  l$dist = l$dist[,-ncol(l$dist)]
  names(l$dist)[1:2] = c("YEAR", "SPP")
  names(l$Zp)[1] = c("YEAR")
  names(l$co)[1] = c("YEAR")
  names(l$lai_v)[1] = c("YEAR")

  l$dat = l$dat_d %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    group_by(YEAR) %>%
    summarize_all(mean)

  n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
  n_year = length(unique(l$dat2$YEAR))

  setwd(wd_back)

  l
}

pf_cat_outputs = function(list1, list2){
  keys <- unique(c(names(list1), names(list2)))
  keys <- keys[keys != "x"]
  l <- lapply(setNames(keys, keys), function(x) {
      bind_rows(list1[[x]], list2[[x]])
    })
  l$x = list1$x
  l
}

# aggregate_annual = function(l){
#   keys <- unique(c(names(l)))
#   keys <- keys[keys != "x"]
#
#   years = unique(l$dat2$YEAR)
#
#   lagg = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% group_by(as.integer(YEAR)) %>% summarize_all(mean)})
#   lagg$x = l$x
#   lagg
# }

pf_subsample_outputs = function(l, interval = 10){
  keys <- unique(c(names(l)))
  keys <- keys[keys != "x"]
  keys <- keys[keys != "dat_d"]

  years = unique(l$dat2$YEAR)

  lsub = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% filter(as.integer(YEAR) %% interval == 0)})
  lsub$x = l$x
  lsub$dat_d = l$dat_d
  lsub
}

pf_slice_time = function(l, ymin, ymax){
  keys <- unique(c(names(l)))
  keys <- keys[keys != "x"]

  years = unique(l$dat2$YEAR)

  lsub = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% filter(YEAR < ymax & YEAR > ymin)})
  lsub$x = l$x
  lsub
}
