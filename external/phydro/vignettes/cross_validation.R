### Cross Validation
require(stringr)
# species_master = unique(dat$Species) %>% str_subset(pattern = "Mediterranean", negate = T) %>% str_subset("Glycine max2", negate=T)

species_master = (par_data %>% filter(dpsi_cal == 1) %>% filter(Species != "Eucalyptus pauciflora"))$Species


K=5

for (spp in 1:length(species_master)){

res_master = data.frame()
  
species = species_master[spp]
species_short = paste(substr(strsplit(species, " ")[[1]][-(5:8)], 0,2), collapse = ".")
data1=filter(dat, Species==species)

data1 = data1[sample(1:nrow(data1), replace = F), ]
N = floor(nrow(data1)/(K))
sizes = rep(N,K)
diff = nrow(data1) - sum(sizes)
if (diff>0){
  for (i in 1:diff){
    sizes[i] = sizes[i]+1
  }
}
start_ids = seq(1, nrow(data1), N)
start_ids = c(1, cumsum(sizes)+1)
cat("Species = ", species, "\n")
cat("splits = ", start_ids, " | ", nrow(data1), "\n")
cat("sizes = ", diff(start_ids), "\n")


ids_list = list()
for (i in 1:K){
  ids_list[[i]] = start_ids[i]:(start_ids[i+1]-1)
}

for (iter in 1:5){
  ids_test = unlist(ids_list[iter])
  ids_train = unlist(ids_list[-iter])
  
  dat_train = data1[ids_train,]
  dat_test  = data1[ids_test,]
  
  x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
  error_fun(x0, dat_train, plot=F, dpsi_file="", dpsi_calib = F, inst=T)
  
  scale = abs(x0)
  
  # out = optimr::optimr(par = x0/scale,
  #                      fn = error_fun,
  #                      scale = scale,
  #                      data=dat_train,
  #                      dpsi_file = species,
  #                      dpsi_calib = T,
  #                      inst = T,
  #                      control = list(par_scale = 1, maxit=300)
  # )
  out = optimr::optimr(par = x0/scale,
                 fn = error_fun,
                 data=dat_train,
                 scale=scale,
                 dpsi_file = species,
                 dpsi_calib = T,
                 inst=T,
                 control = list(par_scale = 1, maxit=100)
  )
  
  out$par = out$par*scale
  names(out$par) = c("K.scalar", "P50", "b", "alpha" ,"yxp50")
  
  res = data.frame(as.list(out$par))
  res$height = mean(dat_train$Ht_used)
  res$species = species
  res$iter = iter
  res$err_train = error_fun(x0, dat_train, plot=F, dpsi_file=species, dpsi_calib = T, inst=T)
  res$err_test  = error_fun(x0, dat_test,  plot=F, dpsi_file=species, dpsi_calib = T, inst=T)
  
  res_master = rbind(res_master, res)
  
  cat("-------\n\n")
}




### All data
x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
# names(x0) = c("K.scalar", "P50", "b", "alpha" ,"yxp50")
error_fun(x0, dat_train, plot=F, dpsi_file=species, dpsi_calib = T, inst=T)

scale = abs(x0)

out = optimr::optimr(par = x0/abs(x0),
                     fn = error_fun,
                     data=dat_train,
                     scale=scale,
                     dpsi_file = species,
                     dpsi_calib = T,
                     inst=T,
                     control = list(par_scale = 1, maxit=100)
)

# error_fun(out$par*scale, data1, plot=F, dpsi_file="", scale=1)

out$par = out$par*scale
names(out$par) = c("K.scalar", "P50", "b", "alpha" ,"yxp50")
res = data.frame(as.list(out$par))
res$height = mean(dat_train$Ht_used)
res$species = species
res$iter = -1
res$err_train = error_fun(x0, dat_train, plot=F, dpsi_file=species, dpsi_calib = T, inst=T)
res$err_test  = error_fun(x0, dat_test,  plot=F, dpsi_file=species, dpsi_calib = T, inst=T)

res_master = rbind(res_master, res)

write.csv(res_master, file = paste0("cross_validation/cross_validation_", species_short,".csv"))

print(summary(res_master))



# iter=3
# error_fun(as.numeric(res_master[iter,1:5]), data1[unlist(ids_list[-iter]),], plot=F, dpsi_file="", scale=1)
# error_fun(as.numeric(res_master[iter,1:5]), data1[unlist(ids_list[iter]),], plot=F, dpsi_file="", scale=1)

# error_fun(as.numeric(out$par), data1, plot=F, dpsi_file="", scale=1)
# x0 = as.numeric((par_data %>% filter(Species == species) )[2:6])
# error_fun(as.numeric(x0), data1, plot=F, dpsi_file="", scale=1)

}



res_all = NULL #read.csv(paste0("cross_validation/cross_validation_", "Qu.il", ".csv"))


for (species in species_master){
  cat(species, "\n")
  species_short = paste(substr(strsplit(species, " ")[[1]][-(5:8)], 0,2), collapse = ".")
  res = read.csv(paste0("cross_validation/cross_validation_", species_short, ".csv"))
  res_all = rbind(res_all, res)
}


fitted = res_all %>% filter(iter == -1) %>% dplyr::select(-X)
write.csv(fitted, file="cross_validation/fitted_avg.csv")

crossval = res_all %>% filter(iter != -1) %>% dplyr::select(-X)
cv = crossval %>% group_by(species) %>% 
  summarize(train_mean = mean(err_train), train_sd = sd(err_train), train_med = median(err_train),
            test_mean = mean(err_test), test_sd = sd(err_test), test_med = median(err_test),
            K_mean = mean(K.scalar), K_sd = sd(K.scalar),
            P50_mean = mean(P50), P50_sd = sd(P50), P50_med = median(P50),
            b_mean = mean(b), b_sd = sd(b),
            alpha_mean = mean(alpha), alpha_sd = sd(alpha),
            yxp50_mean = mean(yxp50), yxp50_sd = sd(yxp50)
  ) 
write.csv(cv, file="cross_validation/crossval.csv")



write.csv(dat %>% filter(Species %in% species_master) %>% group_by(Species) %>% summarize(count = n()), file="n.csv")
