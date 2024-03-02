
# code from https://traitecoevo.github.io/plant/articles/individuals.html

library(plant)

      
# Individual plant    
### TEST of individual plants ###
env <- FF16_fixed_environment(1)
y = grow_plant_to_time(FF16_Individual(s = FF16_Strategy()), seq(0,49,1), env)
# dat=read.delim("/home/jaideep/codes/plant_fate/tests/ind_plant.txt", header=F)
dat=read.delim("/home/jaideep/codes/pspm_package/demo/plant_model/out_single_plant_cm.txt", header=F)
par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(y$state[,"height"]~y$time, type="l")
points(dat$V3~dat$V1, col="blue")
plot(y$state[,"mortality"]~y$time, type="l")
points(dat$V6~dat$V1, col="blue")
plot(y$state[,"fecundity"]~y$time, type="l")
points(dat$V7~dat$V1, col="blue")
plot(y$state[,"area_heartwood"]~y$time, type="l")
points(dat$V9~dat$V1, col="blue")

dat=read.delim("/home/jaideep/codes/pspm_package/demo/plant_model/species_0_X.txt", header=F)
plot(dat$V2~dat$V1, type="o")


lp = read.delim("~/codes/pspm_package/light_profile_ind_plant.txt", header=F)
lp = as.matrix(lp[,-ncol(lp)])
x = exp(seq(log(0.1),log(18), length.out = 100))
matplot(y = x, x= t(lp), type="l", lty=1, col= rainbow(50, end = .75))
  
# points(mortality ~ time, yy, type="l", col="brown")
# points(dat$V3~dat$V1, col="yellow3")

# ### TEST of single cohort ###
# dat2=read.delim("/home/jaideep/codes/plant_fate/tests/cohort.txt", header=F)
# points(dat2$V2~dat2$V1, col="cyan")
# plot(dat2$V3~dat2$V1, col="cyan")
# 
# plot(dat2$V7~dat2$V1, col="green4")
#   

# ### TEST of patch with single cohort ###
# dat2=read.delim("/home/jaideep/codes/plant_fate_/tests/patch.txt", header=F)
# plot(height ~ time, yy, type="l", ylim=c(0,18), main="1")
# lines(dat2$V2~dat2$V1, col="cyan4")
# plot(dat2$V3~dat2$V1, col="cyan")
# 
# plot(dat2$V7~dat2$V1, col="green4")

### TEST of patch with multiple cohorts ###
library(plant)
p0 <- scm_base_parameters("FF16")
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p = p0, mutant = FALSE)
p1$control$environment_rescale_usually = F
p1$seed_rain = 1

ptm <- proc.time()
data1 <- run_scm_collect(p1)
proc.time() - ptm

# ptm <- proc.time()
# equilibrium_seed_rain(p1)
# proc.time() - ptm


t <- data1$time
h <- data1$species[[1]]["height", , ]
m <- data1$species[[1]]["mortality", , ]
ld = data1$species[[1]]["log_density",,]
fe = data1$species[[1]]["seeds_survival_weighted",,]
ha = data1$species[[1]]["area_heartwood",,]
hm = data1$species[[1]]["mass_heartwood",,]
le = data1$env
    
  
matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
  

# heights=read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_hts.txt", header=F, col.names = paste0("V", 1:141))
# lds = read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_lds.txt", header=F, col.names = paste0("V", 1:141))
# vs = read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_vs.txt", header=F, col.names = paste0("V", 1:141))
# ms = read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_ms.txt", header=F, col.names = paste0("V", 1:141))

heights1 = read.delim("~/codes/pspm_package/demo/plant_model/species_0_X.txt", header=F, col.names = paste0("V", 1:141))
lds1 =     read.delim("~/codes/pspm_package/demo/plant_model/species_0_u.txt", header=F, col.names = paste0("V", 1:141))
m1 =       read.delim("~/codes/pspm_package/demo/plant_model/species_0_mort.txt", header=F, col.names = paste0("V", 1:141))
vs1 =      read.delim("~/codes/pspm_package/demo/plant_model/species_0_fec.txt", header=F, col.names = paste0("V", 1:141))
ha1 =      read.delim("~/codes/pspm_package/demo/plant_model/species_0_heart.txt", header=F, col.names = paste0("V", 1:141))
hm1 =      read.delim("~/codes/pspm_package/demo/plant_model/species_0_sap.txt", header=F, col.names = paste0("V", 1:141))
# lds1 = log(lds1)
  
# matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rai  n), xlim=c(0,110), ylim=c(0,18))

par(mfrow=c(3,3), mar=c(2,4,1,1), oma=c(1,1,1,1))

# matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,5), xlim=c(0,10),
#         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
# matplot(heights$V1, heights[,-1], lty=1, col=make_transparent("red", 0.25), type="l",
#         las=1, xlab="Time (y  ears)", ylab="Height (m)", add=T)
# 
# matplot(t, (ld), lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,1),
#         las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain), lwd=2)
# matplot(heights$V1, (lds[,-1]), lty=1, col=make_transparent("red", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Log density", add=T)

matplot(t, h, lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,5), xlim=c(0,10),
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, heights1[,-1], lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", add=T)

matplot(t, (ld), lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,1),
        las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, log(lds1[,-1]), lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Log density", add=T)



matplot(t, ha, lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,5), xlim=c(0,10),
        las=1, xlab="Time (years)", ylab="Heartwood area", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, ha1[,-1], lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Heartwood area", add=T)

matplot(t, hm, lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,1),
        las=1, xlab="Time (years)", ylab="Heartwood mass", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, (hm1[,-1]), lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Heartwood mass", add=T)


matplot(t, m, lty=1, col=make_transparent("black", 0.25), type="l", #ylim=c(0,1),
        las=1, xlab="Time (years)", ylab="Mortality", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, (m1[,-1]), lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Mortality", add=T)


matplot(t, fe, lty=1, col=make_transparent("black", 0.25), type="l", log="",
        las=1, xlab="Time (years)", ylab="Viable Seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, (vs1[,-1]), lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Viable Seeds", add=T)

matplot(t, fe, lty=1, col=make_transparent("black", 0.25), type="l", log="y",
        las=1, xlab="Time (years)", ylab="Viable Seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matplot(heights1$V1, (vs1[,-1]), lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), type="l",
        las=1, xlab="Time (years)", ylab="Viable Seeds", add=T)

      
# # Density distribution
# matplot(x=t(h[seq(1,141,3),]), y=t(log(.0001+exp(ld[seq(1,141,3),]))), type="l", log="x", lty=1, col=rainbow(200/3))
#   
# Light profile 
plot(1e20, ylim=c(0,20), xlim=c(0,1))
for (i in 1:length(le)){
  lines(le[[i]][,"height"]~le[[i]][,"canopy_openness"], col=make_transparent("black"))
  # lines(lez[i,]~leco[i,], col=make_transparent("red"))
}
lp = read.delim("~/codes/pspm_package/demo/plant_model/light_profile_ind_plant.txt", header=F)
lp = as.matrix(lp[,-ncol(lp)])
# x = exp(seq(log(0.01),log(18), length.out = 200))
x = seq(0,20,length.out = 200)
matplot(y = x, x= t(lp), type="l", lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), add=T)

  
  

  

# image(z=as.matrix(heights[,-1]), x=heights[,1], xlim=c(0,105.32), zlim=c(0,20))
# image(z=h, x=t, xlim=c(0,105.32), zlim = c(0,20))
# 
# image(z=log(1+as.matrix(vs[,-1])), x=vs[,1], xlim=c(0,105.32), zlim=c(0,log(1+25000)))
# image(z=log(1+fe), x=t, xlim=c(0,105.32), zlim = c(0,log(1+25000)))

# 
# lez=as.matrix(read.delim("/home/jaideep/codes/plant_fate/tests/le_z.txt", header=F)[-51])
# leco=as.matrix(read.delim("/home/jaideep/codes/plant_fate/tests/le_co.txt", header=F)[-51])

# 
# 
# ### TEST INTERPOLATOR ###
# f1 = read.delim("/home/jaideep/codes/plant_fate/tests/interp.txt", sep=" ", header=F)
# f2 = read.delim("/home/jaideep/codes/plant_fate/tests/interp_orig.txt", sep=" ", header=F)
# 
# plot(f1$V2~f1$V1, type="l", col="red")
# points(f1$V3~f1$V1)
# points(f2$V2~f2$V1, pch=20, col="blue")
# 
# shape = 2
# mean_interval = 30
# scale = (gamma(1.0/shape)/shape/mean_interval)^shape
# p0 = shape*scale^(1.0 / shape) / gamma(1.0 / shape)
# 
# 
