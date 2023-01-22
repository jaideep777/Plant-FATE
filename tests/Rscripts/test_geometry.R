# R script to test:
dat = read.delim("~/codes/tmodel_cpp/geometric_growth_2.txt")
dat$leaf_area = dat$crown_area * dat$lai
dat$heartwood_fraction = 1-dat$sapwood_fraction
dat$hv = (pi*0.25*dat$diameter^2*dat$sapwood_fraction)/dat$leaf_area

par(mfrow=c(4,4), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$height~dat$diameter)
plot(dat$crown_area~I(dat$height*dat$diameter))
plot(dat$crown_area~dat$height)
plot(dat$crown_area~dat$diameter)
plot(dat$leaf_area~dat$crown_area)
plot(dat$sapwood_turnover_rate~dat$height)
# plot(sqrt(4*dat$crown_area/pi)~dat$height)
  
# par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
plot(dat$sapwood_fraction_ode~dat$height)
points(dat$sapwood_fraction~dat$height, type="l", col="red")
plot(dat$total_mass~dat$height, log="")
points(dat$sapwood_mass~dat$height, col="brown")
points(dat$sapwood_mass_ode~dat$height, col="yellow", type="l")
points(dat$heartwood_mass~dat$height, col="grey")
points(dat$heartwood_mass_ode~dat$height, col="blue", type="l")
plot(I(dat$sapwood_mass/dat$leaf_area)~dat$height)
plot(dat$hv)
plot(dat$functional_xylem_fraction~dat$height)
plot(dat$lai)

# plot(dat$total_mass, log="y")
# points(dat$sapwood_mass, col="goldenrod4")
# points(dat$heartwood_mass, col="red4")
# points(dat$root_mass, col="grey")
# points(dat$leaf_mass, col="green3")

plot(dat$litter_mass)
# 
# par(mfrow=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
# plot(dat$height~dat$i)
# plot(dat$diameter~dat$i)
# plot(dat$total_mass~dat$i)
# points(dat$total_prod~dat$i, type="l", col="red")
# 
# 
# 
# wapprox = function(x){
#   if (x > -exp(-sqrt(2))){
#   # Series about 0, Eq. 3.1
#     w = ((((125*x-64)*x+36)*x/24-1)*x+1)*x;
#   }
#   else{
#     # Series about branch point, -1/e, Eq. 4.22
#     p = sqrt(2*(exp(1)*x+1));
#     w = p*(1-(1/3)*p*(1-(11/24)*p*(1-(86/165)*p*(1-(769/1376)*p))))-1;
#   }
#   w
# }
# 
# x = seq(0,100, length.out = 1000)
# plot(x=x, y=pracma::lambertWp(x), type="l")
# points(x=x, y=log(x)-log(log(x))+0.3, type="l", col="red")
# points(x=x, y=0.7*x^0.35, type="l", col="blue")
# points(x=x, y=0.4*x, type="l", col="green")
# 




