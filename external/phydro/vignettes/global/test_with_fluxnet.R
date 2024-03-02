rm(list=ls())
library(FluxDataKit)
library(tidyverse)

df = fdk_convert_lsm(
  site="CH-Dav",
  path="~/Downloads/flux_data_kit_beta/fluxes/",
  fluxnet_format = T,
)

df2 = fdk_convert_lsm(
  site="CH-Dav",
  path="~/Downloads/flux_data_kit_beta/fluxes/",
  fluxnet_format = F,
)

df3 = fdk_downsample_fluxnet(
  df, 
  site = "CH-Dav")

df4 = df2 %>% mutate(date = as.Date(time, tz="GMT")) %>% 
  group_by(date) %>% summarize_all(.funs = mean)

# df4.max = df2 %>% filter(time)


df2 %>% filter(time > as.Date("2015-01-01")) %>% 
  with(plot(GPP~as.POSIXct(time), type="l"))
df3 %>% 
  with(points(GPP_NT_VUT_REF~as.POSIXct(TIMESTAMP, tz = "GMT"), type="l", col="red"))
df4 %>% 
  with(points(GPP~as.POSIXct(date, tz= "GMT"), type="l", col="blue"))

df2 %>% filter(time > as.Date("2020-02-01") & time < as.Date("2020-08-1")) %>% 
  with(plot(GPP~as.POSIXct(time), type="l"))
df3 %>% 
  with(points(GPP_NT_VUT_REF~as.POSIXct(TIMESTAMP, tz = "GMT"), type="l", col="red"))
df4 %>% 
  with(points(GPP~as.POSIXct(date, tz= "GMT"), type="l", col="blue"))

df2 %>% filter(time > as.Date("2020-06-01") & time < as.Date("2020-08-1")) %>% 
  with(plot(GPP~as.POSIXct(time), type="l"))
df3 %>% 
  with(points(GPP_NT_VUT_REF~as.POSIXct(TIMESTAMP, tz = "GMT"), type="l", col="red"))
df4 %>% 
  with(points(GPP~as.POSIXct(date, tz= "GMT"), type="l", col="blue"))

df2 %>% filter(as.Date(time) >= as.Date("2020-06-01") &
               as.Date(time) <= as.Date("2020-06-03") ) %>% 
  with(plot(GPP~as.POSIXct(time), type="l"))
df %>% filter(as.Date(as.POSIXct(TIMESTAMP_START, tz = "GMT", format="%Y%m%d%H%M")) >= as.Date("2020-06-01") &
                 as.Date(as.POSIXct(TIMESTAMP_START, tz = "GMT", format="%Y%m%d%H%M")) <= as.Date("2020-06-03") ) %>% 
  with(points(GPP_NT_VUT_REF~as.POSIXct(TIMESTAMP_START, tz = "GMT", format="%Y%m%d%H%M"), type="l", col="cyan4"))
df3 %>% 
  with(points(GPP_NT_VUT_REF~as.POSIXct(TIMESTAMP, tz = "GMT"), type="l", col="red"))
df4 %>% 
  with(points(GPP~as.POSIXct(date, tz= "GMT"), type="l", col="blue"))


test.3day = df2 %>% filter(as.Date(time) >= as.Date("2020-06-01") &
                 as.Date(time) <= as.Date("2020-06-03") ) 
  
test.3day.daily = df4 %>% filter(as.Date(time) >= as.Date("2020-06-01") &
                                   as.Date(time) <= as.Date("2020-06-03") ) 

test.3day.daily_max = df2 %>% mutate(date = as.Date(time, tz="GMT")) %>% 
  group_by(date) %>% summarize_all(.funs = max)


library(rphydro)

par_cost = list(alpha=0.1, gamma=1);

par_plant = list(conductivity=3e-17, psi50=-2, b=2);

options = list(gs_method = "GS_IGF", 
               et_method = "ET_DIFFUSION",
               ftemp_vj_method = "FV_kumarathunge19",
               ftemp_rd_method = "FR_heskel16",
               ftemp_br_method = "FB_atkin15",
               scale_alpha = T)



test.3day.daily.phydro = test.3day.daily %>% summarize(
  tc = Tair-273.16, 
  ppfd = (SWdown - SWup)*2, 
  vpd = VPD*100, 
  co2 = CO2air, 
  elv = elevation, 
  fapar = FPAR, 
  vwind = Wind 
)

rphydro_analytical_r = function(tc, ppfd, vpd,
                                co2, elv, fapar,
                                kphio, psi_soil, rdark,
                                vwind, par_plant, par_cost, 
                                options){
  rphydro::rphydro_analytical(tc, ppfd, vpd,
                              co2, elv, fapar,
                              kphio, psi_soil, rdark,
                              vwind, par_plant, par_cost, 
                              options)  
}


rphydro_instantaneous_analytical_r = function(vcmax25, jmax25,
                                tc, ppfd, vpd,
                                co2, elv, fapar,
                                kphio, psi_soil, rdark,
                                vwind, par_plant, par_cost, 
                                options){
  rphydro::rphydro_instantaneous_analytical(vcmax25, jmax25,
                              tc, ppfd, vpd,
                              co2, elv, fapar,
                              kphio, psi_soil, rdark,
                              vwind, par_plant, par_cost, 
                              options)  
}


acc = test.3day.daily.phydro %>% slice(1) %>% 
        purrr::pmap_dfr(.f = rphydro_analytical_r, 
                         kphio = 0.087, 
                         psi_soil = 0, 
                         rdark = 0.015, 
                         par_plant = par_plant, 
                         par_cost = par_cost, 
                         options = options)



options = list(gs_method = "GS_IGF", 
               et_method = "ET_DIFFUSION",
               ftemp_vj_method = "FV_kumarathunge19",
               ftemp_rd_method = "FR_heskel16",
               ftemp_br_method = "FB_atkin15",
               scale_alpha = T)

test.3day.phydro = test.3day %>% summarize(
    tc = Tair-273.16, 
    ppfd = (SWdown - SWup)*2, 
    vpd = VPD*100, 
    co2 = CO2air, 
    elv = elevation, 
    fapar = FPAR, 
    vwind = Wind 
  )
    
    
inst = test.3day.phydro %>%
        purrr::pmap_dfr(.f = rphydro_instantaneous_analytical_r,
                        vcmax25 = acc$vcmax25,
                        jmax25 = acc$jmax25,
                        kphio = 0.087, 
                        psi_soil = 0, 
                        rdark = 0.015, 
                        par_plant = par_plant, 
                        par_cost = par_cost, 
                        options = options)
        
test.3day %>% with(plot(GPP~as.POSIXct(time), type="l"))
inst %>% with(points(a~as.POSIXct(test.3day$time), type="l", col="green4"))       
inst %>% with(points(I(vcmax*0.5)~as.POSIXct(test.3day$time), type="l", col="green3"))
inst %>% with(points(I(dpsi*10)~as.POSIXct(test.3day$time), type="l", col="cyan3"))

test.3day.phydro %>% mutate(time = test.3day$time) %>% melt("time") %>% 
  ggplot(aes(y=value, x=as.POSIXct(time))) + 
  geom_line(col="aquamarine4") + 
  geom_vline(xintercept = 25, col="pink") + 
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")

inst %>% mutate(time = test.3day$time) %>% melt("time") %>% 
  ggplot(aes(y=value, x=as.POSIXct(time))) + 
  geom_line(col="aquamarine4") + 
  geom_vline(xintercept = 25, col="pink") + 
  theme_classic() +
  theme(strip.background = element_rect(color = "white", size = 1))+
  facet_wrap(~variable, scales = "free")



sites <- FluxDataKit::fdk_site_info |>
  filter(
    sitename == "CH-Dav"
  ) |>
  mutate(
    data_path = "~/Downloads/flux_data_kit_beta/fluxes/"
  )


# output LSM formatted data
fluxnet <- fdk_convert_lsm(
  site = "CH-Dav",
  path = "~/Downloads/flux_data_kit_beta/fluxes/",
  fluxnet_format = TRUE,
  # out_path = tempdir()
)

# Downsample data
fdk_downsample_fluxnet(
  fluxnet,
  site = "CH-Dav", # a site name
  out_path = tempdir(),
  overwrite = TRUE
)

rsofun_data <- fdk_format_drivers(
  site_info = FluxDataKit::fdk_site_info |>
    filter(sitename == "CH-Dav"),
  path = paste0(tempdir(),"/"),
  verbose = TRUE
)
