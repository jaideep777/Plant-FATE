library(tidyverse)

dat <- readRDS("c:/Users/Jaideep/Downloads/rsofun_driver_data_v3.3.rds")

dat_gfguy <- dat |> 
    dplyr::filter(sitename == "GF-Guy") |> 
    dplyr::pull(forcing) |> 
    purrr::pluck(1) |>
    dplyr::mutate(
        year = lubridate::year(date),
        month = lubridate::month(date),
        decimal_year = lubridate::decimal_date(date),
    ) |> 
    dplyr::select(year, month, decimal_year, temp, vpd, par=ppfd) |>
    group_by(year, month) |>
    summarize(
        decimal_year = mean(decimal_year, na.rm = TRUE),
        temp = mean(temp, na.rm = TRUE),
        vpd = mean(vpd, na.rm = TRUE),
        par = mean(par, na.rm = TRUE),
    ) |>
    dplyr::mutate(
        par = par*1e6,    # convert from mol/m2/s to umol/m2/s
        vpd = vpd/100    # convert from Pa to hPa
    ) |>
    dplyr::mutate(
        par_max = par*4,
        swp = 0.05
    ) |>
    ungroup()

dat_gfguy |> 
    dplyr::select(-year, -month) |>
    tidyr::pivot_longer(-decimal_year) |> 
    ggplot(aes(x=decimal_year, y=value)) +
    geom_line() +
    facet_wrap(~name, scales = "free_y")

dat_gfguy |> 
    readr::write_csv("c:/Users/Jaideep/Downloads/gf-guy_drivers_plantfate.csv")


