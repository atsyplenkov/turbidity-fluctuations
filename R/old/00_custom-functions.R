######################################################################
#
# Tsyplenkov Anatolii
# atsyplenkov@gmail.com
#
#####################################################################


# custom functions --------------------------------------------------------
# The formula to convert temperature to water density
# Reference to [Chen C.T.&Millero F.J., 1986]:
temp_density <- function(t){
  ro = (0.9998395 + 6.7914 * 10^-5 * t - 9.0894 * 10^-6 * t^2 + 1.017 * 10^-7 * t^3 - 1.2846*10^-9*t^4 + 1.1592*10^-11*t^5 - 5.0125 * 10^-14*t^6)*1000
  return(ro)
}

# The function to read data from rp5.ru -- a weather archive
rp5 <- function(csv_path, timezone = "Europe/Moscow") {
  
  meteo <- read.csv(csv_path,
                    header = T, encoding = "UTF-8",
                    sep = ";", skip = 6,
                    row.names = NULL)
  
  names(meteo) <- c("time", names(meteo)[3:ncol(meteo)])
  meteo <- meteo[c("time", "RRR", "T", "Po")]
  colnames(meteo) <- c("time", "RRR", "TTT", "Po")
  
  is.na(meteo$RRR) <- meteo$RRR == "No precipitation" | 
    meteo$RRR == "Trace of precipitation" |
    meteo$RRR == "Следы осадков" |
    meteo$RRR == "Осадков нет" 
  
  meteo %>% 
    mutate(time = as.POSIXct(strptime(time, "%d.%m.%Y %H:%M")),
           time = lubridate::force_tz(time = time,
                                      tzone = timezone),
           RRR = as.numeric(as.character(RRR)),
           days = lubridate::floor_date(time, "day"),
           months = lubridate::floor_date(time, "month"),
           years = lubridate::year(time),
           #Pmm = Po * 13.6 / 1000
           Pmm = Po * 133.322) %>% 
    arrange(time) -> meteo 
  
  return(meteo)
  
}

# GGPLOT2 THEME -----------------------------------------------------------
theme_clean <- function(base_font_family = "Ubuntu",
                        base_font_size = 12,
                        legend = "bottom") {
  
  ggpubr::theme_pubclean(base_family = base_font_family,
                         base_size = base_font_size) +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend,
      strip.background = element_blank()
    )
}

Add_R2 <- function(method = "lm",
                   formula = "y ~ x") {
  
  ggpmisc::stat_poly_eq(aes(label =  paste(stat(eq.label),
                                           stat(rr.label),
                                           sep = "~~~~")),
                        formula = formula,
                        rr.digits = 2,
                        # coef.digits = 2,
                        parse = TRUE)
  
  
}

# Map theme
# SOURCE: https://timogrossenbacher.ch/2018/03/categorical-spatial-interpolation-with-r/
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(family = "Ubuntu", color = "#22211d"),
      # remove all axes
      axis.line = element_blank(),
      # axis.text.x = element_blank(),
      # axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # add a subtle grid
      panel.grid.major = element_line(color = "#dbdbd6", size = 0.6),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#ffffff", color = NA), 
      plot.margin = unit(c(.5, .5, .2, .5), "cm"),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "#ffffff", color = NA), 
      panel.spacing = unit(c(-.1, 0.2, .2, 0.2), "cm"),
      legend.background = element_rect(fill = "#ffffff", color = NA),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11, hjust = 0, color = "#4e4d47"),
      plot.title = element_text(size = 16, hjust = 0.5, color = "#4e4d47"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "#4e4d47", 
                                   margin = margin(b = -0.1, 
                                                   t = -0.1, 
                                                   l = 2, 
                                                   unit = "cm"), 
                                   debug = F),
      plot.caption = element_text(size = 9, 
                                  hjust = .5, 
                                  margin = margin(t = 0.2, 
                                                  b = 0, 
                                                  unit = "cm"), 
                                  color = "#939184"),
      ...
    )
}

# Hydrological events ----------------------------------------------------------
# Defining Local Minimum from HYSEP
# https://github.com/USGS-R/DVstats/blob/master/R/hysep.R
locmin <- function(x, datetime, window = 1) {
  
  timestep <- as.double(signif(difftime(head(datetime)[5],
                                        head(datetime)[4],
                                        units = "hours"), 4))
  
  N2star <- round(window / timestep)
  N2star <- ifelse(N2star %% 2 == 0, N2star  + 1, N2star)
  Nobs <- length(x)
  Ngrp <- ceiling(Nobs / N2star)
  Nfil <- (N2star - 1L) / 2L
  Mid <- as.integer((N2star) / 2)
  LocMin <- sapply(seq(N2star, Nobs), function(i)
    min(x[seq(i - N2star + 1L, i)]) == x[i - Mid]
  )
  LocMin <- c(rep(FALSE, Nfil), LocMin, rep(FALSE, Nfil))
  return(LocMin)
  
}

# HYDROLOGICAL EVENTS -------------------------------------------------------
hydro_events <- function(dataframe,
                         q = NULL,
                         datetime = NULL,
                         window = 1) {
  
  dataframe %>% 
    mutate(q = zoo::na.approx(q, rule = 2)) %>% 
    mutate(LocMin = locmin(x = q,
                           datetime = datetime,
                           window = window),
           test = ifelse(LocMin == FALSE, NA, 1)) -> dataframe
  
  # Remove multiple local minimums
  dataframe$test <- lapply(seq(1, length(dataframe$test)), function(i)
    dataframe$test[i - 1L] == dataframe$test[i])
  
  dataframe$LocMin[dataframe$test == T] <- FALSE
  
  # Name the hydrological events
  dataframe$he <- NA
  dataframe$he[dataframe$LocMin == T] <- 2:(length(dataframe$he[dataframe$LocMin == T])+1)
  dataframe$he[1] <- 1
  
  # Fill the gaps of he's
  dataframe %>%
    mutate(he = zoo::na.locf(he),
           he = as.factor(he)) %>% 
    dplyr::select(-LocMin, - test) %>% 
    as_tibble() -> dataframe 
  
  return(dataframe)
}

# Function to get the difference
delta_ssc <- function(x){
  max(x) - min(x)
}