##########################################################################
#
# Meta-analysis of turbidity fluctuations
# Part 1. Read data
#
# Anatoly Tsyplenkov, Sergey Chalov
# atsyplenkov@gmail.com
#
##########################################################################

library(tidyverse)
library(janitor)
library(openair)
library(magrittr)
library(lubridate)
library(readxl)
library(extrafont)
library(cowplot)
library(atslib)
library(here)
library(signal)
library(zoo)
library(pangaear)
library(data.table)

####################### Djankuat river #############################
# 1) Water discharge
# Download Hydrological data
# https://doi.pangaea.de/10.1594/PANGAEA.895099
djan_q <- pg_data(doi = '10.1594/PANGAEA.895099')[[1]]$data %>% 
  dplyr::select(datetime = 1, q = 2) %>%
  # Convert date/time to a POSIXct
  mutate(datetime = as.POSIXct(strptime(datetime, "%Y-%m-%dT%H:%M"),
                               tz = "Europe/Moscow"),
         datetime = lubridate::force_tz(datetime, "Europe/Moscow")) %>% 
  arrange(datetime)

# 2) Read meteo data
# Download Meteo data
# https://doi.pangaea.de/10.1594/PANGAEA.895696
djan_meteo <- pg_data(doi = '10.1594/PANGAEA.895696')[[1]]$data %>%  
  dplyr::select(date = 1, temp = 2, p = 7) %>% 
  # Convert date/time to a POSIXct
  mutate(date = as_date(date)) %>% 
  arrange(date)

# 3) Turbidity
djan_t <- read_xlsx(here("data", "raw", "djankuat",
                         "turbidity", "00_ANALITE_ALL_DJAN_v2.xlsx")) %>% 
  select(3,4) %>% 
  magrittr::set_colnames(c("datetime", "ntu")) %>% 
  slice(-c(1:14)) %>%  # remove bad entries in the beginning 
  mutate(datetime = force_tz(datetime, "Europe/Moscow"),
         datetime = round_date(datetime, "10 min")) %>% 
  # Remove very high values
  mutate(ntu = case_when(datetime > as.POSIXct("2016-07-11 15:50:00", tz = "Europe/Moscow") &
                           datetime < as.POSIXct("2016-07-15 19:10:00", tz = "Europe/Moscow") ~ NA_real_,
                         TRUE ~ ntu))

# Filter hydrograph data
djan_q %<>% 
  dplyr::filter(datetime >= round_date(first(djan_t$datetime), "1 hour"),
                datetime <= round_date(last(djan_t$datetime), "1 hour"))

# 4) Filter and merge data
djan <- djan_t %>%
  # Split at NA
  # mutate(id = rleid(is.na(ntu)),
  #        id = ifelse(is.na(ntu), NA, id)) %>% 
  # group_by(id) %>% 
  # nest() %>%
  # # Count amount of values in each group
  # mutate(amount = map_int(data, ~nrow(.x))) %>%
  # unnest(c(data)) %>%
  # # Filter itself
  # mutate(amount = ifelse(is.na(id), NA, amount),
  #        ntu_sg = ifelse(is.na(id), NA,
  #                     ifelse(amount > 50, # Filter only big timeseries
  #                            sgolayfilt(ntu, p = 9, n = 15),
  #                            ntu))) %>% 
  # ungroup() %>% 
  dplyr::transmute(datetime, ntu, ntu_sg = ntu) %>% 
  full_join(djan_q, by = "datetime") %>%
  arrange(datetime) %>% 
  mutate(q = zoo::na.approx(q, rule = 2),
         ntu = zoo::na.approx(ntu, rule = 2, maxgap = 2),
         q = runmed(q, k = 5, endrule = "keep")) %>% 
  hydro_events(., q = q, datetime,
               window = 24) 

# 6) Save
save("djan", file = here("data", "tidy", "djan.Rdata"))

# 7) Plot
Sys.setlocale("LC_TIME", "English")

djan_plot <- djan %>% 
  arrange(datetime) %>% 
  mutate(q = zoo::na.approx(q)) %>%
  ggplot() +
  geom_vline(data = djan %>% 
               group_by(he) %>% 
               summarise(start = first(datetime)),
             aes(xintercept = start),
             size = .8,
             color = "grey80",
             linetype = "dashed") +
  geom_line(aes(x = datetime, y = ntu / 1000,
                color = "Turbidity"),
            size = .5) +
  geom_line(aes(x = datetime, y = q,
                color = "Water discharge"),
            size = .5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 1000,
                                         name = "Turbidity, NTU"),
                     expand = c(.001, .001)) +
  scale_x_datetime(breaks = "1 week",
                   date_labels = "%d %b %Y",
                   expand = c(.001, .001)) +
  labs(x = "", y = "Water discharge, m") +
  ggsci::scale_color_nejm(name = "") +
  theme_clean()

djan_rain <-
  djan %>%
  left_join(djan_meteo %>% 
              mutate(datetime = as_datetime(date) + 9*3600),
            by = "datetime") %>% 
  mutate(day = yday(datetime)) %>%
  ggplot() +
  geom_line(aes(x = datetime, y = q / 100),
            alpha = 0) +
  geom_vline(data = djan %>% 
               group_by(he) %>% 
               summarise(start = first(datetime)),
             aes(xintercept = start),
             size = .8,
             color = "grey80",
             linetype = "dashed") +
  geom_rect(data = . %>% group_by(day) %>% 
                 summarise(start = first(datetime),
                           end = last(datetime),
                           p = sum(p, na.rm = T)),
               aes(xmin = start, xmax = end, ymin = 0, ymax = p)) +
  scale_y_reverse() +
  scale_x_datetime(breaks = "2 weeks",
                   date_labels = "%d %b", 
                   expand = c(.001,.001)) +
  theme_clean() +
  labs(x = "", y = "Precipitation, mm") + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0, 0.2, -0.2, 0.2), "cm"))

plot_grid(djan_rain, djan_plot +
            theme(plot.margin = unit(c(-0.2, 0.2, 0, 0.2), "cm")),
          rel_heights = c(1, 3), 
          align = "v",
          nrow = 2) %>% 
  ggsave(filename = here("figures", "supplementary", "suppl2_djan-plot.png"),
         plot = .,
         dpi = 600,
         w = 13, h = 7)
