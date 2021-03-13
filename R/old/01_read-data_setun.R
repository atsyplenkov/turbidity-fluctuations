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

# source("R/00_custom-functions.R")

####################### Setun river #############################
# 1) Water Level
setun_l <- read_csv("data/raw/level/l_setun_20191116-20200130.csv",
                    skip = 1, col_types = "-cd") %>% 
  set_colnames(c("date", "LEVEL")) %>% 
  mutate(date = mdy_hms(date, tz = "Europe/Moscow"),
         date = lubridate::floor_date(date, unit = "minute")
  ) %>% 
  filter(date >= as.POSIXct("2019-11-16 15:00:00"),
         date <= as.POSIXct("2020-01-30 15:00:00")
  ) %>% 
  bind_rows(
    read_csv("data/raw/level/l_setun_20200130-20200319.csv",
             skip = 1, col_types = "-cd") %>% 
      set_colnames(c("date", "LEVEL")) %>% 
      mutate(date = mdy_hms(date, tz = "Europe/Moscow"),
             date = lubridate::floor_date(date, unit = "30 min")
      ) %>% 
      filter(date > as.POSIXct("2020-01-30 15:30:00"),
             date <= as.POSIXct("2020-03-19 15:00:00")
      )
  )

# 2) Read meteo data
moscow <- rp5("data/raw/meteo/Moscow/27605.16.11.2019.19.03.2020.1.0.0.ru.utf8.00000000.csv")

# Interpolate to 30-min interval
moscow_tidy <- moscow %>% 
  as_tibble() %>% 
  dplyr::select(date = time, RRR, TTT, Pmm) %>% 
  openair::timeAverage(., avg.time = "30 min") %>%
  mutate(Pmm = zoo::na.approx(Pmm, rule = 2))

# 3) Convert pressure to height
setun_l %<>%
  left_join(moscow_tidy, by = "date") %>%  
  # mutate(h = signif(LEVEL - Pmm, 3)) %>% 
  mutate(h = LEVEL/9806.65 - Pmm/9806.65)

# 4) Turbidity
setun_t <- read_xlsx("data/raw/turbidity/t_setun_20191116-20200130.xlsx",
                     sheet = 3,
                     range = "A1:B10819",
                     col_names = c("date", "ntu")) %>% 
  mutate(date = force_tz(date, "Europe/Moscow"),
         date = floor_date(date, "10 min")) %>% 
  filter(!is.na(ntu)) %>% 
  filter(date >= as.POSIXct("2019-11-16 15:00:00"),
         date <= as.POSIXct("2020-01-30 15:00:00")
  ) %>% 
  bind_rows(
    read_delim("data/raw/turbidity/setun_19032020.log", delim = "_",
               skip = 27) %>% 
      janitor::clean_names() %>% 
      filter(!is.na(time_h_m_s)) %>% 
      transmute(date = paste(date_y_m_d, time_h_m_s, sep = " "),
                date = as_datetime(date),
                date = force_tz(date, "Europe/Moscow"),
                date = floor_date(date, "minute"),
                ntu = turbidity) %>% 
      slice(-1) %>% # remove wrong value
      filter(date > as.POSIXct("2020-01-30 15:30:00"),
             date <= as.POSIXct("2020-03-19 15:00:00")
      )
  )

# Merge
setun <- setun_t %>% 
  full_join(setun_l, by = "date") %>% 
  mutate(q = zoo::na.approx(h),
         q = runmed(q, k = 5, endrule = "keep"),
         datetime = date) %>% 
  hydro_events(., q = q, datetime,
               window = 48) %>% 
  dplyr::select(he, date, h, ntu, RRR)

setun_tidy <- setun %>% 
  timeAverage(., "20 min") %>% 
  left_join(., setun[, c("date", "he")],
            by = "date") %>% 
  mutate(delta = zoo::rollapply(ntu, 3,
                                FUN = function(x) max(x) - min(x),
                                fill = NA)) %>% 
  group_by(he) %>% 
  mutate(range = max(ntu, na.rm = T) - min(ntu, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(TI = delta/range)

# 5) Hydrologival events database
setun_db <- setun_tidy %>% 
  group_by(he) %>% 
  summarise(
    start = first(date),
    end = last(date),
    length = as.double(signif(difftime(end, start, units = "hours"), 3)), # hours
    ntu_mean = mean(ntu, na.rm = T),
    TI_mean = mean(TI, na.rm = T),
    TI_sd = sd(TI, na.rm = T)
  ) %>%  
  arrange(start)

setun_db %>% 
  summarise(n = n(),
            mean = mean(length),
            sd = sd(length),
            TI = mean(TI_mean),
            TI_sd = sd(TI_mean),
            TI_max = max(TI_mean),
            TI_min = min(TI_mean)
            )

setun_tidy %>% 
  summarise(mean = mean(ntu, na.rm = T),
            range = max(ntu, na.rm = T) - min(ntu, na.rm = T)
            ) %>% 
  as.data.frame()

# PLOTS
# Plot Hydrological events
setun_tidy %>% 
  mutate(h = zoo::na.approx(h)) %>%
  ggplot() +
  geom_vline(data = setun_db,
             aes(xintercept = start),
             size = .8,
             color = "grey80",
             linetype = "dashed") +
  geom_line(aes(x = date, y = ntu / 2000,
                color = "Мутность"),
             size = .5) +
  geom_line(aes(x = date, y = h,
                color = "Уровень"),
            size = .5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 2000,
                                         breaks = seq(0, 2000, 400),
                                         name = "Оптическая мутность, NTU"),
                     expand = c(.001,.001)) +
  scale_x_datetime(breaks = "2 weeks",
                   date_labels = "%d %b",
                   expand = c(.001,.001)) +
  labs(x = "", y = "Уровень, м") +
  ggsci::scale_color_nejm(name = "") +
  theme_clean() -> setun_plot

setun_rain <-
  setun_tidy %>%
  mutate(day = yday(date),
         RRR = ifelse(is.nan(RRR), NA, RRR)) %>% 
  group_by(day) %>% 
  mutate(p = sum(RRR, na.rm = T)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_line(aes(x = date, y = h / 100),
            alpha = 0) +
  geom_vline(data = setun_db,
             aes(xintercept = start),
             size = .8,
             color = "grey80",
             linetype = "dashed") +
  geom_col(aes(x = date, y = p)) +
  scale_y_reverse() +
  scale_x_datetime(breaks = "2 weeks",
                   date_labels = "%d %b",
                   expand = c(.001,.001)) +
  # scale_y_continuous(limits = c(0, 7))+
  theme_clean() +
  labs(x = "", y = "Осадки (Балчуг), мм") + 
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"))

plot_grid(setun_rain, setun_plot +
            theme(plot.margin = unit(c(0, 0.2, 0, 0.2), "cm")),
          rel_heights = c(1, 3), 
          align = "v",
          nrow = 2) %>% 
  ggsave(filename = "figures/setun_mar_plot.png",
         plot = .,
         dpi = 600,
         w = 13, h = 7)

# Explore TI
setun_db %>% 
  filter(length < 300) %>% 
  ggplot(aes(x = length, y = TI_mean)) +
  geom_point() +
  Add_R2(formula = y ~ log(x)) +
  theme_clean()

setun_db %>% 
  mutate(river = "S") %>% 
  ggplot(aes(x = "",y = TI_mean)) +
  geom_boxplot() +
  geom_jitter() +
  theme_clean()
