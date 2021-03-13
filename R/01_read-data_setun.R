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

####################### Setun river #############################
# 1) Water Level
setun_l <- read_csv(here("data", "raw", "setun",
                         "level", "l_setun_20191116-20200130.csv"),
                    skip = 1, col_types = "-cd") %>% 
  set_colnames(c("date", "LEVEL")) %>% 
  mutate(date = mdy_hms(date, tz = "Europe/Moscow"),
         date = lubridate::floor_date(date, unit = "minute")
  ) %>% 
  dplyr::filter(date >= as.POSIXct("2019-11-16 15:00:00", tz = "Europe/Moscow"),
                date <= as.POSIXct("2020-01-30 15:00:00", tz = "Europe/Moscow")
  ) %>% 
  bind_rows(
    read_csv(here("data", "raw", "setun",
                  "level", "l_setun_20200130-20200319.csv"),
             skip = 1, col_types = "-cd") %>% 
      set_colnames(c("date", "LEVEL")) %>% 
      mutate(date = mdy_hms(date, tz = "Europe/Moscow"),
             date = lubridate::floor_date(date, unit = "30 min")
      ) %>% 
      dplyr::filter(date > as.POSIXct("2020-01-30 15:30:00", tz = "Europe/Moscow"),
                    date <= as.POSIXct("2020-03-19 15:00:00", tz = "Europe/Moscow")
      )
  )

# 2) Read meteo data
moscow <- rp5(here("data", "raw", "setun", "meteo",
                   "27605.16.11.2019.19.03.2020.1.0.0.ru.utf8.00000000.csv")) %>% 
  # Interpolate to 30-min interval
  as_tibble() %>% 
  dplyr::select(date = time, RRR, TTT, Pmm) %>% 
  openair::timeAverage(., avg.time = "30 min") %>%
  mutate(Pmm = zoo::na.approx(Pmm, rule = 2))

# 3) Convert pressure to height
setun_l %<>%
  left_join(moscow, by = "date") %>%  
  # mutate(h = signif(LEVEL - Pmm, 3)) %>% 
  mutate(h = LEVEL/9806.65 - Pmm/9806.65)

# 4) Turbidity
setun_t <- read_xlsx(here("data", "raw", "setun",
                          "turbidity", "t_setun_20191116-20200130.xlsx"),
                     sheet = 3,
                     range = "A1:B10819",
                     col_names = c("date", "ntu")) %>% 
  mutate(date = force_tz(date, "Europe/Moscow"),
         date = floor_date(date, "10 min")) %>% 
  dplyr::filter(!is.na(ntu)) %>% 
  dplyr::filter(date >= as.POSIXct("2019-11-16 15:00:00", tz = "Europe/Moscow"),
                date <= as.POSIXct("2020-01-30 15:00:00", tz = "Europe/Moscow")
  ) %>% 
  bind_rows(
    read_delim(here("data", "raw", "setun",
                    "turbidity", "setun_19032020.log"),
               delim = "_",
               skip = 27) %>% 
      janitor::clean_names() %>% 
      dplyr::filter(!is.na(time_h_m_s)) %>% 
      transmute(date = paste(date_y_m_d, time_h_m_s, sep = " "),
                date = as_datetime(date),
                date = force_tz(date, "Europe/Moscow"),
                date = ceiling_date(date, "10 min"),
                ntu = turbidity) %>% 
      slice(-1) %>% # remove wrong value
      dplyr::filter(date > as.POSIXct("2020-01-30 15:30:00", tz = "Europe/Moscow"),
                    date <= as.POSIXct("2020-03-19 15:00:00", tz = "Europe/Moscow")
      )
  )

# 5) Filter and merge data
setun <- setun_t %>% 
  mutate(ntu_sg = sgolayfilt(ntu, p = 9, n = 15)) %>% 
  full_join(setun_l, by = "date") %>% 
  arrange(date) %>% 
  mutate(q = zoo::na.approx(h),
         q = runmed(q, k = 5, endrule = "keep"),
         datetime = date) %>% 
  hydro_events(., q = q, datetime,
               window = 48) %>% 
  dplyr::select(he, date, h, ntu, ntu_sg, RRR)

# 6) Save
save("setun", file = here("data", "tidy", "setun.Rdata"))

# 7) Plot
Sys.setlocale("LC_TIME", "English")

setun_plot <- setun %>% 
  arrange(date) %>% 
  mutate(h = zoo::na.approx(h)) %>%
  ggplot() +
  geom_vline(data = setun %>% 
               group_by(he) %>% 
               summarise(start = first(date)),
             aes(xintercept = start),
             size = .8,
             color = "grey80",
             linetype = "dashed") +
  geom_line(aes(x = date, y = ntu / 2000,
                color = "Turbidity"),
            size = .5) +
  geom_line(aes(x = date, y = h,
                color = "Water stage"),
            size = .5) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 2000,
                                         breaks = seq(0, 2000, 400),
                                         name = "Turbidity, NTU"),
                     expand = c(.001,.001)) +
  scale_x_datetime(breaks = "2 weeks",
                   date_labels = "%d %b %Y",
                   expand = c(.001,.001)) +
  labs(x = "", y = "Water stage, m") +
  ggsci::scale_color_nejm(name = "") +
  theme_clean()

setun_rain <-
  setun %>%
  mutate(day = yday(date),
         RRR = ifelse(is.nan(RRR), NA, RRR)) %>%
  ggplot() +
  geom_line(aes(x = date, y = h / 100),
            alpha = 0) +
  geom_vline(data = setun %>% 
               group_by(he) %>% 
               summarise(start = first(date)),
             aes(xintercept = start),
             size = .8,
             color = "grey80",
             linetype = "dashed") +
  geom_rect(data = . %>% group_by(day) %>% 
                 summarise(start = first(date),
                           end = last(date),
                           p = sum(RRR, na.rm = T)),
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

plot_grid(setun_rain, setun_plot +
            theme(plot.margin = unit(c(-0.2, 0.2, 0, 0.2), "cm")),
          rel_heights = c(1, 3), 
          align = "v",
          nrow = 2) %>% 
  ggsave(filename = here("figures", "supplementary", "suppl1_setun-plot.png"),
         plot = .,
         dpi = 600,
         w = 13, h = 7)
