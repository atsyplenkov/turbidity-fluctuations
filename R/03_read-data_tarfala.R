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
library(xts)
library(dygraphs)

####################### Tarfala river #############################
# 1) Water discharge
tarfala_q <- readr::read_csv(here("data", "raw", "tarfala",
                                  "discharge", "Tarfalajokk_Q_hr_avg_v2.csv")) %>% 
  dplyr::select(datetime = 1, q = 4) %>% 
  mutate(datetime = force_tz(datetime, "GMT"),
         q = replace(q, is.nan(q), 0))

# 2) Turbidity
tarfala_tarfalabron <- read.delim(here("data","raw","tarfala",
                             "turbidity", "tarfalabron",
                             "tarfala_25_08_17.log"),
                      sep = "\t") %>% as_tibble() %>% 
  select(1:3, 5) %>%
  magrittr::set_colnames(c("date", "time", "ntu", "temp")) %>% 
  mutate(datetime = glue::glue("{as.character(date)} {as.character(time)}"),
         datetime = lubridate::dmy_hms(datetime, tz = "GMT"),
         datetime = round_date(datetime, "10 min")) %>% 
  select(datetime, ntu) %>% 
  dplyr::filter(datetime >= as.POSIXct("2017-08-11 19:50:00", tz = "GMT"),
                datetime < as.POSIXct("2017-08-25 11:00:00", tz = "GMT"))

tarfala_rannan <- read_xls(here("data","raw","tarfala",
                                "turbidity", "rannan",
                                "079851_20180624_2311.xls"),
                           skip = 4) %>% 
  transmute(datetime = parse_date_time(Timestamp,
                                       orders = c("d/m/Y H:M:S"),
                                       tz = "GMT"),
            ntu = Turbidity) %>% 
  slice(-1:-8) %>% 
  dplyr::filter(datetime < as.POSIXct("2018-06-24 13:00:00", tz = "GMT"),
                datetime > as.POSIXct("2017-08-25 18:00:00", tz = "GMT")) %>% 
  mutate(ntu = ifelse(ntu < 1, NA_real_, ntu),
         ntu = na.approx(ntu)) 

tarfala_rannan %>% 
  mutate(ntu_sg = sgolayfilt(ntu, p = 9, n = 11)) %>%
  as.xts(x = ., order.by = .$datetime) %>% 
  dygraph()

# 4) Filter and merge data
rannan <- tarfala_rannan %>% 
  mutate(ntu_sg = sgolayfilt(ntu, p = 9, n = 11)) %>%
  left_join(tarfala_q, by = "datetime") %>% 
  arrange(datetime) %>% 
  mutate(q = zoo::na.approx(q, rule = 2),
         q = runmed(q, k = 5, endrule = "keep")) %>% 
  # Split at NA
  mutate(q = replace(q, q == 0, NA)) %>% 
  mutate(id = is.na(q),
         id = rleid(id),
         id = replace(id, is.na(q), NA_integer_)) %>% 
  group_by(id) %>% 
  nest() %>% 
  mutate(data = ifelse(!is.na(id),
                       map(data, ~hydro_events(.x, q, datetime, 24)),
                       data)) %>% 
  unnest(data) %>% 
  ungroup() %>% 
  mutate(he = paste(id, he, sep = "_"),
         he = replace(he, is.na(id), NA_character_),
         ntu = replace(ntu, is.na(id), NA_real_),
         ntu_sg = replace(ntu_sg, is.na(id), NA_real_),) %>% 
  select(he, datetime:q)

bron <- tarfala_tarfalabron %>% 
  mutate(ntu_sg = sgolayfilt(ntu, p = 8, n = 21)) %>% 
  mutate(ntu = ifelse(datetime > as.POSIXct("2017-08-13 9:10:00") &
                        datetime < as.POSIXct("2017-08-13 17:40:00"),
                      NA_real_,
                      ntu),
         ntu_sg = ifelse(datetime > as.POSIXct("2017-08-13 9:10:00") &
                        datetime < as.POSIXct("2017-08-13 17:40:00"),
                      NA_real_,
                      ntu_sg),)

bron %>%
  xts(x = .[, -1], order.by = .$datetime, tzone = "GMT") %>%
  dygraph()



# db_he <- rannan %>% 
#   group_by(he) %>% 
#   summarise(start = first(datetime),
#             end = last(datetime)) %>% 
#   drop_na()
# 
# p <- rannan %>%
#   # filter(he %in% unique(db_he$he)) %>%
#   xts(x = .[, c("q")], order.by = .$datetime) %>%
#   dygraph()
# 
# for (i in 1:nrow(db_he)) {
#   p <- p %>% dyShading(from = db_he$start[i],
#                        to = db_he$end[i],
#                        color = "#FFE6E6")
# }
# 
# print(p)


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
save("setun", file = here("data", "tidy", "setun.Rdata"))

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
