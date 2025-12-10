#----------------------------------------------------------------------------
# 1.) Plot calibration curves for p-values for null simulations
# 2.) Summarise computation times
# Input: output from sim_calibration_co.R and sim_calibration_cc.R
#----------------------------------------------------------------------------

library(tidyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(pbmcapply)
library(ggpubr)
library(stringr)
library(colorspace)

rm(list = ls(all = TRUE))
theme_set(theme_classic())
source("R/simulation_functions.R")

model_cols <- c('black',RColorBrewer::brewer.pal(8, name = "Dark2"))
names(model_cols) <- c("cpam","factor","impulsede2","masigpro","tradeseq","tdeseq",
                       "nbamseq","pairwise-all","pairwise-min")

# select data type case-only (co) or case-control (cc)
type <- "cc"

# load simulation results
pv12_files <- dir("output") %>%
  {str_detect(.,paste0("pv_calib_",type)) & str_detect(.,"nTP_12") & str_detect(.,fixed("rep"))} %>%
  which() %>% dir("output")[.]

pv12_sims <- pv12_files %>%
  map(~ readRDS(paste0("output/",.x)))

pv6_files <- dir("output") %>%
  {str_detect(.,paste0("pv_calib_",type)) & str_detect(.,"nTP_6") & str_detect(.,fixed("rep"))} %>%
  which() %>% dir("output")[.]

pv6_sims <- pv6_files %>%
  map(~ readRDS(paste0("output/",.x)))


if(type == "cc"){
  remove_models <- c("trendcatcher","tradeseq")
} else {
  remove_models <- c("trendcatcher")
}


# plots
p6 <-
  pv6_sims %>%
  map(~discard_at(.x, remove_models)) %>%
  plot_pval_calibration +
  coord_equal() +
  ylim(1,5.3) +
  xlim(-log(0.05,10),NA) +
  scale_colour_manual(values = model_cols, aesthetics = c("colour", "fill"));p6

p12 <-
  pv12_sims %>%
  map(~discard_at(.x, remove_models)) %>%
  plot_pval_calibration +
  coord_equal() +
  ylim(1,5.3) +
  xlim(-log(0.05,10),NA) +
  scale_colour_manual(values = model_cols, aesthetics = c("colour", "fill")) +
  theme(legend.position = "right", legend.direction = "vertical", legend.box.spacing = unit(15,"mm"),
        legend.key.spacing = unit(5,"mm"),
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(3.0, 'cm'),
        axis.title.y = element_text(colour = "white")) +
  guides(col=guide_legend(ncol=1,title = NULL, byrow =T)); p12

ggarrange(p6 + labs(title = NULL, subtitle = "6 time points"),
          NA,
          p12 + labs(title = NULL, subtitle = "12 time points"),
          NA,
          labels = c("A","","B",""),
          nrow = 1,
          widths = c(7,1,7,1),
          common.legend = T, legend = "right", legend.grob = get_legend(p12))

# output for manuscript
#ggsave(paste0("figures/qq_calibration_",type,"_2025_12_02.pdf"), width = 7, height = 5)

time12cc <- pv12_sims %>% map_dfr(~{attributes(.x)$sim_pars$timing}) %>% mutate(ntp = 12, type = "cc")
time12co <- pv12_sims %>% map_dfr(~{attributes(.x)$sim_pars$timing}) %>% mutate(ntp = 12, type = "co")
time6cc <- pv6_sims %>% map_dfr(~{attributes(.x)$sim_pars$timing}) %>% mutate(ntp = 6, type = "cc")
time6co <- pv6_sims %>% map_dfr(~{attributes(.x)$sim_pars$timing}) %>% mutate(ntp = 6, type = "co")

# mean (se) timings
timings <-
  bind_rows(time12cc, time12co, time6cc, time6co) %>%
  dplyr::select(-trendcatcher) %>%
  mutate(across(!all_of(c("ntp","type")), ~ .x/60)) %>%
  group_by(ntp, type) %>%
  summarise(across(everything(), ~ mean(.x) %>% round(2), .names = "mean_{col}"),
            across(everything(), ~ (sd(.x)/sqrt(n()) %>% round(4)), .names = "se_{col}")) %>%
  pivot_longer(-c(ntp,type), names_to = c(".value","method"), names_sep = "_") %>%
  mutate(time = paste0(mean," (",round(se,2),")")) %>%
  filter(method != mean) %>%
  dplyr::select(-mean, -se) %>%
  pivot_wider(names_from = method, values_from = time) %>%
  relocate(type) %>%
  arrange(desc(type)); timings

#timings %>% readr::write_csv("output/simulation_timings_2025_11_28.csv")
#readr::read_csv("output/simulation_timings_2025_11_28.csv")

# write as latex table
timings %>%
  mutate(ntp = as.character(ntp),
         design = as.character(type) %>% toupper(),
         type = NULL
         #ntp = ifelse(type == "co", paste0(ntp," TP (CO)"), paste0(ntp," TP (CC)"))
         ) %>%
  relocate(design) %>%
  rename(`#TP` = ntp) %>%
  knitr::kable(format = "latex", booktabs = T, linesep = "", escape = F)
