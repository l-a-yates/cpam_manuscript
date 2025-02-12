# Plot calibration curves for p-values for null simulations
# Loads output from sim_calibration_co.R and sim_calibration_cc.R

library(tidyr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(purrr)
library(DESeq2)
library(pbmcapply)
library(mgcv)
library(ggpubr)
library(stringr)
library(colorspace)

rm(list = ls(all = TRUE))
theme_set(theme_classic())
source("R/simulation_functions.R")

pv_co_6 <- readRDS("output/pv_calib_co_nSim_30000_nTP_6_nRep_3_nMean_500_lfc_1_seed_443324.rds")
pv_co_12 <- readRDS("output/pv_calib_co_nSim_30000_nTP_12_nRep_3_nMean_500_seed_331269.rds")

pv_cc_6 <- readRDS("output/pv_calib_cc_nSim_30000_nTP_6_nRep_3_nMean_500_seed_972192.rds")
pv_cc_12 <- readRDS("output/pv_calib_cc_nSim_30000_nTP_12_nRep_3_nMean_500_seed_625002.rds")


model_cols <- c('black',RColorBrewer::brewer.pal(8, name = "Dark2"))
names(model_cols) <- c("cpam","factor","impulsede2","masigpro","tradeseq","tdeseq",
                       "nbamseq","pairwise-all","pairwise-min")

type <- "co"
pv6 <- paste0("pv_",type,"_6")
pv12 <- paste0("pv_",type,"_12")

p6 <- pv_co_6 %>%
  {.[names(.) !=  "trendcatcher"]} %>%
  map(~ list(pvalues = .x$pvalue,ecdf = ecdf(.x$pvalue))) %>%
  plot_pval_calibration(ci = 0.95) +
  ylim(1,5.5) +
  scale_colour_manual(values = model_cols) +
  xlim(-log(0.05,10),NA);p6

p12 <- pv_co_12 %>%
  {.[!names(.) %in%  c("trendcatcher")]} %>% #"nbamseq"
  map(~ list(pvalues = .x$pvalue,ecdf = ecdf(.x$pvalue))) %>%
  plot_pval_calibration(ci = 0.95) +
  ylim(1,5.5) +
  scale_colour_manual(values = model_cols) +
  xlim(-log(0.05,10),NA) +
  theme(legend.position = "right", legend.direction = "vertical", legend.box.spacing = unit(15,"mm"),
        legend.key.spacing = unit(5,"mm"),
        legend.text = element_text(size = 10),
        legend.spacing.y = unit(3.0, 'cm'),
        axis.title.y = element_text(colour = "white")) +
  guides(col=guide_legend(ncol=1,title = NULL, byrow =T));p12


#leg <- get_legend(p12) %>% as_ggplot()
# case\u00adonly
ggarrange(p6 + scale_colour_manual(values = model_cols) + labs(title = NULL, subtitle = "6 time points"),
          NA,
          p12 + scale_colour_manual(values = model_cols)+ labs(title = NULL, subtitle = "12 time points"),
          NA,
          labels = c("A","","B",""),
          nrow = 1,
          widths = c(7,1,7,1),
          common.legend = T, legend = "right", legend.grob = get_legend(p12))

# output for manuscript
#ggsave(paste0("figures/fig_pvalues_",type,"_2024_11_11.pdf"), width = 7, height = 5)


