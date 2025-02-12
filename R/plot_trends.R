#----------------------------------------------------------------------------
# Plot ROC curves for simulated case-control and case-only data
# loads outputs from sim_trends_co.R and sim_trends_cc.R
#----------------------------------------------------------------------------

library(tidyr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(purrr)
library(iCOBRA)
library(ggpubr)

rm(list = ls(all = TRUE))
theme_set(theme_classic())

# short series, LFC = 1
cdata_cc_6_1 <- readRDS("output/trends_cc_lfc.1_nRep.3_nTP6_seed.977552_nSim30000_propNull.0.9.rds")
cdata_co_6_1 <- readRDS("output/trends_co_lfc.1_nRep.3_nTP6_seed.33744_nSim30000_propNull.0.9.rds")

# long series, LFC = 1
cdata_cc_12_1 <- readRDS("output/trends_cc_lfc.1_nRep.3_nTP12_seed.273966_nSim30000_propNull.0.9.rds")
cdata_co_12_1 <- readRDS("output/trends_co_lfc.1_nRep.3_nTP12_seed.749033_nSim30000_propNull.0.9.rds")

# long series, LFC = 2
cdata_cc_12_2 <- readRDS("output/trends_cc_lfc.2_nRep.3_nTP12_seed.710977_nSim30000_propNull.0.9.rds")
cdata_co_12_2 <- readRDS("output/trends_co_lfc.2_nRep.3_nTP12_seed.461524_nSim30000_propNull.0.9.rds")

# select sims
ntp = 12
lfc = 2

# prepare data
cdata_cc <- get(paste0("cdata_cc_",ntp,"_",lfc))
cdata_co <- get(paste0("cdata_co_",ntp,"_",lfc))

a_cc <- attributes(cdata_cc)$sim_pars;a_cc
a_co <- attributes(cdata_co)$sim_pars;a_co

plot_cdata_cc <- calculate_performance(cdata_cc, binary_truth = "de") %>%
  prepare_data_for_plot(facetted = F)

plot_cdata_co <- calculate_performance(cdata_co, binary_truth = "de") %>%
  prepare_data_for_plot(facetted = F)


# assign fixed colours to models
model_cols <- c('black',RColorBrewer::brewer.pal(7, name = "Dark2"))
names(model_cols) <- c("cpam","factor","impulsede2","masigpro","tradeseq","tdeseq","nbamseq","pairwise")

# plot ROC curves
d_co <-
  plot_cdata_co@fdrtpr %>%
  as_tibble() %>%
  filter(!method %in% c("pairwise.1","pairwise.3")) %>%
  mutate(method = str_replace(method,"pairwise.2","pairwise")) %>%
  select(thr,FDR,TPR,method)

fdrplot_co <-
  plot_cdata_co@fdrtprcurve %>%
  as_tibble %>%
  filter(!method %in% c("pairwise.1","pairwise.3")) %>%
  mutate(method = str_replace(method,"pairwise.2","pairwise")) %>%
  select(TPR,FDR,method) %>%
  ggplot(aes(FDR,TPR, col = method)) +
  geom_vline(xintercept = c(0.01,0.05,0.1), linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(data = d_co, size = 4, alpha = 0.5) +
  scale_color_manual(values = model_cols) +
  #geom_point(data = d_co, size = 3, col = "white") +
  xlim(0,0.2) +
  ylim(0.5,1) +
  guides(col=guide_legend(ncol=1,title = NULL, byrow =T))


d_cc <-
  plot_cdata_cc@fdrtpr %>%
  as_tibble() %>%
  select(thr,FDR,TPR,method)

fdrplot_cc <-
  plot_cdata_cc@fdrtprcurve %>%
  as_tibble %>%
  select(TPR,FDR,method) %>%
  ggplot(aes(FDR,TPR, col = method)) +
  geom_vline(xintercept = c(0.01,0.05,0.1), linetype = "dashed") +
  geom_line(linewidth = 0.7) +
  geom_point(data = d_cc, size = 4, alpha = 0.5) +
  #geom_point(data = d_cc, size = 2, col = "white") +
  xlim(0,0.2) +
  ylim(0.5,1) +
  scale_color_manual(values = model_cols)

# combined plot
ggarrange(fdrplot_co + labs(title = NULL, subtitle = "case\u00adonly") +
            ylim(0.5,1), # change to ylim(0.85,1) for lfc = 2
          NA,
          fdrplot_cc + labs(title = NULL, subtitle = "case\u00adcontrol") +
            ylim(0.5,1), # change to ylim(0.85,1) for lfc = 2
          NA,
          labels = c("A","","B",""),
          nrow = 1,
          widths = c(7,1,7,1),
          common.legend = T, legend = "right", legend.grob = get_legend(fdrplot_co))

# export for ms
#ggsave(paste0("figures/fig_fdrtpr_",ntp,"_",lfc,"_2025_01_23.pdf"), width = 7, height = 5)


