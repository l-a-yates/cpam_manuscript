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
library(pbmcapply)
library(stringr)

rm(list = ls(all = TRUE))
select <- dplyr::select
theme_set(theme_classic())

plots <- list()

# select sim settings
ntp = 6  # 6 or 12
lfc = 1 # 1 or 2
type = "cc" # "co" or "cc"

# load simulation results
cdata_files <- dir("output") %>%
  {str_detect(.,paste0("trends_",type)) &
      str_detect(.,paste0("lfc.",lfc)) &
      str_detect(.,paste0("nTP",ntp)) &
      str_detect(.,fixed("rep"))} %>%
  which() %>% dir("output")[.]; cdata_files

cdatas <- cdata_files %>%
  map(~ readRDS(paste0("output/",.x)))


# created aggregated/combined COBRAData object
pval_combined <- seq_along(cdatas) %>% map_dfr(~cdatas[[.x]]@pval %>% {`rownames<-`(.,paste0(.x,rownames(.)))})
padj_combined <- seq_along(cdatas) %>% map_dfr(~cdatas[[.x]]@padj %>% {`rownames<-`(.,paste0(.x,rownames(.)))})
truth_combined <- seq_along(cdatas) %>% map_dfr(~cdatas[[.x]]@truth %>% {`rownames<-`(.,paste0(.x,rownames(.)))})

cdata_combined <- COBRAData(pval = pval_combined,
                            padj = padj_combined,
                            truth = truth_combined)

cdplotdata_combined <-
  calculate_performance(cdata_combined,
                        binary_truth = "de",
                        aspects = c("fdrtpr", "fdrtprcurve")) %>%
  prepare_data_for_plot(facetted = F)

d_combined <-
  cdplotdata_combined@fdrtpr %>%
  as_tibble() %>%
  select(thr, FDR, TPR, method) %>%
  mutate(thresh = factor(thr, labels = c("0.01","0.05","0.1")))

dcurve_combined <-
  cdplotdata_combined@fdrtprcurve %>%
  as_tibble() %>%
  select(TPR, FDR, method)

# prepare plot data and save intermediate results for each simulation
savename <- paste0("output/cdplotdata_",type,"_lfc",lfc,"_ntp",ntp,".rds")
if(file.exists(savename)){
  cdplotdata <- readRDS(savename)
} else {
  cdplotdata <- list()
  for(i in 1:length(cdatas)){
    print(paste0("Processing simulation ",i," of ",length(cdatas)))
    cdplotdata[[i]] <- calculate_performance(cdatas[[i]],
                                             binary_truth = "de",
                                             aspects = c("fdrtpr", "fdrtprcurve")) %>%
      prepare_data_for_plot(facetted = F)
  }
  saveRDS(cdplotdata, savename)
}

d <- cdplotdata %>%
  map( ~ .x@fdrtpr %>%
         as_tibble() %>%
         select(thr, FDR, TPR, method)) %>% bind_rows(.id = "sim")

dcurve <- cdplotdata %>%
  map( ~ .x@fdrtprcurve %>%
         as_tibble() %>%
         select(TPR, FDR, method)) %>% bind_rows(.id = "sim")

# assign fixed colours to models for consistent plotting
model_cols <- c('black',RColorBrewer::brewer.pal(7, name = "Dark2"))
names(model_cols) <- c("cpam","factor","impulsede2","masigpro","tradeseq","tdeseq","nbamseq","pairwise")

# remove failed methods and adjust names according to type
if(type == "co"){
  d <- d %>%
    filter(!method %in% c("pairwise.1","pairwise.3","trendcatcher")) %>%
    mutate(method = str_replace(method,"pairwise.2","pairwise"))
  dcurve <- dcurve %>%
    filter(!method %in% c("pairwise.1","pairwise.3","trendcatcher")) %>%
    mutate(method = str_replace(method,"pairwise.2","pairwise")) %>%
    filter(TPR > 0.5 & FDR < 0.2)
  d_combined <- d_combined %>%
    filter(!method %in% c("pairwise.1","pairwise.3","trendcatcher")) %>%
    mutate(method = str_replace(method,"pairwise.2","pairwise"))
  dcurve_combined <- dcurve_combined %>%
    filter(!method %in% c("pairwise.1","pairwise.3","trendcatcher")) %>%
    mutate(method = str_replace(method,"pairwise.2","pairwise")) %>%
    filter(TPR > 0.5 & FDR < 0.2)
} else{
  d <- d %>%
    filter(!method %in% c("trendcatcher","tradeseq","masigpro"))
  dcurve <- dcurve %>%
    filter(!method %in% c("trendcatcher","tradeseq","masigpro")) %>%
    filter(TPR > 0.5 & FDR < 0.2)
  d_combined <- d_combined %>%
    filter(!method %in% c("trendcatcher","tradeseq","masigpro"))
  dcurve_combined <- dcurve_combined %>%
    filter(!method %in% c("trendcatcher","tradeseq","masigpro")) %>%
    filter(TPR > 0.5 & FDR < 0.2)
}

# set y-axis minimum based on lfc
if(lfc == 1){
  ylim_min <- 0.5
} else {
  ylim_min <- 0.85
}

# plot FDR-TPR curves
plots[[paste0(type,"_ntp",ntp,"_lfc",lfc)]] <-
  dcurve %>%
  ggplot(aes(FDR,TPR)) +
  geom_vline(xintercept = c(0.01,0.05,0.1), linetype = "dashed") +
  geom_line(linewidth = 0.2, alpha = 0.1, aes(group = interaction(method,sim),col = method)) +
  geom_point(aes(shape = thresh,col = method), data = d_combined, size = 4, alpha = 0.6) +
  geom_line(aes(group = method, col = method), data = dcurve_combined, linewidth = 0.7) +
  scale_color_manual(values = model_cols) +
  xlim(0,0.15) +
  ylim(ylim_min,1) +
  guides(col=guide_legend(ncol=1,title = NULL, byrow =T,order = 1)) +
  labs(shape = NULL)

# inspect plot
plots[[paste0(type,"_ntp",ntp,"_lfc",lfc)]]

# manuscript plots
# after running the above code across all combinations of type, ntp, lfc, generate combined plots:

ggarrange(plots$co_ntp12_lfc1 + labs(title = NULL, subtitle = "case\u00adonly"), NA,
          plots$cc_ntp12_lfc1 + labs(title = NULL, subtitle = "case\u00adcontrol"), NA,
          labels = c("A","","B",""),
          nrow = 1,
          widths = c(7,1,7,1),
          common.legend = T, legend = "right")

#ggsave(paste0("figures/fig_fdrtpr_ntp12_lfc1_2025_12_04_light.pdf"), width = 7, height = 5)

ggarrange(plots$co_ntp6_lfc1, NA,
          plots$cc_ntp6_lfc1, NA,
          labels = c("A","","B",""),
          nrow = 1,
          widths = c(7,1,7,1),
          common.legend = T, legend = "right")

#ggsave(paste0("figures/fig_fdrtpr_ntp6_lfc1_2025_12_04_light.pdf"), width = 7, height = 5)

ggarrange(plots$co_ntp12_lfc2 + labs(title = NULL, subtitle = "case\u00adonly"), NA,
          plots$cc_ntp12_lfc2 + labs(title = NULL, subtitle = "case\u00adonly"), NA,
          labels = c("A","","B",""),
          nrow = 1,
          widths = c(7,1,7,1),
          common.legend = T, legend = "right")

#ggsave(paste0("figures/fig_fdrtpr_ntp12_lfc2_2025_12_04_light.pdf"), width = 7, height = 5)

