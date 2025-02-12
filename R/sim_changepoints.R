# Simulate data with known changepoints
# Estimate changepoints using cpam and pairwise methods
# Generate plots to compare the performance of the two methods

library(tidyr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(purrr)
library(pbmcapply)
library(mgcv)
library(ggpubr)
library(stringr)

theme_set(theme_classic())

rm(list = ls(all = TRUE))
source("R/simulation_functions.R")
source("R/models_co.R")
map <- purrr::map
select <- dplyr::select
num_cores <- 6

#---- simulation parameters
nSim <- 3e4
nRep <- 3
nTP <- 6
nMean <- 500
seed <- sample(1e6,1); print(seed); set.seed(seed)

bss = c("micv","cv","mdcx","cx","tp")
prop_null <- 0.9
nNull <- round(prop_null*nSim)
nDE <- nSim - nNull
lfc <- 2 # use lfc values from 0.5 to 2 in 0.25 increments

# fix mean expression pattern of the LFC scale
mu <- c(0,0,1,2,3,4)

sim_mu <-
  simulate_pairs_from_empirical(nSim) %>%
  mutate(lfc=lfc) %>%
  mutate(mu0 = 2e7 * mu0 / sum(mu0)) %>%  # normalise to a fixed total counts
  arrange(mu0 < 20) %>% # prevent targets with very low counts being used for DE samples
  mutate(nTP = nTP,
         bs = c(rep("DE",nDE),rep("null",nNull)),
         target_id = paste0("g",formatC(1:n(), width = 5, flag = "0")),
         mu = if_else(bs == "DE", list(mu),list(rep(0,6))))

count_matrix <-
  pbmcmapply(function(mu,mu0,size) {2^(lfc*mu)} %>%
               {mu0*.} %>%
               rep(each = nRep) %>%
               map_dbl(rcount, size = size, nMean = nMean) %>%
               `names<-`(paste0("X",1:(nTP*nRep))),
             sim_mu$mu, sim_mu$mu0,sim_mu$size, mc.cores = 6, SIMPLIFY = T) %>% t %>%
  `row.names<-`(sim_mu$target_id)


name_suffix <- paste0("nSim_",nSim,"_nTP_",nTP,"_nRep_",nRep,"_nMean_",nMean,"_seed_",seed)

# save simulation data
write.csv(count_matrix, paste0("output/count_matrix_cp_",name_suffix,".csv"))

# design
ed <- tibble(sample = paste0("X",1:(nTP*nRep)),
             time = rep(1:nTP-1, each = nRep),
             rep = rep(1:nRep,nTP))

# fit cpam
m.cpam = cpam::prepare_cpam(
  exp_design = ed,
  count_matrix = count_matrix,
  gene_level = T,
  aggregate_to_gene = F,
  normalize = F,
  num_cores = num_cores
)
m.cpam = cpam::compute_p_values(m.cpam) # 1:13
m.cpam = cpam::estimate_changepoint(m.cpam) # 0:51
m.cpam = cpam::select_shape(m.cpam) # 1:39
cp_cpam <- m.cpam$changepoints %>% select(target_id,cp_min,cp_1se)

# fit pairwise
ed.pairwise <- data.frame(sample = colnames(count_matrix),
                          time = factor(rep(1:nTP-1, each = nRep)))

pairwise <-
    {2:nTP-1} %>%
    set_names %>%
    purrr::map(function(t2){
      print(paste0(t2,"/",nTP))
      ed = ed.pairwise %>%
        filter(time %in% c(0, t2)) %>%
        mutate(time = factor(time))

      cts = count_matrix[, ed$sample]

      DESeqDataSetFromMatrix(countData = cts,
                             colData = ed,
                             design = ~ time) %>%
        DESeq %>% DESeq2::results()
    }) %>% imap(~ .x %>% as_tibble(rownames = "target_id") %>%
                  mutate(tp = .y) %>%
                  select(target_id,pvalue,tp)) %>%
    bind_rows()

# extract pairwise changepoints
cp_pw <-
  pairwise %>%
  group_by(tp) %>%
  mutate(p.adj= p.adjust(pvalue)) %>%
  filter(p.adj < 0.05) %>%
  ungroup %>%
  group_by(target_id) %>%
  summarise(cp_pw = min(as.numeric(tp)))

# store all results for the current lfc value
cp_res <-
  sim_mu %>%
  mutate(cp_true = mu %>% map_dbl(~ {.x == .x[1]} %>% which %>% max %>% {.-1})) %>%
  filter(cp_true != 0 & bs != "null") %>%
  left_join(cp_cpam, by = "target_id") %>%
  left_join(cp_pw, by = "target_id") %>%
  select(starts_with("cp_")) %>%
  mutate(lfc = lfc,
         nRep = nRep)


# save results for the current lfc value
readr::write_csv(cp_res, paste0("output/cps_res_LFC_",lfc,"_nRep_",nRep,".csv"))


stop("rerun models with next lfc value")


#----------- plots

# load results for all lfc values
cp_results <-
  c("0.5","0.75","1","1.25","1.5","1.75","2") %>%
  map_dfr(~ readr::read_csv(paste0("output/cps_res_LFC_",.x,"_nRep_3.csv"))) %>%
  mutate(cp_pw = cp_pw - 1)


res <-
  cp_results %>%
  pivot_longer(cols = all_of(c("cp_min","cp_1se","cp_pw")),
               names_to = "method", values_to = "cp") %>%
  group_by(lfc,method) %>%
  summarise(bias = mean(cp - cp_true, na.rm = T),
            bias_se = sd(cp - cp_true, na.rm = T)/sqrt(sum(!is.na(cp))), #
            p0 = mean(cp_true == cp, na.rm = T),
            p0_se = sd(cp_true == cp, na.rm = T)/sqrt(1),
            p0_se_analytic = sqrt(p0*(1-p0)/sum(!is.na(cp))),
            p1 = mean(abs(cp_true - cp) <= 1, na.rm = T),
            p1_se = sd(abs(cp_true - cp) <= 1, na.rm = T)/sqrt(1),
            p1_se_analytic = sqrt(p1*(1-p1)/sum(!is.na(cp))))



scols <- scale_color_manual(values = c("cp_min" = "#08519C",
                                       "cp_1se" = "#CB181D",
                                       "cp_pw" = "grey20"),
                            labels = c("cp_min" = "cpam (min)",
                                       "cp_1se" = "cpam (1se)",
                                       "cp_pw" = "pairwise"))

plot_bias <-
  res %>%
  filter(lfc > 0.25) %>%
  ggplot(aes(lfc,bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "longdash", col = "grey10") +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin = bias - 1*bias_se, ymax = bias + 1*bias_se), width = 0.05) +
  labs(title = NULL,
       x = latex2exp::TeX("log$_2$-fold change (LFC)"),
       y = latex2exp::TeX("estimated CP $-$ true CP"),
       col = NULL) +
  #theme_minimal() +
  scale_x_continuous(breaks = seq(0.5, 2, 0.25),
                     limits = c(0.5, 2)) +
  theme(legend.position = "bottom") +
  scols

plot_p0 <-
  res %>%
  filter(lfc > 0.25) %>%
  ggplot(aes(lfc,p0, col = method)) +
  geom_hline(yintercept = 1, linetype = "longdash", col = "grey10") +
  geom_line(size = 0.5, alpha = 0.5) +
  geom_point(size = 2) +
  #geom_errorbar(aes(ymin = p0 - 2*p0_se_analytic, ymax = p0 + 2*p0_se_analytic), width = 0.05) +
  labs(y = "proportion correct",
       x = latex2exp::TeX("log$_2$-fold change (LFC)"),
       col = NULL) +
  #theme_minimal() +
  scale_x_continuous(breaks = seq(0.5, 2, 0.25),
                     limits = c(0.5, 2)) +
  theme(legend.position = "bottom") +
  scols


ggarrange(plot_bias, NULL, plot_p0, nrow = 1, common.legend = T, legend = "right",
          labels = c("A","","B"), widths = c(3,0.3,3))

# export figure for manuscript
#ggsave("figures/fig_estimate_cp_2024_11_20.pdf", width = 7, height = 3)



