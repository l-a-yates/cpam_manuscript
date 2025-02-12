# Simulates data with null and DE trends in the case-control setting
# Fits various models to the data
# Computes p-values and their adjustments
# Saves the results

library(tidyr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(pbmcapply)
library(mgcv)
library(scam)

theme_set(theme_classic())

rm(list = ls(all = TRUE))
source("R/simulation_functions.R")
source("R/models_cc.R")
map <- purrr::map
select <- dplyr::select
num_cores <- 4

#---- simulation parameters
nSim <- 3e4
nRep <- 3
nTP <- 12
nMean <- 500
seed <- sample(1e6,1); print(seed); set.seed(seed)

bss = c("micv","cv","mdcx","cx","tp")
prop_null <- 0.9
nNull <- round(prop_null*nSim)
nDE <- nSim - nNull
lfc <- 2

# randomly sample mean-dispersion pairs and shapes
sim_mu0 <-
  simulate_pairs_from_empirical(nSim) %>%
  mutate(lfc=lfc) %>%
  mutate(mu0 = 2e7*mu0/sum(mu0)) %>% # normalise to a fixed total counts
  arrange(mu0 < 20) %>% # prevent targets with very low counts being used for DE samples
  mutate(target_id = paste0("g",1:n()),
         nTP = nTP,
         bs_control = sample(bss,nSim,replace=T), # sample shapes for control
         bs_treatment = c(sample(bss,nDE, replace = T),rep("null",nNull))) # sample shapes for case (treatment)


# simulate mean counts across the time series
sim_mu <-
  sim_mu0 %>%
  mutate(mean_control = pbmcmapply(simulate_mean, sim_mu0$nTP,sim_mu0$bs_control, mc.cores = 6, SIMPLIFY = F),
         mean_treatment_diff = pbmcmapply(simulate_mean, sim_mu0$nTP,sim_mu0$bs_treatment, mc.cores = 6, SIMPLIFY = F)) %>%
  rowwise() %>%
  mutate(mean_treatment = list(mean_control + lfc*mean_treatment_diff)) %>%
  ungroup

# simulate replicate counts for control
count_matrix_control <-
  pbmcmapply(function(mu,mu0,size) (2^mu) %>% {mu0*.} %>%
               rep(each = nRep) %>%
               map_dbl(rcount, size = size, nMean = nMean) %>%
               `names<-`(paste0("X",1:(nTP*nRep))),
             sim_mu$mean_control, sim_mu$mu0,sim_mu$size, mc.cores = 6, SIMPLIFY = T) %>% t %>%
  `row.names<-`(sim_mu0$target_id)

# simulate replicate counts for treatment
count_matrix_treatment <-
  pbmcmapply(function(mu,mu0,size) (2^mu) %>% {mu0*.} %>%
               rep(each = nRep) %>%
               map_dbl(rcount, size = size, nMean = nMean) %>%
               `names<-`(paste0("X",(1:(nTP*nRep))+nTP*nRep)),
             sim_mu$mean_treatment, sim_mu$mu0,sim_mu$size, mc.cores = 6, SIMPLIFY = T) %>% t %>%
  `row.names<-`(sim_mu0$target_id)

# combine into a single matrix
count_matrix <- cbind(count_matrix_control,count_matrix_treatment)

# plot a random DE target
sample(nSim*(1-prop_null),1) %>%
  {tibble(control = count_matrix_control[.,],
          treatment = count_matrix_treatment[.,],
          time = rep(1:nTP, each = nRep)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time,value, col = name)) +
      geom_smooth() +
      geom_point() +
      labs(subtitle = .)}


# set simulation truth
de.truth <- sim_mu0 %>% transmute(target_id, de = as.numeric(bs_treatment!="null"))
name_suffix <- paste0("nSim_",nSim,"_nTP_",nTP,"_nRep_",nRep,"_nMean_",nMean,"_seed_",seed)

# save simulation data
write.csv(count_matrix, paste0("output/trend_cc_matrix_",name_suffix,".csv"))

# experimental design
ed <- tibble(sample = paste0("X",1:(2*nTP*nRep)),
             time = rep(rep(1:nTP-1, each = nRep),2),
             condition = c(rep("control",nTP*nRep),rep("treatment",nTP*nRep)))
ed


#---- fit models
pv_cc <- list()
pv_cc[["factor"]] <- pv_cc_factor(count_matrix,ed)
pv_cc[["cpam"]] <- pv_cc_cpam(count_matrix,ed,num_cores = num_cores) # 10 mins
pv_cc[["impulsede2"]] <- pv_cc_impulsede2(count_matrix,ed) # 37 mins
pv_cc[["masigpro"]] <- pv_cc_masigpro(count_matrix,ed,nRep = nRep)
pv_cc[["tradeseq"]] <- pv_cc_tradeseq(count_matrix,ed)
pval_pairwise <- pv_cc_pairwise(count_matrix,ed)
pv_cc[["pairwise"]] <- pval_pairwise %>% group_by(target_id) %>% summarise(pvalue = min(pvalue))

pval <-
  pv_cc %>% imap(~ .x %>% mutate("{.y}":= pvalue, pvalue= NULL)) %>%
  purrr::reduce(full_join, by = "target_id") %>%
  {data.frame(select(.,-target_id),row.names = .$target_id)}

# adjust p-values
padj = pval %>% mutate(across(everything(), ~ p.adjust(.x, method = "BH")))
padj_pairwise_2 <-
  pval_pairwise %>% mutate(pvalue = p.adjust(pvalue, method = "BH")) %>%
  group_by(target_id) %>% summarise(pairwise_2 = min(pvalue))
padj[["deseq_pairwise_2"]] <- padj_pairwise_2$pairwise_2


truth <- de.truth %>% {data.frame(select(.,-target_id),row.names = .$target_id)}
cdata <- iCOBRA::COBRAData(pval = pval, padj = padj,truth = truth)
attributes(cdata)$sim_pars <- rstan::nlist(nSim,nRep,nTP,nMean,seed,bss,prop_null,lfc)
cdata_name <- paste0("lfc.",lfc,"_nRep.",nRep,"_nTP",nTP,"_seed.",seed,"_nSim",nSim,"_propNull.",prop_null)

saveRDS(cdata,paste0("output/trends_cc_",cdata_name,".rds"))


