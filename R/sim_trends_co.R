# Simulates data with null and DE trends in the case-only setting
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
source("R/models_co.R")
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


sim_mu <-
  simulate_pairs_from_empirical(nSim) %>%
  mutate(lfc=lfc) %>%
  mutate(mu0 = 2e7 * mu0 / sum(mu0)) %>%  # normalise to a fixed total counts
  arrange(mu0 < 20) %>% # prevent targets with very low counts being used for DE samples
  mutate(target_id = paste0("g", 1:n()))

sim_counts <-
  sim_mu %>%
  mutate(bs = c(sample(bss,nDE, replace = T),rep("null",nNull))) %>%
  rowwise() %>%
  mutate(mu = list((2^(lfc*simulate_mean(nTP = nTP, bss = bs))) %>% {mu0*./mean(.)}),
         counts = list(mu %>% rep(each = nRep) %>%
                         map_dbl(rcount, size = size, nMean = nMean) %>%
                         `names<-`(paste0("X",1:(nTP*nRep)))))

count_matrix <- sim_counts$counts %>% bind_rows %>% as.matrix %>%
  `row.names<-`(sim_counts %>% pull(target_id))

# set simulation truth
de.truth <- sim_counts %>% transmute(target_id, de = as.numeric(bs!="null"))
name_suffix <- paste0("nSim_",nSim,"_nTP_",nTP,"_nRep_",nRep,"_nMean_",nMean,"_seed_",seed)

# save simulation data
write.csv(count_matrix, paste0("output/trend_co_matrix_",name_suffix,".csv"))


#---- fit models

pv <- list()
pv[["factor"]] <- pv_co_factor(count_matrix,nTP,nRep,nSim,nMean,seed)
pv[["cpam"]] <- pv_co_cpam(count_matrix,nTP,nRep,num_cores = num_cores)
pv[["impulsede2"]] <- pv_co_impulsede2(count_matrix,nTP,nRep, num_cores = num_cores)
pv[["masigpro"]] <- pv_co_masigpro(count_matrix,nTP,nRep)
pv[["tradeseq"]] <- pv_co_tradeseq(count_matrix,nTP,nRep)
pv[["tdeseq"]] <- pv_co_tdeseq(count_matrix,nTP,nRep, num_cores = num_cores)
pv[["nbamseq"]] <- pv_co_nbamseq(count_matrix,nTP,nRep,num_cores = num_cores)

pval_pairwise <- pv_co_pairwise(count_matrix,nTP,nRep)
pv[["pairwise.3"]] <- pv[["pairwise.2"]] <- pv[["pairwise.1"]] <-
  pval_pairwise %>% group_by(target_id) %>% summarise(pvalue = min(pvalue))

# p-value adjustment strategies for the pairwise models
# pairwise-1 select min p-value at each time point and then adjust this set
# pairwise-2 adjust p-values for each timepoint, then select min for each timepoint
# pairwise-3 adjust all p-values for all timepoints together, then select mean for each timepoint
padj_pairwise_2 <-
  pval_pairwise %>%
  mutate(pvalue = p.adjust(pvalue, method = "BH"), .by = tp) %>%
  summarise(pvalue = min(pvalue), .by = target_id)

padj_pairwise_3 <-
  pval_pairwise %>%
  mutate(pvalue = p.adjust(pvalue, method = "BH")) %>%
  summarise(pvalue = min(pvalue), .by = target_id)

pval <-
  pv %>% imap(~ .x %>% mutate("{.y}":= pvalue, pvalue= NULL)) %>% purrr::reduce(full_join, by = "target_id") %>%
  {data.frame(select(.,-target_id),row.names = .$target_id)}

# deal with negative p-values in nbamseq
pval$nbamseq[pval$nbamseq<0 & !is.na(pval$nbamseq)] <- 3e-20

padj = pval %>% mutate(across(everything(), ~ p.adjust(.x, method = "BH")))
padj[["pairwise.2"]] <- padj_pairwise_2$pvalue
padj[["pairwise.3"]] <- padj_pairwise_3$pvalue

truth <- de.truth %>% {data.frame(select(.,-target_id),row.names = .$target_id)}

cdata <- iCOBRA::COBRAData(pval = pval, padj = padj,truth = truth)
attributes(cdata)$sim_pars <- rstan::nlist(nSim,nRep,nTP,nMean,seed,bss,prop_null,lfc)
cdata_name <- paste0("lfc.",lfc,"_nRep.",nRep,"_nTP",nTP,"_seed.",seed,"_nSim",nSim,"_propNull.",prop_null)
saveRDS(cdata,paste0("output/trends_co_",cdata_name,".rds"))



