# Simulates null data in the case-only setting
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
prop_null <- 1
nNull <- round(prop_null*nSim)
nDE <- nSim - nNull
lfc <- 1

sim_mu <-
  simulate_pairs_from_empirical(nSim) %>%
  mutate(lfc=lfc) %>%
  mutate(mu0 = 2e7 * mu0 / sum(mu0)) %>%  # normalise to a fixed total counts
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

name_suffix <- paste0("nSim_",nSim,"_nTP_",nTP,"_nRep_",nRep,"_nMean_",nMean,"_seed_",seed)


write.csv(count_matrix, paste0("output/null_co_matrix_",name_suffix,".csv"))

#---- fit models
pv <- list()
pv[["factor"]] <- pv_co_factor(count_matrix,nTP,nRep,nSim,nMean,seed)
pv[["cpam"]] <- pv_co_cpam(count_matrix,nTP,nRep,num_cores = num_cores)
pv[["impulsede2"]] <- pv_co_impulsede2(count_matrix,nTP,nRep, num_cores = num_cores)
pv[["masigpro"]] <- pv_co_masigpro(count_matrix,nTP,nRep)
pv[["tradeseq"]] <- pv_co_tradeseq(count_matrix,nTP,nRep)
pv[["trendcatcher"]] <- pv_co_trendcatcher(count_matrix,nTP,nRep,num_cores)
pv[["tdeseq"]] <- pv_co_tdeseq(count_matrix,nTP,nRep, num_cores = num_cores)
pv[["nbamseq"]] <- pv_co_nbamseq(count_matrix,nTP,nRep,num_cores = num_cores)

pval_pairwise <- pv_co_pairwise(count_matrix,nTP,nRep)
pv[["pairwise-min"]] <- pval_pairwise %>% group_by(target_id) %>% summarise(pvalue = min(pvalue))
pv[["pairwise-all"]] <- pval_pairwise %>% select(-tp)

# save
attributes(pv)$sim_pars <- rstan::nlist(nSim,nRep,nTP,nMean,seed,bss,prop_null,lfc)
s
saveRDS(pv,paste0("output/pv_calib_co_",name_suffix,".rds"))


