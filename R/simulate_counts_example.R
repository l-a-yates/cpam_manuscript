# Simulates the small data set which is included in the cpam package

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
nSim <- 1e3
nRep <- 5
nTP <- 6
nMean <- 500
#seed <- sample(1e6,1) #440578
set.seed(440578)

bss = c("micv","cv","mdcx","cx","tp")
prop_null <- 0.9
nNull <- round(prop_null*nSim)
nDE <- nSim - nNull
lfc <- 1

sim_mu <-
  simulate_pairs_from_empirical(nSim) %>%
  mutate(lfc=lfc) %>%
  mutate(mu0 = 2e7 * mu0 / sum(mu0)) %>%  # normalise to a fixed total counts
  arrange(mu0 < 20) %>% # prevent targets with very low counts being used for DE samples
  mutate(target_id = paste0("g", sprintf("%03d", sample(1:n()))))

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

# plot random DE gene
sample((1-prop_null)*nSim,1) %>%
  {tibble(counts = count_matrix[.,],
          time = rep(1:nTP, each = nRep))} %>%
      ggplot(aes(time,counts)) +
      geom_smooth() +
      geom_point()

# plot random null gene
((1-prop_null)*nSim + sample(prop_null*nSim,1)) %>%
  {tibble(counts = count_matrix[.,],
          time = rep(1:nTP, each = nRep))} %>%
  ggplot(aes(time,counts)) +
  geom_smooth() +
  geom_point()


# sort count matrix by rownames
count_matrix <- count_matrix[order(rownames(count_matrix)),]

# save simulation data
saveRDS(count_matrix, "data/counts_example.rds")


ed <- tibble(sample = colnames(count_matrix),
           time = rep(1:nTP, each = nRep))


# fit cpam
cpo <- cpam::prepare_cpam(exp_design = ed,
                    count_matrix = count_matrix,
                    gene_level = T,
                    num_cores = 4)
cpo <- cpam::compute_p_values(cpo) # 5 seconds
cpo <- cpam::estimate_changepoint(cpo) # 5 seconds
cpo <- cpam::select_shape(cpo) # 5 seconds
saveRDS(cpo,"data/cpo_example.rds")

#visualize(cpo)



