# functions compute p-values for various case-only models

pv_co_factor <- function(count_matrix,nTP,nRep){
  message("fitting factor model with DESeq2")
  ed.deseq2.factor = data.frame(sample = colnames(count_matrix),
                                time = factor(rep(1:nTP, each = nRep)))

  deseq2_factor =
    DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,
                           colData = ed.deseq2.factor,
                           design = ~ time) %>%
    DESeq2::DESeq(test = "LRT",
          full = ~ time,
          reduced = ~ 1) %>%
    DESeq2::results() %>%
    as_tibble(rownames = "target_id") %>%
    select(target_id,pvalue)

}

pv_co_pairwise <- function(count_matrix,nTP,nRep){
  message("fitting pairwise with DESeq2")
  ed.deseq2 <- data.frame(sample = colnames(count_matrix),
                          time = factor(rep(1:nTP, each = nRep)))

  deseq_pairwise <-
    2:nTP %>%
    set_names %>%
    map(function(t2){
      print(paste0(t2,"/",nTP))
      ed = ed.deseq2 %>%
        filter(time %in% c(1, t2)) %>%
        mutate(time = factor(time))

      cts = count_matrix[, ed$sample]

      DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                     colData = ed,
                                     design = ~ time) %>%
        DESeq2::DESeq() %>% DESeq2::results()
    }) %>% imap(~ .x %>% as_tibble(rownames = "target_id") %>%
                  mutate(tp = .y) %>%
                  select(target_id,pvalue,tp)) %>%
    bind_rows()

}


pv_co_cpam <- function(count_matrix,nTP,nRep,num_cores){
  message("fitting cpam")

  ed.cpam <- data.frame(sample = colnames(count_matrix),
                        time = (rep(1:nTP, each = nRep)))

  m.cpam = cpam::prepare_cpam(
    exp_design = ed.cpam,
    count_matrix = count_matrix,
    gene_level = T,
    aggregate_to_gene = F,
    num_cores = num_cores
  )
  m.cpam = cpam::compute_p_values(m.cpam)

  m.cpam$p_table %>% select(target_id, pvalue = p_val_target)

}


pv_co_impulsede2 <- function(count_matrix,nTP,nRep,num_cores){
  message("fitting impulsede2")

  ed.impulsede2 <-
    data.frame(
      Sample = colnames(count_matrix),
      Time = (rep(1:nTP, each = nRep)),
      Condition = "case"
    )

  ImpulseDE2::runImpulseDE2(count_matrix %>% as.matrix,ed.impulsede2, scaNProc = num_cores)$dfImpulseDE2Results %>%
    as_tibble(rownames = "target_id") %>% select(target_id, pvalue = p)

}


pv_co_masigpro <- function(count_matrix,nTP,nRep,nSim,nMean,seed){
  message("fitting masigpro")
  ed.masigpro = data.frame(
    Time = (rep(1:nTP, each = nRep)),
    Group = rep(1, nTP * nRep),
    Replicates = 1:nRep
  ) %>%
    `row.names<-`(colnames(count_matrix)) %>%
    maSigPro::make.design.matrix()

  maSigPro::p.vector(
    count_matrix,
    ed.masigpro,
    Q = 0.05,
    MT.adjust = "BH",
    min.obs = 6,
    counts = T
  )$p.vector %>%
    as_tibble(rownames = "target_id") %>% dplyr::rename(pvalue = p.value)

}

pv_co_tradeseq <- function(count_matrix,nTP,nRep){
  message("fitting tradeseq")

  time <- rep(1:nTP, each = nRep) %>% as.matrix
  weights <- rep(1,ncol(count_matrix)) %>% as.matrix
  rownames(time) <- colnames(count_matrix)

  tradeSeq::fitGAM(count_matrix, pseudotime=time, cellWeights=weights, nknots=5, parallel = F) %>%
    tradeSeq::associationTest() %>%
    as_tibble(rownames = "target_id") %>%
    select(target_id, pvalue)
}


pv_co_trendcatcher <- function(count_matrix,nTP,nRep,num_cores){
  message("fitting trendcatcher")

  # install_github("jaleesr/TrendCatcher", dependencies = TRUE, build_vignettes = FALSE)
  # https://github.com/jaleesr/TrendCatcher
  # https://jaleesr.github.io/TrendCatcher/index.html

  ts <- rep(1:nTP-1, each = nRep)
  reps <- rep(1:nRep,nTP)
  colnames(count_matrix) <- paste0("A_",ts,"_Rep",reps)
  write.csv(count_matrix, paste0("output/cm_trendCatcher_seed.",seed,".csv"))
  library(TrendCatcher)
  TrendCatcher::run_TrendCatcher(
    count.table.path = paste0("output/cm_trendCatcher_seed.", seed, ".csv"),
    baseline.t = 0,
    time.unit = "h",
    min.low.count = 1,
    para.core.n = num_cores,
    dyn.p.thres = 0.05,
    show.verbose = F
  )$master.table %>%
    as_tibble %>%
    dplyr::select(target_id = Gene, pvalue = dyn.p.val) %>%
    dplyr::mutate(pvalue = as.numeric(pvalue))

}


pv_co_tdeseq <- function(count_matrix,nTP,nRep,nSim,nMean,seed,num_cores){

  message("fitting tdeseq")
  # https://github.com/fanyue322/TDEseq
  # https://fanyue322.github.io/TDEseq

  ed = data.frame(sample = colnames(count_matrix),
                  stage = rep(1:nTP-1, each = nRep),
                  group = "A") %>%
    `row.names<-`(colnames(count_matrix))

  TDEseq::CreateTDEseqObject(counts = count_matrix, meta.data=ed) %>%
    TDEseq::tdeseq(object = ., num.core = num_cores) %>%
    TDEseq::GetTDEseqAssayData('tde') %>%
    as_tibble(rownames = "target_id") %>%
    select(target_id, pvalue)

}

pv_co_nbamseq <- function(count_matrix,nTP,nRep, num_cores){
  ed = data.frame(time = (rep(1:nTP, each = nRep)), row.names = colnames(count_matrix))

  count_matrix[(count_matrix %>% rowSums() != 0),] %>%
    #`mode<-`("integer") %>%
    NBAMSeq::NBAMSeqDataSet(colData = ed, design = ~ s(time)) %>%
    NBAMSeq::NBAMSeq(BPPARAM = BiocParallel::MulticoreParam(num_cores), parallel = T) %>%
    NBAMSeq::results("time") %>%
    as_tibble(rownames = "target_id") %>%
    select(target_id,pvalue)

}
