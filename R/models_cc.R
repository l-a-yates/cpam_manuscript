# functions compute p-values for various case-control models

pv_cc_factor <- function(count_matrix,ed){
  message("fitting time-as-a-factor model")

  ed.factor <- ed %>% mutate(time = factor(time), condition = factor(condition))
  deseq2_factor =
    DESeq2::DESeqDataSetFromMatrix(countData = count_matrix,
                                   colData = ed.factor,
                                   design = ~ time*condition) %>%
    DESeq2::DESeq(test = "LRT",
                  full = ~ time*condition,
                  reduced = ~ time) %>%
    DESeq2::results() %>%
    as_tibble(rownames = "target_id") %>%
    select(target_id,pvalue)

}

pv_cc_pairwise <- function(count_matrix,ed){
  message("fitting pairwise model")

  deseq_pairwise <-
    (1:nTP-1) %>%
    purrr::set_names() %>%
    purrr::map(function(t2){
      print(paste0(t2+1,"/",nTP))
      ed1 = ed %>%
        filter(time == t2) %>%
        mutate(condition = factor(condition))

      DESeq2::DESeqDataSetFromMatrix(countData = count_matrix[, ed1$sample],
                                     colData = ed1,
                                     design = ~ condition) %>%
        DESeq2::DESeq() %>% DESeq2::results()
    }) %>% imap(~ .x %>% as_tibble(rownames = "target_id") %>%
                  mutate(tp = .y) %>%
                  select(target_id,pvalue,tp)) %>%
    bind_rows()

}

pv_cc_cpam <- function(count_matrix,ed, num_cores = 1, intercept_cc = "1"){
  message("fitting cpam")

  cpam::prepare_cpam(
    exp_design = ed,
    count_matrix = count_matrix,
    gene_level = T,
    model_type = "case-control",
    aggregate_to_gene = F,
    num_cores = num_cores,
    intercept_cc = intercept_cc
  ) %>%
    cpam::compute_p_values() %>%
    {.$p_table} %>%
    select(target_id, pvalue = p_val_target)

}

pv_cc_impulsede2 <- function(count_matrix,ed){
  message("fitting impulsede2")

  ed.impulsede2 <-
    ed %>%
    dplyr::rename(Sample = sample,
                  Time = time,
                  Condition = condition) %>%
    mutate(Condition = str_replace(Condition,"treatment","case"))

  ImpulseDE2::runImpulseDE2(count_matrix,ed.impulsede2,boolCaseCtrl = T,scaNProc = num_cores)$dfImpulseDE2Results %>%
    as_tibble(rownames = "target_id") %>% select(target_id, pvalue = p)

}

pv_cc_masigpro <- function(count_matrix,ed,nRep){
  message("fitting masigpro")

  ed.masigpro = data.frame(
    Time = ed$time,
    Replicates = 1:nRep,
    Control = as.numeric(ed$condition=="control"),
    Treatment = as.numeric(ed$condition=="treatment")
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

pv_cc_tradeseq <- function(count_matrix,ed){
  message("fitting tradeseq")

  time <- ed$time %>% as.matrix
  weights <- rep(1,ncol(count_matrix)) %>% as.matrix
  rownames(time) <- colnames(count_matrix)

  tradeSeq::fitGAM(count_matrix, pseudotime=time, cellWeights=weights,
                   conditions = factor(ed$condition), nknots=5, parallel = F) %>%
    tradeSeq::associationTest() %>%
    as_tibble(rownames = "target_id") %>%
    select(target_id, pvalue)
}
