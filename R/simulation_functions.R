simulate_pairs_from_empirical <- function(nSim){
  message("Loading Cheung data set")
  cheung <-
    read.table("data/cheung_count_table.txt",
               header=TRUE,
               stringsAsFactors=FALSE,
               row.names=1) %>%
    {.[rowSums(.)>0,]}

  message("Computing empirical mean-dispersion pairs")
  disp <- DGEList(cheung) %>%
    estimateGLMCommonDisp %>%
    estimateGLMTrendedDisp %>%
    estimateGLMTagwiseDisp

  pairs <- tibble(mu=round(rowMeans(disp$counts)),
                  d=disp$trended.dispersion) %>%
    summarise(d = mean(d), .by = mu) %>%
    arrange(mu)

  mu_emp <- pairs$mu %>% unname
  d_emp <- pairs$d %>% unname

  message("Simulating pairs from empirical distribution")
  tibble(mu0 = rnbinom(n=nSim,mu=mean(rowMeans(disp$counts)),size=disp$common.dispersion),
         size = sapply(mu0, function(x) d_emp[which.min(abs(mu_emp-min(x,max(mu_emp))))]))

}


# simulate a shape-constrained trend
# set a fix shape using `bss` or otherwise a random shape
# set number of time points using nTP
simulate_mean <- function(nTP, bss = c("micv","cv","mdcx","cx","tp"), cp = T){

  if(length(bss)==1) if(bss=="null") return(rep(0,nTP))

  if(!all(bss %in% c("micv","cv","mdcx","cx","tp"))) stop("bad shape")

  # set patterns
  m0 <- tibble(x2 = c(1,1,1,0,0,0,1),
                 x3 = c(1,1,-1,1,1,0,-1),
                 x4 = c(1,-1,-1,1,-1,1,1)) %>% as.matrix()

  m <-
    rbind(m0,-m0) %>%
    `rownames<-`(c("micv","cv","cv","micvcp1","cvcp1","micvcp2","tp1",
                   "mdcx","cx","cx","mdcxcp1","cxcp1","mdcxcp2","tp2"))

  keep <- rownames(m) %>% map_lgl(~{any(str_detect(.x,bss))})
  m <- m[keep,]

  bs <- sample(bss,1)

  if(bs == "tp"){
    e <- sample(1:2,1)
  } else if(cp){
    if(bs %in% c("micv","mdcx")){
      e <- sample(c("","cp1","cp2"),1)
    } else {
      e <- sample(c("","cp1"),1)
    }
  } else e <- ""

  type <- paste0(bs,e); type
  r <- c(0,cumsum(m[type,]*runif(3)))

  if(bs == "tp"){
    fit <- gam(r ~ s(t, k = 4, fx = TRUE), data = data.frame(t = 0:3/3, r = r))
  } else {
    fit <- scam::scam(r ~ s(t, k = 4, bs = bs), data = data.frame(t = 0:3/3, r = r))
  }

  ts = round(0:(nTP-1)/(nTP-1),2)

  rs = predict(fit, newdata = data.frame(t = ts), type = "response")

  if(str_detect(type,"cp1")){
    rs[1:floor(nTP/4)] <- rs[floor(nTP/4)]
  } else if(str_detect(type,"cp2")){
    rs[1:floor(nTP/2)] <- rs[floor(nTP/2)]
  }

  rs %>%
    {. - .[1]} %>%
    {./max(abs(.))}

}

# simulate count data based on Spies et al method
rcount <- function(mu,size,nMean,nRep){
  round(mean(rnbinom(n=min(round(max(mu/4,1)),nMean),mu=mu,size=size)))
}

