TMT_Contamination <- c(126.1277, 127.1247, 127.1310, 128.1281, 128.1344, 129.1314, 129.1377, 130.1348, 130.1411, 131.1381, 131.1444, 
                       132.1415, 132.1478, 133.1448, 133.1511, 134.1482)





fitG <- function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }
rSquared_tmt <- function(tmt_fit, y, x){
  
  ssTot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  ssRes <- sum((y - tmt_fit[3]*dnorm(x,tmt_fit[1],tmt_fit[2]))^2, na.rm = TRUE)
  r2 <- 1 - ssRes/ssTot
  return(r2)
}






analyze_Immonium_TMT <- function (file, width = 0.001, ions = immoniumIons, fixSigma = T){ 

  message(sprintf("Reading file [%s]", file))
  msrun <- openMSfile(file, backend = "Ramp")
  hd <- header(msrun)
  immscans <- which(hd$collisionEnergy > 45 & (hd$msLevel > 
                                                 1))
  message(sprintf("\t%d immonium scans found", length(immscans)))
  result <- data.frame()
  message("Processing:")
  result <- hd %>% filter(seqNum %in% immscans) %>% group_by(seqNum) %>% 
    do({
      ires <- data.frame()
      dd = .
      ii <- dd$seqNum[[1]]
      ss <- peaks(msrun, ii)
      ss_range <- diff(range(ss[, 1]))
      for (ion_ in names(ions)) {
        mz = monoMass(ions[[ion_]])
        peaks <- c("13C")
        if ("N" %in% names(ions[[ion_]]) & ions[[ion_]]["N"] >
            0)
          peaks <- c(peaks, "15N")
        if ("H" %in% names(ions[[ion_]]) & ions[[ion_]]["H"] >
            0)
          peaks <- c(peaks, "2H")
        if ("O" %in% names(ions[[ion_]]) & ions[[ion_]]["O"] >
            0)
          peaks <- c(peaks, "18O")
        if (ss_range > 5) {
          cres <- get_isopeaks(ss, mz, width = width,
                               npoint = 6, fixSigma = fixSigma, peaks = peaks)
        }
        else {
          cres <- get_isopeaks_nomono(ss, mz, width = width,
                                      npoint = 6, fixSigma = fixSigma)
        }
        if (nrow(cres) >= 1)
          ires <- ires %>% bind_rows(cres %>% mutate(rt = dd$retentionTime[[1]],
                                                     tic = dd$totIonCurrent[[1]], bpI = dd$basePeakIntensity[[1]],
                                                     bpMZ = dd$basePeakMZ[[1]], ion = ion_, PrecursorMZ = dd$precursorMZ[[1]]))
      }
      tmt_count = 0
      max_tmt = 0
      sum_tmt = 0
      tmt_RSq = c()
      for (j_tmt in TMT_Contamination){
        p_int <- findInterval(c(j_tmt - 0.008, j_tmt + 0.008), ss[, 1])
        # if (diff(p_int) < npoint) {
        #   warning(sprintf("No peak found at mz=%.4f with tol=%.4f", 
        #                   j_tmt, tol))
        #   return(data.frame())
        # }
        sub_ss <- ss[p_int[1]:p_int[2], ]
        istart = iend = which.max(sub_ss[, 2])
        while (sub_ss[istart, 2] > 0 & istart > 1) istart <- istart - 
          1
        while (sub_ss[iend, 2] > 0 & iend < nrow(sub_ss)) iend <- iend + 
          1
        # if ((iend - istart) < npoint) {
        #   warning(sprintf("Less then %d datapoints in peak at mz=%.4f",
        #                   npoint, j_tmt))
        #   return(data.frame())
        # }
        sub_ss <- sub_ss[istart:iend, ]
        I = sub_ss[which.max(sub_ss[, 2]), 2] * sd(sub_ss[, 
                                                          1])/2
        mu = sub_ss[which.max(sub_ss[, 2]), 1]
        
        sigma = sd(sub_ss[, 
                          1])
        fit2P = fitG(sub_ss[, 1], sub_ss[, 2],mu,sigma, I)
        p = fit2P$par
        tmt_RSq <- c(tmt_RSq, rSquared_tmt(p, y = sub_ss[,2], x = sub_ss[,1]))
        
        if (abs(rSquared_tmt(p, y = sub_ss[,2], x = sub_ss[,1])) >= 0.9){
          tmt_count = tmt_count+1
          sum_tmt = sum_tmt + max(sub_ss[, 2])
          
          if (max(sub_ss[,2])>=max_tmt){
            max_tmt = max(sub_ss[,2])
          }
          
        }
        
      }
      # ires <- ires %>% bind_rows(ires %>% mutate(
      #                            Detected_TMT = tmt_count, Max_TMT_Intensity = max_tmt,
      #                            Total_TMT_Intensity = sum_tmt, TMT_R2 = median(tmt_RSq)))
      #                            
      ires <- ires %>% mutate(Detected_TMT = tmt_count, Max_TMT_Intensity = max_tmt, Total_TMT_Intensity = sum_tmt, TMT_R2 = median(tmt_RSq))
      
      ires$n <- apply(ires, 1, function(x) ions[[x["ion"]]][sub("\\d+", 
                                                                "", x["peak"])])
      ires
    })
  result
}


mzMLtoCSV_Contamination = function (pattern = "*.mzML(.gz)?", width = 0.001) 
{
  library(isoms)
  library(parallel)
  ncpu = detectCores()
  args_ <- commandArgs(trailingOnly = TRUE)
  dir_ = "."
  if (length(args_) > 0) 
    dir_ <- args_[[1]]
  if (length(args_) > 1) 
    pattern <- args_[[2]]
  if (!dir.exists(dir_)) {
    message("First argument has to be the directory containing mzML files")
    return(0)
  }
  ff <- list.files(path = dir_, pattern = pattern)
  if (length(ff) < 1) {
    message("Now mzML files found in directory [", dir_, 
      "] using pattern [", pattern, "].")
    return(0)
  }
  message(sprintf("Processing %d files. You can go and have some beer meanwhile", 
    length(ff)))
  if (Sys.info()["sysname"] == "Windows") {
    cl <- makeCluster(ncpu)
    clusterExport(cl, c("dir_", "ff"), environment())
    result <- bind_rows(parLapply(cl, ff, function(ii) {
      library(isoms)
      in_f <- file.path(dir_, ii)
      outf <- sub(".mzML(.gz)?$", "_fit.csv", in_f)
      if (!file.exists(outf)) {
        xxx <- analyze_Immonium_TMT(in_f, ions = immoniumIons, 
          width = width, fixSigma = T)
        write_csv(xxx, outf)
      }
      else xxx <- suppressMessages(read_csv(outf))
      xxx %>% mutate(file = in_f) %>% group_by(file) %>% 
        summarise(nScans = length(unique(seqNum)))
    }))
    stopCluster(cl)
  }
  else {
    result <- bind_rows(mclapply(ff, function(ii) {
      in_f <- file.path(dir_, ii)
      outf <- sub(".mzML(.gz)?$", "_fit.csv", in_f)
      if (!file.exists(outf)) {
        xxx <- analyze_Immonium_TMT(in_f, ions = immoniumIons, 
          width = width, fixSigma = T)
        write_csv(xxx, outf)
      }
      else xxx <- suppressMessages(read_csv(outf))
      xxx %>% mutate(file = in_f) %>% group_by(file) %>% 
        summarise(nScans = length(unique(seqNum)))
    }, mc.cores = ncpu, mc.preschedule = TRUE))
  }
  return(result)
}


