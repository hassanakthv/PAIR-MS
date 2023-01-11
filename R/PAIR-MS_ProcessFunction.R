
iion <- c("P", "Hyp", "L/I", "R", "R2", "K", "Hyl", "A", "H", "F", "S","V", "Y", "N","N2","D","Q","E","T","Pyr","m57","K2","m89",
          "m73","m87","M","C","camC","xC","W")


PAIRMS <- function (file = "experimentDesign.csv", data_tab, outdir, ioi = list(iointrest), 
          ioi_ = list(iontoshow), correct = TRUE) 
{
  iion <- c("P", "Hyp", "L/I", "R", "R2", "K", "Hyl", "A", "H", "F", "S","V", "Y", "N","N2","D","Q","E","T","Pyr","m57","K2","m89",
            "m73","m87","M","C","camC","xC","W")
  library(isoms)
  args_ <- commandArgs(trailingOnly = TRUE)
  if (length(args_) > 0) {
    expDesignFile <- args_[[1]]
    if (!file.exists(expDesignFile)) {
      if (file.exists(file)) {
        warning(sprintf("File [%s] doesn't exist, using [%s] instead", 
                        expDesignFile, file))
        expDesignFile <- file
      }
      else {
        message("You must provide corrent experimentDesign file location")
        return(0)
      }
    }
  }
  else {
    if (file.exists(file)) 
      expDesignFile <- file
    else {
      message("You must provide corrent experimentDesign file location")
      return(0)
    }
  }
  expDesign <- suppressMessages(read_csv(expDesignFile))
  message("Samples to analyze: ", toString(unique(expDesign$Sample)))
  message("Number of controls: ", nrow(expDesign %>% filter(Sample == 
                                                              "Control")))

  outdir <- outdir
  data_ <- expDesign %>% rowwise() %>% do({
    xx = .
    fd <- suppressWarnings(suppressMessages(read_csv(as.character(data_tab$Path[which(data_tab$File==xx$File)]), 
                                                     col_types = cols(tic = col_double(), bpI = col_double(), 
                                                                      ltic = col_double()))))
    if ("totIonCurrent" %in% names(fd)) 
      fd <- fd %>% mutate(tic = totIonCurrent) %>% select(-totIonCurrent)
    fd %>% mutate(I = abs(I), sigma = abs(sigma)) %>% mutate(file = xx$File, 
                                                             group = xx$Sample, loading = xx$Loading) %>% filter(rt > 
                                                                                                                   (xx$Start * 60) & (rt < xx$End * 60))
  }) %>% ungroup() %>% mutate(I0 = I/isoratio) %>% mutate(dI = I0/tic, 
                                                          ldI = log10(dI)) %>% mutate(ltic = log10(tic))
 
  x_html <- SummarizefileOnline(data = data_, group = "group", resultPath = outdir, 
                 correct = correct, IOI = ioi, Info = expDesign[1, 10:ncol(expDesign)], 
                 IOI_ = ioi_)
  return(x_html)
}


SummarizefileOnline <- function (data = NA, files = NA, group = ifelse(is.na(data), 
                                                                       "file", NA), resultPath = "./isoMS_result", correct = T, 
                                 IOI = NA, Info = NA, IOI_ = NA) 
{
  if (is.na(data)) {
    if (class(files) == "character") 
      data <- bind_rows(lapply(files, function(f) {
        if (file.exists(f)) 
          read_csv(f) %>% mutate(file = f)
        else data.frame()
      }))
  }
  if (nrow(data) < 2) {
    warning("You have to provide either data frame or vecor of file names to proceed")
    return(0)
  }
  if (!dir.exists(resultPath)) 
    dir.create(resultPath, recursive = TRUE)
  mass_tol <- 3e-04
  if (!("file" %in% names(data))) {
    data <- data %>% mutate(file = "file")
  }
  if (is.na(group)) 
    group = "group"
  if (!(group %in% names(data))) {
    data <- data %>% mutate(group = "group")
  }
  else if (group != "group") {
    data <- data %>% mutate_(group = sprintf("`%s`", group))
  }
  message("Files to process: ")
  print(unique(data$file))
  message("Groups to process: ")
  print(unique(data$group))
  spectra <- unique(data$seqNum)
  mono_spectra <- unique((data %>% filter(peak == "0"))$seqNum)
  nomono_spectra <- setdiff(spectra, mono_spectra)
  data <- data %>% filter(ion %in% IOI)
  ranges <- data %>% group_by(file, peak, ion) %>% summarise(meanerror = mean(masserror), 
                                                             medianerror = median(masserror), meansq = mean(masserror^2), 
                                                             gg = median(isoratio/n, na.rm = T), ir = median(isoratio, 
                                                                                                             na.rm = T), md = mad(isoratio/n, na.rm = T), nhits = n(), 
                                                             mmin = gg - 3 * md, mmax = gg + 3 * md)
  data <- data %>% mutate(g_ = isoratio/n)
  hf <- unique(data$file)
  hh <- data.frame()
  for (i in hf) {
    htemp <- data %>% filter(file == i)
    hi <- unique(htemp$ion)
    for (j in hi) {
      htemp_ <- htemp %>% filter(ion == j)
      hp <- unique(htemp_$peak)
      for (k in hp) {
        ref <- ranges %>% filter(file == i & ion == 
                                   j & peak == k)
        hh_ <- htemp_ %>% filter(peak == k) %>% filter(peak == 
                                                         "0" | g_ > ref$mmin & g_ < ref$mmax & R2 > 
                                                         0.9)
        hh <- rbind(hh, hh_)
      }
    }
  }
  data <- hh
  rranges <- ranges %>% filter(nhits >= 2 & gg > 0) %>% summarize(mgg = median(gg, 
                                                                               na.rm = T), madgg = mad(gg, na.rm = T), mmd = min(md))
  ranges <- ranges %>% ungroup() %>% left_join(rranges, by = c("file", 
                                                               "peak")) %>% mutate(isgood = peak == "0" | ((md/ir < 
                                                                                                              1) & (meansq < 1e-06) & (abs(meanerror) < mass_tol))) %>% 
    select(-mgg, -madgg, -mmd)
  data <- data %>% ungroup() %>% filter(abs(masserror) < mass_tol) %>% 
    left_join(ranges, by = c("file", "ion", "peak")) %>% 
    filter(!is.na(isgood)) %>% filter(isgood & (peak == 
                                                  "0" | (isoratio/n > mmin & isoratio/n < mmax)))
  data_ldI <- data %>% group_by(peak, ion) %>% filter(ldI > 
                                                        (median(ldI) - 2 * mad(ldI))) %>% mutate(ldI = ldI - 
                                                                                                   median(ldI)) %>% mutate(gg = isoratio/n) %>% ungroup()
  ldI_model <- data_ldI %>% filter(peak != "0") %>% group_by(peak, 
                                                             ion) %>% do({
                                                               dd <- .
                                                               peak <- dd$peak[[1]]
                                                               ion <- dd$ion[[1]]
                                                               if (length(unique(dd$group)) > 1) {
                                                                 tidy(lm(gg ~ ldI + group, data = .))
                                                               }
                                                               else tidy(lm(gg ~ ldI, data = .))
                                                             }) %>% filter(term == "ldI")
  data_ldI <- data_ldI %>% left_join(ldI_model, by = c("peak", 
                                                       "ion")) %>% mutate(gg = gg - ldI * estimate)
  if (correct) {
    data <- data_ldI %>% mutate(gamma = gg)
  }
  else {
    data <- data %>% mutate(gamma = isoratio/n)
  }
  # if (file.exists(file.path(resultPath, "all_aa.RData"))) {
  #   message("all_aa.RData file exists, will not recalculate it.")
  #   load(file.path(resultPath, "all_aa.RData"))
  # }
  #else {
    message("Computing 'any' ion: ")
    gpb = dplyr::group_by
    rw <- dplyr::rowwise
    if ("multidplyr" %in% installed.packages()) {
      library(multidplyr)
      cl <- create_cluster()
      cluster_library(cl, "dplyr")
      cluster_library(cl, "isoms")
      set_default_cluster(cl)
      gpb <- partition
      rw <- partition
    }
    any_aa <- data %>% gpb(group, file, seqNum, rt) %>% 
      do({
        dd <- .
        tic_ <- dd$tic[[1]]
        ltic_ <- dd$ltic[[1]]
        res <- data.frame()
        peaks <- setdiff(unique(dd$peak), "0")
        for (el in peaks) {
          good_aa <- (dd %>% filter(peak == el & I > 
                                      0) %>% distinct(ion))$ion
          i0s <- (dd %>% filter(ion %in% good_aa & peak == 
                                  "0"))$I
          i1s <- (dd %>% filter(ion %in% good_aa & peak == 
                                  el))$isoratio
          i1s <- i1s * i0s
          ns <- (dd %>% filter(ion %in% good_aa & peak == 
                                 el))$n
          rs <- i1s/i0s
          i0s <- i0s[rs < 1]
          i1s <- i1s[rs < 1]
          ns <- ns[rs < 1]
          rs <- rs[rs < 1]
          rsn <- rs/ns
          rs_good <- (abs(rsn - weighted.mean(rsn, i0s)) < 
                        2.5 * mad(rsn))
          r <- sum(i1s[rs_good]/ns[rs_good])/sum(i0s[rs_good])
          if (!is.na(r)) 
            res <- res %>% bind_rows(data.frame(peak = el, 
                                                gamma = r, logI = log2(sum(i0s[rs_good])), 
                                                tic = tic_, dI = sum(i0s[rs_good])/tic_, 
                                                ldI = log10(sum(i0s[rs_good])/tic_), ltic = ltic_))
        }
        if (nrow(res) > 0) {
          res$ion = "any"
        }
        res
      }) %>% collect() %>% ungroup()
    message("Puting ions together: ")
    all_aa <- data %>% filter(n > 0 & I > 0) %>% mutate(logI = log2(I)) %>% 
      select(group, file, peak, gamma, logI, tic, ltic, 
             dI, ldI, I0, seqNum, rt, ion) %>% bind_rows(any_aa) %>% 
      ungroup()
    save(any_aa, all_aa, Info, file = file.path(resultPath, 
                                                "all_aa.RData"))

  message("Grouping results and writing output tables: ")
  gg <- all_aa %>% distinct(group)
  message("Rendering HTML report: ")
  resultPath <- normalizePath(resultPath)
  rmarkdown::render(system.file("rmd", "isoMS_report.Rmd", package = "isoms"),
                    envir = sys.frame(sys.nframe()), output_file = file.path(resultPath,
                                                                             "isoMS_report.html"))
  g_plot <- ggplot(data = data %>% dplyr::filter(peak != 0 & ion=="Hyp"), aes(x = log10(tic), y = gamma*100, fill = ion)) + geom_point(shape = 21, color = "black") +
    theme_bw() + facet_wrap(ion~peak, scales = "free")
 
  return(data)
  
}

