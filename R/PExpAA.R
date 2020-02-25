PExpAA <- function (file = "experimentDesign.csv", ioi = list(iointrest), ioi_ = list(iontoshow),
          correct = TRUE) 
{
  library(isoms)
  args_ <- commandArgs(trailingOnly = TRUE)
    if(length(args_)>0){
      expDesignFile <- args_[[1]]
      if (!file.exists(expDesignFile)){
        if(file.exists(file)){
          warning(sprintf("File [%s] doesn't exist, using [%s] instead", expDesignFile, file))
          expDesignFile <- file
        } else {
          message("You must provide corrent experimentDesign file location")
          return(0)
        }
      }
    } else{
      if(file.exists(file))
        expDesignFile <- file
      else{
        message("You must provide corrent experimentDesign file location")
        return(0)
      }
    }
          
          
  expDesign <- suppressMessages(read_csv(expDesignFile))
  message("Samples to analyze: ", toString(unique(expDesign$Sample)))
  message("Number of controls: ", nrow(expDesign %>% filter(Sample == 
                                                              "Control")))
  outdir <- file.path(dirname(normalizePath(expDesignFile)), 
                      sub(".csv$", "", basename(normalizePath(expDesignFile))))
  data_ <- expDesign %>% rowwise() %>% do({
    xx = .
    fd <- suppressWarnings(suppressMessages(read_csv(xx$File, 
                                                     col_types = cols(tic = col_double(), bpI = col_double(), 
                                                                      ltic = col_double()))))
    if ("totIonCurrent" %in% names(fd)) 
      fd <- fd %>% mutate(tic = totIonCurrent) %>% select(-totIonCurrent)
    fd %>% mutate(I = abs(I), sigma = abs(sigma)) %>% mutate(file = xx$File, 
                                                             group = xx$Sample, loading = xx$Loading) %>% filter(rt > 
                                                                                                                   (xx$Start * 60) & (rt < xx$End * 60))
  }) %>% ungroup() %>% mutate(I0 = I/isoratio) %>% mutate(dI = I0/tic, 
                                                          ldI = log10(dI)) %>% mutate(ltic = log10(tic)) %>%
                                                          filter(group==ion)
  SummarizeAA(data = data_, group = "group", resultPath = outdir, 
                     correct = correct, IOI = ioi, Info = expDesign[1, 10:ncol(expDesign)], IOI_ = ioi_)
}
