mzMLtoCSV <- function(pattern = "*.mzML(.gz)?", width=0.001) {
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
        message("Now mzML files found in directory [", dir_, "] using pattern [",
            pattern, "].")
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
                xxx <- analyze_immoniums(in_f, ions = immoniumIons, width = width,
                  fixSigma = T)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }))
        stopCluster(cl)
    } else {
        result <- bind_rows(mclapply(ff, function(ii) {
            in_f <- file.path(dir_, ii)
            outf <- sub(".mzML(.gz)?$", "_fit.csv", in_f)
            if (!file.exists(outf)) {
                xxx <- analyze_immoniums(in_f, ions = immoniumIons, width = width,
                  fixSigma = T)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }, mc.cores = ncpu, mc.preschedule = TRUE))
    }
    return(result)
    # for(f in ff){ in_f <- file.path(dir_, f) outf <- sub('.mzML$','_fit.csv', in_f)
    # if(!file.exists(outf)){ message('Converting file [', f, ']')
    # analyze_immoniums(file = in_f, width = 0.0015) %>% write_csv(outf) }else{
    # message('File already processed: [', f, '].') } }
}

mzMLtoCSV2 <- function(pattern = "*.mzML(.gz)?", width=0.0015) {
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
        message("Now mzML files found in directory [", dir_, "] using pattern [",
            pattern, "].")
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
            outf <- sub(".mzML(.gz)?$", "_fit2.csv", in_f)
            if (!file.exists(outf)) {
                xxx <- analyze_immoniums2(in_f, ions = immoniumIons, width = width)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }))
        stopCluster(cl)
    } else {
        result <- bind_rows(mclapply(ff, function(ii) {
            in_f <- file.path(dir_, ii)
            outf <- sub(".mzML(.gz)?$", "_fit2.csv", in_f)
            if (!file.exists(outf)) {
                xxx <- analyze_immoniums2(in_f, ions = immoniumIons, width = width)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }, mc.cores = ncpu, mc.preschedule = TRUE))
    }
    return(result)
    # for(f in ff){ in_f <- file.path(dir_, f) outf <- sub('.mzML$','_fit.csv', in_f)
    # if(!file.exists(outf)){ message('Converting file [', f, ']')
    # analyze_immoniums(file = in_f, width = 0.0015) %>% write_csv(outf) }else{
    # message('File already processed: [', f, '].') } }
}

mzMLtoCSV_ <- function(pattern = "*.mzML(.gz)?", width=0.001) {
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
        message("Now mzML files found in directory [", dir_, "] using pattern [",
            pattern, "].")
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
                xxx <- analyze_immoniums_(in_f, ions = AAIons, width = width,
                  fixSigma = T)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }))
        stopCluster(cl)
    } else {
        result <- bind_rows(mclapply(ff, function(ii) {
            in_f <- file.path(dir_, ii)
            outf <- sub(".mzML(.gz)?$", "_fit.csv", in_f)
            if (!file.exists(outf)) {
                xxx <- analyze_immoniums_(in_f, ions = AAIons, width = width,
                  fixSigma = T)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }, mc.cores = ncpu, mc.preschedule = TRUE))
    }
    return(result)
}

mzMLtoCSV_2 <- function(pattern = "*.mzML(.gz)?", width=0.001) {
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
        message("Now mzML files found in directory [", dir_, "] using pattern [",
            pattern, "].")
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
                xxx <- analyze_immoniums_2(in_f, ions = AAIons, width = width,
                  fixSigma = T)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }))
        stopCluster(cl)
    } else {
        result <- bind_rows(mclapply(ff, function(ii) {
            in_f <- file.path(dir_, ii)
            outf <- sub(".mzML(.gz)?$", "_fit.csv", in_f)
            if (!file.exists(outf)) {
                xxx <- analyze_immoniums_2(in_f, ions = AAIons, width = width,
                  fixSigma = T)
                write_csv(xxx, outf)
            } else # message(sprintf('File already processed: %s', outtile))
            xxx <- suppressMessages(read_csv(outf))
            xxx %>% mutate(file = in_f) %>% group_by(file) %>% summarise(nScans = length(unique(seqNum)))
        }, mc.cores = ncpu, mc.preschedule = TRUE))
    }
    return(result)
}
