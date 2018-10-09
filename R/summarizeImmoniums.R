#' summarizeImmoniums - generate the summary of the immonium analysis
#' Given the results of \code{analyzeImmoniums} function, produces a nice html
#' report
#'
#' @param data a dataframe returned by \code{analyzeImmoniums} function.
#' If \code{NA} then using \code{files}
#' @param files vector of csv-files containing results of \code{analyzeImmoniums}
#' @param group name of the column in data
#'
#' @return
#' @export
#'
#' @examples
summarizeImmoniums <- function(data = NA,
                               files = NA,
                               group = ifelse(is.na(data),
                                              "file",
                                              NA),
                               resultPath = "./isoMS_result",
                               correct=T) {
    if (is.na(data)) {
        if (class(files) == "character")
            data <- bind_rows(lapply(files, function(f) {
                if (file.exists(f))
                  read_csv(f) %>% mutate(file = f) else data.frame()
            }))
    }
    if (nrow(data) < 2) {
        warning("You have to provide either data frame or vecor of file names to proceed")
        return(0)
    }
    if (!dir.exists(resultPath))
        dir.create(resultPath, recursive = TRUE)

    if ("totIonCurrent" %in% names(data))
        data <- data %>% mutate(tic = totIonCurrent)
    if ("tic" %in% names(data)) {
        message("TIC adjustment")
        TIClimits <- data %>% group_by(file) %>% summarize(minTIC = min(tic, na.rm = T),
            maxTIC = max(tic, na.rm = T)) %>% summarize(low = max(c(2e+07, minTIC)),
            high = min(maxTIC))
        data <- data %>% filter(tic > TIClimits$low[1]/4 & tic < TIClimits$high[1] *
            4)
    } else {
        data$tic <- 10
    }
    mass_tol <- 3e-04
    if (!("file" %in% names(data))) {
        data <- data %>% mutate(file = "file")
    }
    if (is.na(group))
        group = "group"
    if (!(group %in% names(data))) {
        data <- data %>% mutate(group = "group")
    } else if (group != "group") {
        data <- data %>% mutate_(group = sprintf("`%s`", group))
    }
    message("Files to process: ")
    print(unique(data$file))
    message("Groups to process: ")
    print(unique(data$group))
    spectra <- unique(data$seqNum)
    mono_spectra <- unique((data %>% filter(peak == "0"))$seqNum)
    nomono_spectra <- setdiff(spectra, mono_spectra)



    # removing outliers
    ranges <- data %>%
      group_by(file, peak, ion) %>%
      summarise(meanerror = mean(masserror),
                medianerror = median(masserror),
                meansq = mean(masserror^2),
                gg = median(isoratio/n, na.rm = T),
                ir = median(isoratio, na.rm = T),
                md = mad(isoratio/n, na.rm = T),
                nhits = n(), mmin = gg - 3 * md, mmax = gg + 3 * md)

    rranges <- ranges %>%
      filter(nhits > 10 & gg > 0) %>%
      summarize(mgg = median(gg, na.rm = T),
                madgg = mad(gg, na.rm = T),
                mmd = min(md))
    ranges <- ranges %>%
      ungroup() %>%
      left_join(rranges, by = c("file", "peak")) %>%
      mutate(isgood = peak == "0" | ((md < 8 * mmd) & (md/ir < 1) & (meansq < 1e-06) &
        (abs(meanerror) < mass_tol))) %>% select(-mgg, -madgg, -mmd)

    data <- data %>%
      ungroup() %>%
      filter(abs(masserror) < mass_tol) %>%
      left_join(ranges, by = c("file", "ion", "peak")) %>%
      filter(!is.na(isgood)) %>% filter(isgood &
        (peak == "0" | (isoratio/n > mmin & isoratio/n < mmax)))
    data %>%
      group_by(peak, ion) %>%
      filter(ldI > (median(ldI) - 2*mad(ldI))) %>%
      mutate(ldI = ldI-median(ldI)) %>%
      mutate(gg=isoratio/n) %>% ungroup() -> data_ldI


    data_ldI %>%
      filter(peak != '0') %>%
      group_by(peak, ion) %>%
      do({
        dd <- .
        peak <- dd$peak[[1]]
        ion <- dd$ion[[1]]
        if(length(unique(dd$group))>1){
          tidy(lm(gg~ldI+group, data=.))
        } else
          tidy(lm(gg~ldI, data=.))
      }) %>% filter(term=='ldI') -> ldI_model

    data_ldI <- data_ldI %>%
      left_join(ldI_model, by=c('peak', 'ion')) %>%
      mutate(gg=gg-ldI*estimate)

    if(correct){
      data <- data_ldI %>%  # mutate(isoratio = gg*n) %>%
        mutate(gamma = gg)
    } else {
      data <- data %>% mutate(gamma = isoratio/n)
    }

    if (file.exists(file.path(resultPath, "all_aa.RData"))) {
        message("all_aa.RData file exists, will not recalculate it.")
        load(file.path(resultPath, "all_aa.RData"))
    } else {
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

        any_aa <- data %>% gpb(group, file, seqNum, rt) %>% do({
            dd <- .
            tic_ <- dd$tic[[1]]
            ltic_ <- dd$ltic[[1]]

            res <- data.frame()
            peaks <- setdiff(unique(dd$peak), "0")
            for (el in peaks) {
                good_aa <- (dd %>% filter(peak == el & I > 0) %>% distinct(ion))$ion
                i0s <- (dd %>% filter(ion %in% good_aa & peak == "0"))$I
                i1s <- (dd %>% filter(ion %in% good_aa & peak == el))$isoratio
                i1s <- i1s * i0s
                ns <- (dd %>% filter(ion %in% good_aa & peak == el))$n

                rs <- i1s/i0s
                i0s <- i0s[rs < 1]
                i1s <- i1s[rs < 1]
                ns <- ns[rs < 1]
                rs <- rs[rs < 1]
                rsn <- rs/ns
                rs_good <- (abs(rsn - weighted.mean(rsn, i0s)) < 2.5 * mad(rsn))
                r <- sum(i1s[rs_good]/ns[rs_good])/sum(i0s[rs_good])
                if (!is.na(r))
                  res <- res %>% bind_rows(data.frame(peak = el, gamma = r, logI = log2(sum(i0s[rs_good])),
                    tic = tic_, dI = sum(i0s[rs_good])/tic_, ldI = log10(sum(i0s[rs_good])/tic_),
                    ltic = ltic_))
            }
            if (nrow(res) > 0) {
                # res$seqNum = dd$seqNum[[1]] res$rt = dd$rt[[1]]
                res$ion = "any"
            }
            res
        }) %>% collect() %>% ungroup()

        message("Puting ions together: ")

        all_aa <- data %>% filter(n > 0 & I > 0) %>% mutate(logI = log2(I)) %>%
            select(group, file, peak, gamma, logI, tic, ltic, dI, ldI, I0, seqNum,
                rt, ion) %>% bind_rows(any_aa) %>% ungroup()
        save(any_aa, all_aa, file = file.path(resultPath, "all_aa.RData"))
        message("Performing LOESS method: ")
        cluster_copy(cl, all_aa)
        loess_C <- all_aa %>% ungroup() %>% distinct(group, file, ion) %>% arrange(ion, group,
            file) %>% rowwise() %>% # rw() %>%
          do({
            ion_ = .$ion[[1]]
            file_ = .$file[[1]]
            group_ = .$group[[1]]
            quantifyIsoRatio(data = all_aa, control = "Control", file = file_, ion = ion_,
                peak = "13C") %>% mutate(file = file_, ion = ion_, group = group_)
        }) %>% collect() %>% ungroup()
        loess_N <- all_aa %>% ungroup() %>% distinct(group, file, ion) %>% arrange(ion,
            group, file) %>% rowwise() %>% do({
            ion_ = .$ion[[1]]
            file_ = .$file[[1]]
            group_ = .$group[[1]]
            quantifyIsoRatio(data = all_aa, control = "Control", file = file_, ion = ion_,
                peak = "15N") %>% mutate(file = file_, ion = ion_, group = group_)
        }) %>% collect() %>% ungroup()
        loess_H <- all_aa %>% ungroup() %>% distinct(group, file, ion) %>% arrange(ion,
            group, file) %>% rowwise() %>% do({
            ion_ = .$ion[[1]]
            file_ = .$file[[1]]
            group_ = .$group[[1]]
            quantifyIsoRatio(data = all_aa, control = "Control", file = file_, ion = ion_,
                peak = "2H") %>% mutate(file = file_, ion = ion_, group = group_)
        }) %>% collect() %>% ungroup()
        loess_O <- all_aa %>% ungroup() %>% distinct(group, file, ion) %>% arrange(ion,
            group, file) %>% rowwise() %>% do({
            ion_ = .$ion[[1]]
            file_ = .$file[[1]]
            group_ = .$group[[1]]
            quantifyIsoRatio(data = all_aa, control = "Control", file = file_, ion = ion_,
                peak = "18O") %>% mutate(file = file_, ion = ion_, group = group_)
        }) %>% collect() %>% ungroup()
        save(any_aa, all_aa, loess_C, loess_N, loess_H, loess_O, file = file.path(resultPath,
            "all_aa.RData"))
    }

    message("Grouping results and writing output tables: ")

    gg <- all_aa %>% distinct(group)
    message("Rendering HTML report: ")
    resultPath <- normalizePath(resultPath)
    rmarkdown::render(system.file("Rmd/isoMS_report.Rmd", package = getPackageName()),
        envir = sys.frame(sys.nframe()), output_file = file.path(resultPath, "isoMS_report.html"))
    rmarkdown::render(system.file("Rmd/isoMS_loess.Rmd", package = getPackageName()),
        envir = sys.frame(sys.nframe()), output_file = file.path(resultPath, "isoMS_loess.html"))
}

summarizeImmoniums2 <- function(data = NA, files = NA, group = ifelse(is.na(data),
    "file", NA), resultPath = "./isoMS_result") {
    if (is.na(data)) {
        if (class(files) == "character")
            data <- bind_rows(lapply(files, function(f) {
                if (file.exists(f))
                  read_csv(f) %>% mutate(file = f) else data.frame()
            }))
    }
    if (nrow(data) < 2) {
        warning("You have to provide either data frame or vecor of file names to proceed")
        return(0)
    }
    if (!dir.exists(resultPath))
        dir.create(resultPath, recursive = TRUE)

    if ("totIonCurrent" %in% names(data)) {
        message("TIC adjustment")
        TIClimits <- data %>% group_by(file) %>% summarize(minTIC = min(totIonCurrent,
            na.rm = T), maxTIC = max(totIonCurrent, na.rm = T)) %>% summarize(low = max(c(2e+07,
            minTIC)), high = min(maxTIC))
        data <- data %>% filter(totIonCurrent > TIClimits$low & totIonCurrent < TIClimits$high)
    }
    mass_tol <- 3e-04

    if (!("file" %in% names(data))) {
        data <- data %>% mutate(file = "file")
    }
    if (is.na(group))
        group = "group"
    if (!(group %in% names(data))) {
        data <- data %>% mutate(group = "group")
    } else if (group != "group") {
        data <- data %>% mutate_(group = sprintf("`%s`", group))
    }
    message("Files to process: ")
    print(unique(data$file))
    message("Groups to process: ")
    print(unique(data$group))

    # removing outliers
    ranges <- data %>% group_by(file, element, ion) %>% summarize(gg = median(estimate/n,
        na.rm = T), ir = median(estimate, na.rm = T), md = mad(estimate, na.rm = T),
        nhits = n(), mmin = ir - 3 * md, mmax = ir + 3 * md)
    rranges <- ranges %>% filter(nhits > 10 & gg > 0) %>% summarize(mgg = median(gg,
        na.rm = T), madgg = mad(gg, na.rm = T), mmd = min(md))
    ranges <- ranges %>% ungroup() %>% left_join(rranges, by = c("file", "element")) %>%
        # mutate(isgood = ((abs(gg-mgg) < 5*madgg)&(md < 4*mmd)&(md/ir<1))) %>%
    mutate(isgood = ((abs(log2(gg/mgg)) < 2) & (md < 6 * mmd) & (md/ir < 1))) %>%
        select(-mgg, -madgg, -mmd)
    data <- data %>% ungroup() %>% left_join(ranges, by = c("file", "ion", "element")) %>%
        filter(!is.na(isgood)) %>% filter(isgood & (estimate > mmin & estimate <
        mmax))

    message("Linear-model based method ")
    lm_res_aa <- data %>% filter(n > 0 & I > 10) %>% mutate(I0 = I/100, I = I0 *
        estimate/n * 100) %>% group_by(element, group, file, ion) %>% do({
        mm <- try(MASS::rlm(I ~ I0, data = .))
        if (class(mm) != "try-error") {
            tidy(mm)
        } else {
            data.frame()
        }

    }) %>% mutate(peak = element) %>% mutate(peak = sub("C", "13C", peak)) %>% mutate(peak = sub("N",
        "15N", peak)) %>% mutate(peak = sub("H", "2H", peak))
    lm_res_aa %>% write_csv(file.path(resultPath, sprintf("linearmodel2_summary_aa.csv")))
    lm_res <- data %>% filter(n > 0 & I > 10) %>% mutate(I0 = I/100, I = I0 * estimate/n *
        100) %>% group_by(element, group, file) %>% do({
        mm <- try(MASS::rlm(I ~ I0 + ion, data = ., method = "MM"))
        if (class(mm) != "try-error") {
            tidy(mm)
        } else {
            data.frame()
        }
    }) %>% mutate(peak = element) %>% mutate(peak = sub("C", "13C", peak)) %>% mutate(peak = sub("N",
        "15N", peak)) %>% mutate(peak = sub("H", "2H", peak))
    lm_res %>% write_csv(file.path(resultPath, sprintf("linearmodel2_summary.csv")))


    if (file.exists(file.path(resultPath, "all_aa2.RData"))) {
        message("all_aa.RData file exists, will not recalculate it.")
        load(file.path(resultPath, "all_aa2.RData"))
    } else {
        message("Computing 'any' ion: ")
        gpb = dplyr::group_by
        if ("multidplyr" %in% installed.packages()) {
            library(multidplyr)
            cl <- create_cluster()
            cluster_library(cl, "dplyr")
            set_default_cluster(cl)
            gpb <- partition
        }

        any_aa <- data %>% gpb(group, file, seqNum, rt) %>% do({
            dd <- .
            res <- data.frame()
            elements <- setdiff(unique(dd$element), "0")
            for (el in elements) {
                good_aa <- (dd %>% filter(element == el & I > 0) %>% distinct(ion))$ion
                i0s <- (dd %>% filter(ion %in% good_aa & element == el))$I
                i1s <- (dd %>% filter(ion %in% good_aa & element == el) %>% mutate(xx = I *
                  estimate))$xx
                ns <- (dd %>% filter(ion %in% good_aa & element == el))$n

                rs <- i1s/i0s
                i0s <- i0s[rs < 1]
                i1s <- i1s[rs < 1]
                ns <- ns[rs < 1]
                rs <- rs[rs < 1]
                rs_good <- (abs(rs - median(rs)) < 1.5 * mad(rs))
                r <- sum(i1s[rs_good]/ns[rs_good])/sum(i0s[rs_good])
                if (!is.na(r))
                  res <- res %>% bind_rows(data.frame(element = el, gamma = r, logI = log2(sum(i1s))))
            }
            if (nrow(res) > 0) {
                # res$seqNum = dd$seqNum[[1]] res$rt = dd$rt[[1]]
                res$ion = "any"
            }
            res
        }) %>% collect() %>% ungroup()

        message("Puting ions together: ")

        all_aa <- data %>% filter(n > 0 & I > 0) %>% mutate(gamma = estimate/n, logI = log2(I)) %>%
            select(group, file, element, gamma, logI, seqNum, rt, ion) %>% bind_rows(any_aa) %>%
            ungroup()
        save(any_aa, all_aa, file = file.path(resultPath, "all_aa2.RData"))
    }
    message("Grouping results and writing output tables: ")

    for (pp in unique(all_aa$element)) {
        vals <- (all_aa %>% filter(element == pp))$gamma
        vmin <- quantile(vals, 0.05)
        vmax <- quantile(vals, 0.95)
        mmad <- all_aa %>% filter(element == pp) %>% group_by(ion, file) %>% summarize(m = median(gamma),
            md = mad(gamma))
        aa_data <- all_aa %>% filter(element == pp) %>% left_join(mmad %>% ungroup(),
            by = c("ion", "file")) %>% filter(md > 0) %>% filter(gamma > (m - 3 *
            md) & gamma < (m + 3 * md))
        aa_data$ion <- factor(aa_data$ion, levels = (mmad %>% summarise(md = median(md)) %>%
            arrange(md))$ion)
        res_f <- aa_data %>% group_by(ion, group, file) %>% summarize(mean = mean(gamma) *
            100, median = median(gamma) * 100, wmean = weighted.mean(gamma, logI) *
            100, sd = sd(gamma) * 100, se = sd(gamma)/sqrt(n()) * 100, n = n()) %>%
            mutate(CV = se/wmean * 100) %>% arrange(ion)
        res_f %>% write_csv(file.path(resultPath, sprintf("file_summary_%s.csv",
            pp)))
        res_f %>% summarize(gmean = mean(mean), gmedian = mean(median), gwmean = mean(wmean),
            sd_mean = sd(mean), sd_median = sd(median), sd_wmean = sd(wmean), n = n()) %>%
            mutate(cv_mean = sd_mean/gmean/sqrt(n) * 100, cv_median = sd_median/gmedian/sqrt(n) *
                100, cv_wmean = sd_wmean/gwmean/sqrt(n) * 100) %>% arrange(cv_wmean) %>%
            write_csv(file.path(resultPath, sprintf("group_summary2_%s.csv", pp)))
    }
    all_aa <- all_aa %>% mutate(peak = element) %>% mutate(peak = sub("C", "13C",
        peak)) %>% mutate(peak = sub("N", "15N", peak)) %>% mutate(peak = sub("H",
        "2H", peak))

    gg <- all_aa %>% distinct(group)
    message("Rendering HTML report: ")
    resultPath <- normalizePath(resultPath)
    rmarkdown::render(system.file("Rmd/isoMS_report.Rmd", package = getPackageName()),
        envir = sys.frame(sys.nframe()), output_file = file.path(resultPath, "isoMS_report.html"))
}

summarizeImmoniumsLoess <- function(data = NA, files = NA, group = ifelse(is.na(data),
    "file", NA), resultPath = "./isoMS_result") {
    if (is.na(data)) {
        if (class(files) == "character")
            data <- bind_rows(lapply(files, function(f) {
                if (file.exists(f))
                  read_csv(f) %>% mutate(file = f) else data.frame()
            }))
    }
    if (nrow(data) < 2) {
        warning("You have to provide either data frame or vecor of file names to proceed")
        return(0)
    }
    if (!dir.exists(resultPath))
        dir.create(resultPath, recursive = TRUE)

    if ("totIonCurrent" %in% names(data))
        data <- data %>% mutate(tic = totIonCurrent)
    if ("tic" %in% names(data)) {
        message("TIC adjustment")
        TIClimits <- data %>% group_by(file) %>% summarize(minTIC = min(tic, na.rm = T),
            maxTIC = max(tic, na.rm = T)) %>% summarize(low = max(c(2e+07, minTIC)),
            high = min(maxTIC))
        data <- data %>% filter(tic > TIClimits$low/2 & tic < TIClimits$high * 2)
    } else {
        data$tic <- 10
    }
    mass_tol <- 3e-04

    if (!("file" %in% names(data))) {
        data <- data %>% mutate(file = "file")
    }
    if (is.na(group))
        group = "group"
    if (!(group %in% names(data))) {
        data <- data %>% mutate(group = "group")
    } else if (group != "group") {
        data <- data %>% mutate_(group = sprintf("`%s`", group))
    }
    message("Files to process: ")
    print(unique(data$file))
    message("Groups to process: ")
    print(unique(data$group))
    spectra <- unique(data$seqNum)
    mono_spectra <- unique((data %>% filter(peak == "0"))$seqNum)
    nomono_spectra <- setdiff(spectra, mono_spectra)

    # mass error filtering removing outliers

    ranges <- data %>% group_by(file, peak, ion) %>% summarise(meanerror = mean(masserror),
        medianerror = median(masserror), meansq = mean(masserror^2), gg = median(isoratio/n,
            na.rm = T), ir = median(isoratio, na.rm = T), md = mad(isoratio, na.rm = T),
        nhits = n(), mmin = ir - 3 * md, mmax = ir + 3 * md)
    rranges <- ranges %>% filter(nhits > 10 & gg > 0) %>% summarize(mgg = median(gg,
        na.rm = T), madgg = mad(gg, na.rm = T), mmd = min(md))
    ranges <- ranges %>% ungroup() %>% left_join(rranges, by = c("file", "peak")) %>%
        # mutate(isgood = peak=='0' | ((abs(gg-mgg) < 5*madgg)&(md < 4*mmd)&(md/ir<1)))
    # %>%
    mutate(isgood = peak == "0" | ((md < 8 * mmd) & (md/ir < 1))) %>% select(-mgg,
        -madgg, -mmd)
    data <- data %>% ungroup() %>% filter(abs(masserror) < mass_tol) %>% left_join(ranges,
        by = c("file", "ion", "peak")) %>% filter(!is.na(isgood)) %>% filter(isgood &
        (peak == "0" | (isoratio > mmin & isoratio < mmax)))

    message("Linear-model based method ")
    lm_res_aa <- data %>% filter(n > 0 & I > 10 & isoratio > 0 & tic > 0) %>% mutate(I0 = I/isoratio/100,
        I = I/n, ltic = log10(tic)) %>% group_by(peak, group, file, ion) %>% do({
        mm <- try(MASS::rlm(I ~ I0 + ltic, data = .))
        if (class(mm) != "try-error") {
            tidy(mm)
        } else {
            data.frame()
        }

    })
    lm_res_aa %>% write_csv(file.path(resultPath, sprintf("linearmodel_summary_aa.csv")))
    lm_res <- data %>% filter(n > 0 & I > 10 & tic > 0 & isoratio > 0) %>% mutate(I0 = I/isoratio/100,
        I = I/n, ltic = log10(tic)) %>% group_by(peak, group, file) %>% do({
        mm <- try(MASS::rlm(I ~ I0 + ion + ltic, data = ., method = "MM"))
        if (class(mm) != "try-error") {
            tidy(mm)
        } else {
            data.frame()
        }
    })
    lm_res %>% write_csv(file.path(resultPath, sprintf("linearmodel_summary.csv")))


    if (file.exists(file.path(resultPath, "all_aa.RData"))) {
        message("all_aa.RData file exists, will not recalculate it.")
        load(file.path(resultPath, "all_aa.RData"))
    } else {
        message("Computing 'any' ion: ")
        gpb = dplyr::group_by
        if ("multidplyr" %in% installed.packages()) {
            library(multidplyr)
            cl <- create_cluster()
            cluster_library(cl, "dplyr")
            set_default_cluster(cl)
            gpb <- partition
        }

        any_aa <- data %>% gpb(group, file, seqNum, rt) %>% do({
            dd <- .
            tic_ <- dd$tic[[1]]
            res <- data.frame()
            peaks <- setdiff(unique(dd$peak), "0")
            for (el in peaks) {
                good_aa <- (dd %>% filter(peak == el & I > 0) %>% distinct(ion))$ion
                i0s <- (dd %>% filter(ion %in% good_aa & peak == "0"))$I
                i1s <- (dd %>% filter(ion %in% good_aa & peak == el))$I
                ns <- (dd %>% filter(ion %in% good_aa & peak == el))$n

                rs <- i1s/i0s
                i0s <- i0s[rs < 1]
                i1s <- i1s[rs < 1]
                ns <- ns[rs < 1]
                rs <- rs[rs < 1]
                rs_good <- (abs(rs - median(rs)) < 1.5 * mad(rs))
                r <- sum(i1s[rs_good]/ns[rs_good])/sum(i0s[rs_good])
                if (!is.na(r))
                  res <- res %>% bind_rows(data.frame(peak = el, gamma = r, logI = log2(sum(i1s)),
                    tic = tic_))
            }
            if (nrow(res) > 0) {
                # res$seqNum = dd$seqNum[[1]] res$rt = dd$rt[[1]]
                res$ion = "any"
            }
            res
        }) %>% collect() %>% ungroup()

        message("Puting ions together: ")

        all_aa <- data %>% filter(n > 0 & I > 0) %>% mutate(gamma = isoratio/n, logI = log2(I)) %>%
            select(group, file, peak, gamma, logI, tic, seqNum, rt, ion) %>% bind_rows(any_aa) %>%
            ungroup()
        save(any_aa, all_aa, file = file.path(resultPath, "all_aa.RData"))
    }
    message("Grouping results and writing output tables: ")

    for (pp in unique(all_aa$peak)) {
        vals <- (all_aa %>% filter(peak == pp))$gamma
        vmin <- quantile(vals, 0.05)
        vmax <- quantile(vals, 0.95)
        mmad <- all_aa %>% filter(peak == pp) %>% group_by(ion, file) %>% summarize(m = median(gamma),
            md = mad(gamma))
        aa_data <- all_aa %>% filter(peak == pp) %>% left_join(mmad %>% ungroup(),
            by = c("ion", "file")) %>% filter(md > 0) %>% filter(gamma > (m - 3 *
            md) & gamma < (m + 3 * md))
        aa_data$ion <- factor(aa_data$ion, levels = (mmad %>% summarise(md = median(md)) %>%
            arrange(md))$ion)
        res_f <- aa_data %>% group_by(ion, group, file) %>% summarize(mean = mean(gamma) *
            100, median = median(gamma) * 100, wmean = weighted.mean(gamma, logI) *
            100, sd = sd(gamma) * 100, se = sd(gamma)/sqrt(n()) * 100, n = n(), logI = median(logI)) %>%
            mutate(CV = se/wmean * 100) %>% arrange(ion)
        res_f %>% write_csv(file.path(resultPath, sprintf("file_summary_%s.csv",
            pp)))
        res_f %>% summarize(gmean = mean(mean), gmedian = mean(median), gwmean = mean(wmean),
            sd_mean = sd(mean), sd_median = sd(median), sd_wmean = sd(wmean), n = n()) %>%
            mutate(cv_mean = sd_mean/gmean/sqrt(n) * 100, cv_median = sd_median/gmedian/sqrt(n) *
                100, cv_wmean = sd_wmean/gwmean/sqrt(n) * 100) %>% arrange(cv_wmean) %>%
            write_csv(file.path(resultPath, sprintf("group_summary_%s.csv", pp)))
    }
    if (length(nomono_spectra) > 0) {
        cN <- (data %>% ungroup() %>% filter(seqNum %in% nomono_spectra) %>% filter(peak ==
            "13C") %>% slice(1:1))$n
        gC <- 0.01
        data %>% filter(seqNum %in% nomono_spectra) %>% # filter(abs(masserror)<mass_tol) %>%
        narrow_data <- filter(n > 0 & I > 0) %>% mutate(gamma = isoratio/n * cN *
            gC, logI = log2(I))
        narrowmmad <- narrow_data %>% group_by(peak, file) %>% summarize(m = median(gamma),
            md = mad(gamma))

        narrowres_f <- narrow_data %>% left_join(narrowmmad %>% ungroup(), by = c("peak",
            "file")) %>% filter(md > 0) %>% filter(gamma > (m - 3 * md) & gamma <
            (m + 3 * md)) %>% group_by(ion, peak, group, file) %>% summarize(mean = mean(gamma) *
            100, median = median(gamma) * 100, wmean = weighted.mean(gamma, logI) *
            100, sd = sd(gamma) * 100, se = sd(gamma)/sqrt(n()) * 100, n = n()) %>%
            mutate(CV = se/wmean * 100) %>% arrange(ion)
        narrowres_f %>% write_csv(file.path(resultPath, sprintf("narrow_file_summary.csv")))
        narrowres_f %>% summarize(gmean = mean(mean), gmedian = mean(median), gwmean = mean(wmean),
            sd_mean = sd(mean), sd_median = sd(median), sd_wmean = sd(wmean), n = n()) %>%
            mutate(cv_mean = sd_mean/gmean/sqrt(n) * 100, cv_median = sd_median/gmedian/sqrt(n) *
                100, cv_wmean = sd_wmean/gwmean/sqrt(n) * 100) %>% arrange(cv_wmean) %>%
            write_csv(file.path(resultPath, sprintf("narrow_group_summary.csv")))

    }



    gg <- all_aa %>% distinct(group)
    message("Rendering HTML report: ")
    resultPath <- normalizePath(resultPath)
    rmarkdown::render(system.file("Rmd/isoMS_report.Rmd", package = getPackageName()),
        envir = sys.frame(sys.nframe()), output_file = file.path(resultPath, "isoMS_report.html"))
}
