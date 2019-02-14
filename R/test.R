if(F){
  library(isoms)
  library(parallel)
  setwd("/Volumes/Seagate 2T/Immonium/Fusion/HeLa/Conc/")
  ff <- list.files(".", pattern="^201704.*.mzML$")
  ncpu <- 8
  message(sprintf("Processing %d files. You can go and have some beer meanwhile", length(ff)))
  result <- bind_rows(mclapply(
    ff, function(ii){
      outfile <- sub(".mzML","_gaussian.csv", ii)
      if(!file.exists(outfile)){
        xxx <- analyze_immoniums(ii, ions=immoniumIons, width=0.0015, fixSigma=T)
        write_csv(xxx, outfile)
      }else
        #message(sprintf("File already processed: %s", outtile))
        xxx <- suppressMessages(read_csv(outfile))
      xxx %>% mutate(
        group=sub(".*_(\\d+ng).*","\\1",ii),
        file=sub("^.*?(\\d+_F_(\\d+)_.*).mzML","\\1",ii))
    }, mc.cores=ncpu, mc.preschedule=TRUE))

  result <- read_csv("/Volumes/Seagate 2T/Immonium/Fusion/Peptome/20170727/20170727_F_02_PeptoneCaseinTryptic_03_Imm_gaussian.csv") %>%
    mutate(file='f',group='g')

  result %>%
    filter(n>0 & I>0 & rt>600 ) %>%
  #  filter(ion %in% c("P","L","V")) %>%
    mutate(I0 = I/isoratio/100, I02 = I0*I0, I03 = sqrt(I0), I=I/n) %>%
    filter(I0>0) %>%
    group_by(peak, group, file) %>%
    do({
      mm <- try(MASS::rlm(I~I0 + ion + I02+0, data=., method='MM'))
      if(class(mm)!='try-error'){
        tidy(mm)
      }else{
        data.frame()
      }
    }) -> lm_res

  lm_res %>%
    filter(term=='I0') %>%
    ggplot(aes(x=file, y=estimate, color=group)) + geom_point(size=3) + facet_grid(peak~., scales = 'free') + theme_bw()+theme(legend.position = "",axis.text.x = element_text(angle = 90, hjust = 1))

  result %>%
    filter(n>0 & I>0 & ion=='P' & file=='05' & peak=='13C') %>%
    mutate(I0 = I/isoratio/100, I=I/n) %>% filter(rt>600) %>%
    MASS::rlm("I~I0", data=.) -> mm

  result %>%
    filter(n>0 & I>0 & ion=='P' & peak=='13C' & group %in% c('200ng','400ng','800ng')) %>%
    mutate(I0 = (I/isoratio/100), I=(I/n)) %>% filter(rt>600&I<1.5*I0) %>%
    ggplot(aes(x=sqrt(I0), y=I/I0, color=group) ) + geom_point(size=0.1, alpha=0.8) + theme_minimal() #+ scale_x_log10()


  result %>%
    filter(n>0 & I>0 & ion=='P' & file=='06' & peak=='13C') %>%
    mutate(I0 = I/isoratio/100, I02 = I0*I0, I=I/n) %>%
    ggplot(aes(x=I0, y=I)) + geom_point() + scale_x_log10() + scale_y_log10()


  lm_res %>%
    select(peak, group, file, term, estimate) %>%
    spread(term, estimate) %>%
    ggplot(aes(x=`ionP`, y=I0, color=group)) + geom_point() + facet_grid(peak~., scales = 'free')
  immoniumIons = list(
     'P' = c('C' = 4, 'H' = 8, 'N' = 1), # 70.06587
      'L' = c('C' = 5, 'H' = 12, 'N' = 1), # 86.09707
       'R' = c('C' = 5, 'H' = 13, 'N' = 4),  # 129.11358
     #  'K' = c('C' = 5, 'H' = 13, 'N' = 2),  # 101.10784
       'H' = c('C'= 5, 'H' = 8, 'N' = 3),   # 110.07174
     'F' = c('C'= 8, 'H' = 10, 'N'= 1),   # 120.08115
     #st(mclapply(1:32, function(x) sum(rnorm(1e7))))  'S' = c('C'= 2, 'H' = 6, 'N'= 1, 'O' = 1),  # 60.04514
       'V' = c('C' = 4, 'H' = 10, 'N' = 1), # 72.08151
       'Y' = c('C' = 8, 'H' = 10, 'N' = 1, 'O' = 1),    # 136.07593
     #  'N' = c('C' = 3, 'H' = 7, 'N' = 2, 'O' = 1),
     #  'D' = c('C' = 3, 'H' = 6, 'N' = 1, 'O' = 2),
     #  'Q' = c('C' = 4, 'H' = 9, 'N' = 2, 'O' = 1),
       'E' = c('C' = 4, 'H' = 8, 'N' = 1, 'O' = 2),
       'T' = c('C' = 3, 'H' = 8, 'N' = 1, 'O' = 1),
       'Pyr' = c('C' = 4, 'H' = 6, 'N' = 1, 'O' = 1),
     #  'm57' = c('C' = 4, 'H' = 9),
       'm84' = c('C' = 5, 'H' = 10, 'N' = 1)
     #  'm89' = c('C' = 4, 'H' = 9, 'O' = 2),
     #  'm73' = c('C' = 3, 'H' = 5, 'O' = 2),
     #  'm87' = c('C' = 4, 'H' = 7, 'O' = 2)
       )

  ii <- '/Volumes/Seagate 2T/Immonium/Fusion/Peptome/20170727/fullP/20170727_F_03_PeptoneCaseinTryptic_Imm_400.mzML'
  outfile <- sub(".mzML","_gaussian2.csv", ii)
  if(!file.exists(outfile)){
    xxx <- analyze_immoniums2(ii, ions=immoniumIons, width=0.0015, fixSigma=T)
    write_csv(xxx, outfile)
    xxx %>% filter(R2>0.95 & rt<4000) %>% filter(element=='H' & ion %in%c('P','L','V','H')) %>%
      mutate(gamma = estimate/n) %>%
      ggplot(aes(x=rt, y=gamma, color=ion)) + geom_point() + geom_smooth() + scale_y_log10()
    xxx %>% filter(R2>0.95 & rt<4000 & element != 'DM') %>%
      mutate(gamma = estimate/n) %>%
      group_by(ion, element) %>%
      summarize(mmedian = median(gamma), mmean=mean(gamma), nn=n(), se=sd(gamma)/sqrt(nn+1))
  }

  data %>%
    group_by(file, element, ion) %>%
    summarize(gg = median(estimate/n, na.rm=T),ir = median(estimate, na.rm=T), md = mad(estimate, na.rm=T), nhits=n(), mmin = ir-3*md, mmax=ir+3*md) -> ranges
  ranges %>%
    filter(nhits>10 & gg>0) %>%
    summarize(mgg = median(gg, na.rm=T), madgg =mad(gg, na.rm=T), mmd = min(md)) -> rranges
  ranges <- ranges %>%
    ungroup() %>%
    left_join(rranges, by=c('file','element')) %>%
    mutate(isgood = ((abs(gg-mgg) < 5*madgg)&(md < 4*mmd)&(md/ir<1))) %>%
    select(-mgg, -madgg, -mmd)
  data <- data %>%
    ungroup()%>%
    left_join(ranges, by=c('file','ion','element')) %>%
    filter(!is.na(isgood)) %>%
    filter(isgood & (estimate > mmin & estimate < mmax))

  message("Linear-model based method ")
  data %>%
    filter(n>0 & I>10) %>%
    mutate(I0 = I/100, I=I0*estimate/n*100) %>%
    group_by(element, group, file, ion) %>%
    do({
      mm <- try(MASS::rlm(I~I0, data=.))
      if(class(mm)!='try-error'){
        tidy(mm)
      }else{
        data.frame()
      }

    })  %>%
    mutate(peak=element) %>%
    mutate(peak=sub("C","13C",peak)) %>%
    mutate(peak=sub("N","15N",peak)) %>%
    mutate(peak=sub("H","2H",peak))-> lm_res_aa
  lm_res_aa %>%write_csv(file.path(resultPath, sprintf("linearmodel2_summary_aa.csv")))
  data %>%
    filter(n>0 & I>10) %>%
    mutate(I0 = I/100, I=I0*estimate/n*100) %>%
    group_by(element, group, file) %>%
    do({
      mm <- try(MASS::rlm(I~I0+ion, data=., method='MM'))
      if(class(mm)!='try-error'){
        tidy(mm)
      }else{
        data.frame()
      }
    })  %>%
    mutate(peak=element) %>%
    mutate(peak=sub("C","13C",peak)) %>%
    mutate(peak=sub("N","15N",peak)) %>%
    mutate(peak=sub("H","2H",peak))-> lm_res
  lm_res %>%write_csv(file.path(resultPath, sprintf("linearmodel2_summary.csv")))

  data %>%
    filter(n>0 & I>10 & rt>600 & rt<3600) %>%
    mutate(I0 = I/100, I=I0*estimate/n*100) %>%
    group_by(element) %>%
    do({
      mm <- try(MASS::rlm(I~I0+ion+log10(group), data=., method='MM'))
      if(class(mm)!='try-error'){
        tidy(mm)
      }else{
        data.frame()
      }
    })  %>%
    mutate(peak=element) %>%
    mutate(peak=sub("C","13C",peak)) %>%
    mutate(peak=sub("N","15N",peak)) %>%
    mutate(peak=sub("H","2H",peak))-> lm_res

}
