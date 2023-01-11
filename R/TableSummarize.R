TableSummarize <- function(data = df, elem = NULL, AA_Comp = F, grp =F){
  #aa_data <- result %>% dplyr::filter(peak==elem)
  aa_data <- df %>% dplyr::filter(peak==elem)
  if(AA_Comp==F){
  aa_data %>% 
      mutate(logI = log10(I)) %>% 
    dplyr::group_by(ion, group,file) %>% 
    dplyr::summarise(
      mean=mean(gamma, na.rm = T)*100, 
      median=median(gamma, na.rm = T)*100, 
      #wmean = weighted.mean(gamma, logI, na.rm = T)*100, 
      sd = sd(gamma, na.rm = T)*100, 
      se = sd(gamma, na.rm = T)/sqrt(n())*100, 
      n=n(), 
      tic=median(log10(tic), na.rm=T)) %>% 
    arrange(ion) -> res_group
  res_group %>% 
    ungroup() %>% 
    arrange(ion, group, file) %>% 
    datatable(extensions = 'Buttons', options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    )) %>%
    formatSignif(c('mean','median','wmean','sd'), 6) %>%  formatSignif(c('se','tic'), 2) %>%
    formatStyle(
      'mean',
      background = styleColorBar(res_group$mean, 'azure'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    ) %>%
    formatStyle(
      'median',
      background = styleColorBar(res_group$median, 'azure'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center') -> res_group
    # ) %>%
    # formatStyle(
    #   'wmean',
    #   background = styleColorBar(res_group$wmean, 'azure'),
    #   backgroundSize = '100% 90%',
    #   backgroundRepeat = 'no-repeat',
    #   backgroundPosition = 'center'
    # ) -> res_group
  res_group -> fileres_C
  if (grp == T){
    aa_data %>% 
      mutate(logI = log10(I)) %>% 
      dplyr::group_by(ion, group,file) %>% 
      dplyr::summarise(
        mean=mean(gamma, na.rm = T)*100, 
        median=median(gamma, na.rm = T)*100, 
        #wmean = weighted.mean(gamma, logI, na.rm = T)*100, 
        sd = sd(gamma, na.rm = T)*100, 
        se = sd(gamma, na.rm = T)/sqrt(n())*100, 
        n=n(), 
        tic=median(log10(tic), na.rm=T)) %>% 
      arrange(ion) %>% 
      summarize(gmean=mean(mean), gmedian=mean(median), sd_mean=sd(mean), sd_median=sd(median),  n=n()) %>%
      mutate(
        cv_mean = sd_mean/gmean/sqrt(n)*100,
        cv_median = sd_median/gmedian/sqrt(n)*100,
      ) -> res
    res %>% 
      datatable(extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      )) %>%  
      formatSignif(c('gmean','gmedian','sd_mean','sd_median'), 6) %>%
      formatSignif(c('cv_mean','cv_median'), 2) %>%
      formatStyle(
        'gmean',
        background = styleColorBar(res$gmean, 'azure'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      formatStyle(
        'gmedian',
        background = styleColorBar(res$gmedian, 'azure'),
        backgroundSize = '100% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      ) -> res
    return(res)
    
    
  }
  if (grp ==F){
  return(fileres_C)
  }
  }
 
  
  
}


Table_Formatter <- function(dt_tab){
  dt_tab %>%
    formatSignif(c('mean','median','wmean','sd'), 6) %>%  formatSignif(c('se','tic'), 2) %>%
    formatStyle(
      'mean',
      background = styleColorBar(res_group$mean, 'azure'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    ) %>%
    formatStyle(
      'median',
      background = styleColorBar(res_group$median, 'azure'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    ) %>%
    formatStyle(
      'wmean',
      background = styleColorBar(res_group$wmean, 'azure'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    )
  
}