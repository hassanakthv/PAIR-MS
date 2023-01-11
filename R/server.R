require(shiny)
require(shinydashboard)
require(shinyDirectoryInput)
require(shinyFiles)
require(shinythemes)
require(isoms)
require(shinycssloaders)
require(shinyjs)
require(shinyWidgets)
require(shinyanimate)
require(base64enc)
#require(plotly)

function(input, output, session) {
  
  options(shiny.maxRequestSize=1e6*1024^2)
  volumes <- getVolumes()
 
  
  observe({
    files <- list.files(pattern = ".csv")
    
    updateSelectizeInput(session = session, inputId = 'select_input', choices = files)
  })
  observe({
    files <- list.files(pattern = ".csv")
    
    updateSelectizeInput(session = session, inputId = 'select_de', choices = files)
  })
  

  shinyDirChoose(
    input,
    'dir',
    roots = volumes(),
    filetypes = c('', 'txt', 'bigWig', "tsv", "csv", "bw")
  )
  
  global <- reactiveValues(datapath = getwd())
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    global$datapath
  })
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 if (!"path" %in% names(dir())) return()
                 home <- normalizePath("~")
                 global$datapath <-
                   file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
               })
  

  output$lists <- renderText(paste0("Selected files for analysis:","") , quoted = F)

  toListen <- reactive({
    list(input$select_input,input$select_de,input$Process)
  })
  observeEvent(toListen(),{
    #startAnim(session, 'title', 'bounce')
  data_tab <- data.frame()
  for (i in input$select_input$datapath){
    print(paste("Here is the path", i))
    df <- read.csv(i)
    
    output$ff <- renderPrint(print(head(df)))
    csv_tag <- data.frame("File" = input$select_input$name, "Path" = input$select_input$datapath)
    data_tab <- rbind(data_tab , csv_tag)
    
    
    
    
  }
  data_tab <- data_tab %>% dplyr::distinct()
  output$rawtable <-renderTable(
    data_tab$File, align = 'l', hover = T
  )
  reactive({
    if(input$Convert!=0){
      mzml_names <- c()
      for (i in input$mzml_input$datapath){
        mzml_names <- c(mzml_names, i)
        }
      
      
      output$converting_text = "Conversion started"
      if (input$mzmlType == "AA"){
        mzMLtoCSV_(mzml_names,width = 0.0015)
      }
      if (input$mzmlType == "Imm"){
        mzMLtoCSV(mzml_names,width = 0.0015)
      }
      
    }
    
    
  })
  
  #reactive({
  if(input$Process!=0){
    runjs(paste0('$("#Process").css("animation","")'))
    startAnim(session, 'title', 'bounce')
    
    withProgress(message = "Processing Data",{
      dd <- c("P", "Hyp", "L/I", "R", "R2", "K", "Hyl", "A", "H", "F", "S","V", "Y", "N","N2","D","Q","E","T","Pyr","m57","K2","m89",
          "m73","m87","M","C","camC","xC","W")
      v <- PAIRMS(file = input$select_de$datapath, data_tab = data_tab, outdir = global$datapath, ioi = dd, ioi_ = input$ShowIon, correct = input$Correct)
    })
    
    withProgress(message = "Making Plots",{
    ##### 
    ### Carbon
    output$c_plot_rt <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="13C" & ion %in% input$ShowIon),aes(x = rt/60, y = gamma*100, fill = group, text = paste0("RT: ",round(rt/60,2),'\n', "Ratio: ", round(gamma*100, 4), 
                                                                                                            '\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_point(size = 3 , shape = 21, color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("13"), "C/",tags$sup("12"), "C, %")), x = "Retention Time, min")
        labs(y = bquote(''^13*'C'/''^12*'C, (%)'), x = "Retention Time, min")
     #g1 %>% ggplotly(tooltip = c("text")) -> g2
     #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$c_plot_tic <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="13C" & ion %in% input$ShowIon),aes(x = log10(tic), y = gamma*100, fill = group, text = paste0("Log10(TIC): ",round(log10(tic),2),'\n', "Ratio: ", round(gamma*100, 4), 
                                                                                                              '\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_point(size = 3 , shape = 21, color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("13"), "C/",tags$sup("12"), "C, %")), x = HTML(paste0("TIC, Log", tags$sub("10"))))
        labs(y = bquote(''^13*'C'/''^12*'C, (%)'), x = bquote(TIC ~ (log[10])))
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$c_plot_box <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="13C" & ion %in% input$ShowIon),aes(x = file, y = gamma*100, fill = group, text = paste0("SampleID: ",file,'\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        theme(axis.text.x = element_blank()) + 
        geom_boxplot(color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("13"), "C/",tags$sup("12"), "C, %")), x = "Sample ID")
        labs(y = bquote(''^13*'C'/''^12*'C, (%)'), x = "Sample ID")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      ##subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$c_plot_density <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="13C" & ion %in% input$ShowIon),aes(x = gamma*100, fill = group, text = paste0("SampleID: ",file,'\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_density(color = "black", alpha = 0.6) + #labs(x = HTML(paste0(tags$sup("13"), "C/",tags$sup("12"), "C, %")), y = "Density")
        labs(x = bquote(''^13*'C'/''^12*'C, (%)'), y = "Density")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      ##subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    }) 
    
    
    output$c_plot_mean <- renderPlot({
      
      v %>% 
        dplyr::filter(peak == "13C" & ion %in% input$ShowIon) %>%
        mutate(logI = log10(I)) %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*100, 
          median=median(gamma)*100, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*100, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic),na.rm=T)) %>% 
        arrange(tic) %>% 
        ggplot(aes(x=tic, y=mean,  color=group)) +
        geom_point(size=3) + geom_errorbar(aes(ymin=mean-3*se, ymax=mean+3*se),size=2)  +
        geom_smooth(method = "lm", se=F) + 
        theme_bw() + facet_grid(.~ion, scales='free') + 
        labs(y = bquote('Mean'^13*'C'/''^12*'C, (%)'), x = bquote(TIC ~ (log[10]))) -> g_tic
        #labs(y = HTML(paste0('Mean', tags$sup("13"), "C/",tags$sup("12"), "C, %")), x = HTML(paste0("TIC, Log", tags$sub("10")))) -> g_tic
      
      
      #g_tic %>% ggplotly(tooltip = c("text"))
      
      g_tic
    })
    
    output$c_table_sample <- renderDT(server = F, {
      
      v %>% 
        dplyr::filter(peak=="13C" & ion != 'any') %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*100, 
          median=median(gamma)*100, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*100, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic), na.rm=T)) %>% 
        arrange(ion) -> res_group
      
      res_group %>% 
        ungroup() %>% 
        arrange(ion, group, file) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
          caption = htmltools::tags$caption("Sample Summary")
        )) %>%
        formatSignif(c('mean','median','sd'), 6) %>%  formatSignif(c('se','tic'), 2) %>%
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
          backgroundPosition = 'center')
        # ) %>%
        # formatStyle(
        #   'wmean',
        #   background = styleColorBar(res_group$wmean, 'azure'),
        #   backgroundSize = '100% 90%',
        #   backgroundRepeat = 'no-repeat',
        #   backgroundPosition = 'center'
        # )
        # 
      
      
        
    })
    
    output$c_table_group <- renderDT(server = F, {
      v %>% 
        dplyr::filter(peak=="13C" & ion != 'any') %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*100, 
          median=median(gamma)*100, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*100, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic), na.rm=T)) %>% 
        arrange(ion) %>% 
        ungroup() %>% 
        group_by(ion, group) %>% 
        summarize(
          gg_mean = mean(mean, na.rm = T), 
          gg_med = mean(median, na.rm = T), 
          sd_mean = sd(mean, na.rm = T), 
          sd_med = sd(median, na.rm = T), 
          cv_mean = sd(mean, na.rm = T)/sqrt(n())*100, 
          cv_med = sd(median, na.rm = T)/sqrt(n())*100, 
          tic = mean(tic, na.rm = T),
          n_rep = n(), 
          n_total = sum(n, na.rm = T)) -> res_gg
      res_gg %>% ungroup() %>% 
        arrange( desc(n_total),ion, group,) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
          caption = htmltools::tags$caption("Group Summary")
        )) %>%
        formatSignif(c('gg_mean','gg_med','sd_mean', "sd_med"), 6) %>%  formatSignif(c('cv_mean','cv_med','tic'), 2) %>%
        formatStyle(
          'gg_mean',
          background = styleColorBar(res_gg$gg_mean, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        ) %>%
        formatStyle(
          'gg_med',
          background = styleColorBar(res_gg$gg_med, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center')
      
        
    })
    
    ##### 
    ### Nitrogen
    
    output$n_plot_rt <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="15N" & ion %in% input$ShowIon),aes(x = rt/60, y = gamma*100, fill = group, text = paste0("RT: ",round(rt/60,2),'\n', "Ratio: ", round(gamma*100, 4), 
                                                                                                                                       '\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_point(size = 3 , shape = 21, color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("15"), "N/",tags$sup("14"), "N, %")), x = "Retention Time, min")
        labs(y = bquote(''^15*'N'/''^14*'N, (%)'), x = "Retention Time, min")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$n_plot_tic <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="15N" & ion %in% input$ShowIon),aes(x = log10(tic), y = gamma*100, fill = group, text = paste0("Log10(TIC): ",round(log10(tic),2),'\n', "Ratio: ", round(gamma*100, 4), 
                                                                                                                                            '\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_point(size = 3 , shape = 21, color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("15"), "N/",tags$sup("14"), "N, %")), x = HTML(paste0("TIC, Log", tags$sub("10"))))
        labs(y = bquote(''^15*'N'/''^14*'N, (%)'), x = bquote(TIC~(log[10])))
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$n_plot_box <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="15N" & ion %in% input$ShowIon),aes(x = file, y = gamma*100, fill = group, text = paste0("SampleID: ",file,'\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        theme(axis.text.x = element_blank()) + 
        geom_boxplot(color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("15"), "N/",tags$sup("14"), "N, %")), x = "Sample ID")
        labs(y = bquote(''^15*'N'/''^14*'N, (%)'), x = "Sample ID")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$n_plot_density <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="15N" & ion %in% input$ShowIon),aes(x = gamma*100, fill = group, text = paste0("SampleID: ",file,'\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_density(color = "black", alpha = 0.6) + #labs(x = HTML(paste0(tags$sup("15"), "N/",tags$sup("14"), "N, %")), y = "Density")
        labs(x = bquote(''^15*'N'/''^14*'N, (%)'), y = "Density")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    }) 
    
    
    output$n_plot_mean <- renderPlot({
      
      v %>% 
        filter(peak == "15N" & ion %in% input$ShowIon) %>%
        mutate(logI = log10(I)) %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*100, 
          median=median(gamma)*100, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*100, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic),na.rm=T)) %>% 
        arrange(tic) %>% 
        ggplot(aes(x=tic, y=mean,  color=group)) +
        geom_point(size=3) + geom_errorbar(aes(ymin=mean-3*se, ymax=mean+3*se),size=2)  +
        #geom_smooth(method = "lm", se=F) + 
        theme_bw() + facet_grid(.~ion, scales='free') + 
        #labs(y = HTML(paste0('Mean', tags$sup("15"), "N/",tags$sup("14"), "N, %")), x = HTML(paste0("TIC, Log", tags$sub("10")))) -> g_tic
        labs(y = bquote('Mean '^15*'N'/''^14*'N, (%)'), x = bquote(TIC~(log[10]))) -> g_tic
        g_tic
      #g_tic %>% ggplotly(tooltip = c("text"))
      
      
    })
    
    output$n_table_sample <- renderDT({
      v %>% 
        dplyr::filter(peak=="15N" & ion != 'any') %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*100, 
          median=median(gamma)*100, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*100, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic), na.rm=T)) %>% 
        arrange(ion) -> res_group
      
      res_group %>% 
        ungroup() %>% 
        arrange(ion, group, file) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
          caption = htmltools::tags$caption("Sample Summary")
        )) %>%
        formatSignif(c('mean','median','sd'), 6) %>%  formatSignif(c('se','tic'), 2) %>%
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
          backgroundPosition = 'center')
      
    })
    
    
    output$n_table_group <- renderDT({
      v %>% 
        dplyr::filter(peak=="15N" & ion != 'any') %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*100, 
          median=median(gamma)*100, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*100, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic), na.rm=T)) %>% 
        arrange(ion) %>% 
        ungroup() %>% 
        group_by(ion, group) %>% 
        summarize(
          gg_mean = mean(mean, na.rm = T), 
          gg_med = mean(median, na.rm = T), 
          sd_mean = sd(mean, na.rm = T), 
          sd_med = sd(median, na.rm = T), 
          cv_mean = sd(mean, na.rm = T)/sqrt(n())*100, 
          cv_med = sd(median, na.rm = T)/sqrt(n())*100, 
          tic = mean(tic, na.rm = T),
          n_rep = n(), 
          n_total = sum(n, na.rm = T)) -> res_gg
      res_gg %>% ungroup() %>% 
        arrange( desc(n_total),ion, group,) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
          caption = htmltools::tags$caption("Group Summary")
        )) %>%
        formatSignif(c('gg_mean','gg_med','sd_mean', "sd_med"), 6) %>%  formatSignif(c('cv_mean','cv_med','tic'), 2) %>%
        formatStyle(
          'gg_mean',
          background = styleColorBar(res_gg$gg_mean, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        ) %>%
        formatStyle(
          'gg_med',
          background = styleColorBar(res_gg$gg_med, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center')
      
      
      
      
    })
    
    #####
    ### Hydrogen
    
    output$h_plot_rt <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="2H" & ion %in% input$ShowIon),aes(x = rt/60, y = gamma*1e6, fill = group, text = paste0("RT: ",round(rt/60,2),'\n', "Ratio: ", round(gamma*1e6, 1), 
                                                                                                                                       '\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_point(size = 3 , shape = 21, color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("2"), "H/",tags$sup("1"), "H, ppm")), x = "Retention Time, min")
        labs(y =  bquote(''^2*'H'/''^1*'H, (ppm)'), x = "Retention Time, min")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$h_plot_tic <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="2H" & ion %in% input$ShowIon),aes(x = log10(tic), y = gamma*1e6, fill = group, text = paste0("Log10(TIC): ",round(log10(tic),2),'\n', "Ratio: ", round(gamma*1e6, 1), 
                                                                                                                                            '\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_point(size = 3 , shape = 21, color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("2"), "H/",tags$sup("1"), "H, ppm")), x = HTML(paste0("TIC, Log", tags$sub("10"))))
        labs(y = bquote(''^2*'H'/''^1*'H, (ppm)'),x= bquote(TIC~(log[10])))
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$h_plot_box <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="2H" & ion %in% input$ShowIon),aes(x = file, y = gamma*1e6, fill = group, text = paste0("SampleID: ",file,'\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        theme(axis.text.x = element_blank()) + 
        geom_boxplot(color = "black", alpha = 0.6) + #labs(y = HTML(paste0(tags$sup("2"), "H/",tags$sup("1"), "H, ppm")), x = "Sample ID")
        labs(y = bquote(''^2*'H'/''^1*'H, (ppm)'), x = "Sample ID")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    })
    
    output$h_plot_density <- renderPlot({
      
      g1 <- ggplot(v %>% dplyr::filter(peak=="2H" & ion %in% input$ShowIon),aes(x = gamma*1e6, fill = group, text = paste0("SampleID: ",file,'\n', "Group: ", group))) +
        theme_bw() + facet_wrap(ion~., scales = "free") +
        theme(strip.placement = "inside") +
        geom_density(color = "black", alpha = 0.6) + #labs(x = HTML(paste0(tags$sup("2"), "H/",tags$sup("1"), "H, ppm")), y = "Density")
        labs(x = bquote(''^2*'H'/''^1*'H, (ppm)'), y ="Density")
      #g1 %>% ggplotly(tooltip = c("text")) -> g2
      #subplot(list(g2, g2), nrows = 2, margin = .1)
      #g2
      g1
    }) 
    
    
    output$h_plot_mean <- renderPlot({
      
      v %>% 
        filter(peak == "2H" & ion %in% input$ShowIon) %>%
        mutate(logI = log10(I)) %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*1e6, 
          median=median(gamma)*1e6, 
          #wmean = weighted.mean(gamma, logI)*1e6, 
          sd = sd(gamma)*1e6, 
          se = sd(gamma)/sqrt(n())*100, 
          n=n(), 
          tic=median(log10(tic),na.rm=T)) %>% 
        arrange(tic) %>% 
        ggplot(aes(x=tic, y=mean,  color=group)) +
        geom_point(size=3) + geom_errorbar(aes(ymin=mean-3*se, ymax=mean+3*se),size=2)  +
        geom_smooth(method = "lm", se=F) + 
        theme_bw() + facet_grid(.~ion, scales='free_y') + 
        #labs(y = HTML(paste0('Mean', tags$sup("2"), "H/",tags$sup("1"), "H, ppm")), x = HTML(paste0("TIC, Log", tags$sub("10")))) -> g_tic
        labs(y = bquote('Mean '^2*'H'/''^1*'H, (ppm)'), x = bquote(TIC~(log[10]))) -> g_tic
      g_tic
      #g_tic %>% ggplotly(tooltip = c("text"))
      
      
    })
    
    output$h_table_sample <- renderDT({
      v %>% 
        dplyr::filter(peak=="2H" & ion != 'any') %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*1e6, 
          median=median(gamma)*1e6, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*1e6, 
          se = sd(gamma)/sqrt(n())*1e6, 
          n=n(), 
          tic=median(log10(tic), na.rm=T)) %>% 
        arrange(ion) -> res_group
      
      res_group %>% 
        ungroup() %>% 
        arrange(ion, group, file) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
          caption = htmltools::tags$caption("Sample Summary")
        )) %>%
        formatSignif(c('mean','median','sd'), 2) %>%  formatSignif(c('se','tic'), 2) %>%
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
          backgroundPosition = 'center')
    })
    
    output$h_table_group <- renderDT({
      v %>% 
        dplyr::filter(peak=="2H" & ion != 'any') %>% 
        group_by(ion, group,file) %>% 
        summarize(
          mean=mean(gamma)*1e6, 
          median=median(gamma)*1e6, 
          #wmean = weighted.mean(gamma, logI)*100, 
          sd = sd(gamma)*1e6, 
          se = sd(gamma)*1e6/sqrt(n()), 
          n=n(), 
          tic=median(log10(tic), na.rm=T)) %>% 
        arrange(ion) %>% 
        ungroup() %>% 
        group_by(ion, group) %>% 
        summarize(
          gg_mean = mean(mean, na.rm = T), 
          gg_med = mean(median, na.rm = T), 
          sd_mean = sd(mean, na.rm = T), 
          sd_med = sd(median, na.rm = T), 
          cv_mean = sd(mean, na.rm = T)/sqrt(n())*100, 
          cv_med = sd(median, na.rm = T)/sqrt(n())*100, 
          tic = mean(tic, na.rm = T),
          n_rep = n(), 
          n_total = sum(n, na.rm = T)) -> res_gg
      res_gg %>% ungroup() %>% 
        arrange( desc(n_total),ion, group,) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
          caption = htmltools::tags$caption("Group Summary")
        )) %>%
        formatSignif(c('gg_mean','gg_med','sd_mean', "sd_med"), 2) %>%  formatSignif(c('cv_mean','cv_med','tic'), 2) %>%
        formatStyle(
          'gg_mean',
          background = styleColorBar(res_gg$gg_mean, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        ) %>%
        formatStyle(
          'gg_med',
          background = styleColorBar(res_gg$gg_med, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center')
      
    })
    
    
    #####
    output$AminoAcid_plot <- renderPlot({
      
      v %>% 
        filter(peak=="13C" & ion != 'any') %>% 
        group_by(group, ion) %>% 
        summarize(abundance=sum(I0, na.rm=T)) %>% 
        mutate(abundance = abundance/sum(abundance)*100) -> aa_abundances
      aa_abundances$ion <- factor(aa_abundances$ion, 
                                  levels=aa_abundances %>% 
                                    group_by(ion) %>% 
                                    summarize(s=sum(abundance)) %>% 
                                    arrange(desc(s)) %>% 
                                    select(ion) %>% unlist())
      aa_abundances %>% 
        ggplot(aes(x=ion, y=abundance, fill=group)) + geom_col(position='dodge') + 
        theme_bw() + ylab("Abundance, %") + xlab("Amino Acid") + ggtitle("Amino Acid Abundance") + 
        theme(plot.title = element_text(size = 10, face = "bold", color = "black")) -> AA_plot
      #ggplotly(AA_plot)
      AA_plot
      
      
    })
    
    output$AminoAcid_table <- renderDT(server = F,{
      v %>% 
        filter(peak=="13C" & ion != 'any') %>% 
        group_by(group, ion) %>% 
        summarize(abundance=sum(I0, na.rm=T)) %>% 
        mutate(abundance = abundance/sum(abundance)*100) -> aa_abundances
      
      aa_abundances %>% 
        arrange(group, desc(abundance),ion) %>% 
        datatable(extensions = 'Buttons', options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
        ), caption = htmltools::tags$caption("Amino Acid Composition")) %>% 
        formatRound("abundance", 2) %>%
        formatStyle(
          'abundance',
          background = styleColorBar(aa_abundances$abundance, 'azure'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
      
    })
    
    })
  }
  })
 
  
  

  output$conversion <- renderText({
    if (input$mzmlType=="AA"){
      "mzMLtoCSV_() Rendering to Create amino acid files"
    }
    if (input$mzmlType=="Imm"){
      "mzMLtoCSV() Rendering to Create immonium files"
    }
  })
  
  
  
}