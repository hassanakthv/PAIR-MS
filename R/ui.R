library(markdown)
#b64 <- base64enc::dataURI(file="PAIR-MS logo_wback.png", mime="image/png")
bs64 <- "https://github.com/hassanakthv/PAIR-MS/blob/master/PAIR-MS%20logo_wback.png"
navbarPage(theme = shinytheme("flatly"),
           tags$head(
             tags$style(HTML('.navbar-nav > li > a, .navbar-brand {
                            padding-top:4px !important; 
                            padding-bottom:4px !important;
                            height: 50px;
                            }
                           .navbar {min-height:50x !important;}'))
           ),
           title = img(src=b64,
                         height = '50px',
                         width ='100px',align = "left-top"
                         ),
           fluid = T,
           
           tabPanel("Conversion - mzml to csv file",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("mzmlType", "Molecular Level",
                                     c("Amino Acid Ion"="AA", "Immonimum Ion "="Imm")
                        ),
                        fileInput(inputId = 'mzml_input', label = 'Choose your mzml files', multiple = TRUE),
                        shinyDirButton("dir", "Select Directory for output Converted Files", "Upload"),
                      ),
                      mainPanel(
                        actionButton(inputId = "Convert", label = "Start Conversion"),
                        verbatimTextOutput("converting_text")
                      )
                    )
           ),
           tabPanel("Creating Report with csv file",
                    
                    sidebarLayout(
                      sidebarPanel(
                    
                    textOutput('folder'),
                    shinyDirButton("dir", "Select Result Output Directory", "Upload"),
                    verbatimTextOutput("dir", placeholder = TRUE)  ,
                   
                    fileInput(inputId = 'select_input', label = 'Choose your csv files', accept = ".csv", multiple = TRUE),
                    fileInput(inputId = 'select_de', label = 'Choose your Experiment Design file', accept = ".csv", multiple = FALSE),
                    
                    verbatimTextOutput("lists"),
                 
                    tableOutput("rawtable"),
                    
                    radioButtons(inputId = "MSLevel",
                                       label = "Select class molecules of interest:",
                                       choices = c("Amino Acid ion" = "AA", "Immonium ion" = "Imm"),
                                       selected = "Imm"),
                    # checkboxGroupInput(inputId = "MSLevel",
                    #                    label = "Select class molecules of interest:",
                    #                    choices = names(immoniumIons),
                    #                    selected = names(immoniumIons), inline = T),
                    checkboxGroupInput(inputId = "ShowIon",
                                       label = "Select the Ions of interest:",
                                       choices = names(immoniumIons),
                                       selected = c("P", "V", "L/I", "Hyp"), inline = T), 
                    radioButtons(inputId = "Correct",
                                       label = "Would you like to correct readouts based on the TIC values:",
                                       choices = c("Yes" = T, "No" = F),
                                       selected = T, inline = T),
                    radioButtons(inputId = "Control", 
                                 label = "Do you have specified samples named 'Control'?",
                                 choices = c("Yes" = T, "No" = F), 
                                 selected = F, inline = T),
                    withAnim(),
                    tags$div(id = 'title', h1('Processing:')),
                    actionButton(inputId = "Process", label = "Start Analysis"),
                    
                              
                    
                      
                    
                    ), 
                    mainPanel(
                      #shinycssloaders::withSpinner({
                    #uiOutput('markdown')
                      # fluidRow(
                      #   column(6,box("Amino Acid Composition",plotOutput("AminoAcid"))),
                      #   column(3,box("Ratio vs Log10 TIC",plotlyOutput("gplot")))),
                      # 
                    #### Amino Acid
                    h1("Relative Amino Acid Composition"),
                    plotOutput('AminoAcid_plot', inline = F, width = 1200, height = 800),
                    h1("Amino Acid Table"),
                    DT::dataTableOutput('AminoAcid_table'),
                    
                    #### Carbon
                    h1(HTML(paste0(tags$sup("13"), "C/",tags$sup("12"),"C - Plots"))),
                    h4("Ratio vs Retention Time"),
                    plotOutput('c_plot_rt', inline = F, width = 1200, height = 800),
                    h4("Ratio vs Total Ion Current"),
                    plotOutput('c_plot_tic', inline = F, width = 1200, height = 800),
                    h4("Ratio Distribution for Samples"),
                    plotOutput('c_plot_box', inline = F, width = 1200, height = 800),
                    h4("Ratio Density Plot"),
                    plotOutput('c_plot_density', inline = F, width =1200, height = 800),
                    h4("Mean Ratio vs TIC"),
                    plotOutput('c_plot_mean', inline = F, width = 1200, height = 800),
                    
                    h1(HTML(paste0(tags$sup("13"), "C/",tags$sup("12"),"C - Tables"))),
                    h4("File Summary Table"),
                    DT::dataTableOutput('c_table_sample'),
                    h4("Group Summary Table"),
                    DT::dataTableOutput('c_table_group'),
                    
                    #### Nitrogen
                    h1(HTML(paste0(tags$sup("15"), "N/",tags$sup("14"),"N - Plots"))),
                    h4("Ratio vs Retention Time"),
                    plotOutput('n_plot_rt', inline = F, width = 1200, height = 800),
                    h4("Ratio vs Total Ion Current"),
                    plotOutput('n_plot_tic', inline = F, width = 1200, height = 800),
                    h4("Ratio Distribution for Samples"),
                    plotOutput('n_plot_box', inline = F, width = 1200, height = 800),
                    h4("Ratio Density Plot"),
                    plotOutput('n_plot_density', inline = F, width = 1200, height =800),
                    h4("Mean Ratio vs TIC"),
                    plotOutput('n_plot_mean', inline = F, width = 1200, height = 800),
                    
                    h1(HTML(paste0(tags$sup("15"), "N/",tags$sup("14"),"N - Tables"))),
                    h4("File Summary Table"),
                    DT::dataTableOutput('n_table_sample'),
                    h4("Group Summary Table"),
                    DT::dataTableOutput('n_table_group'),
                    
                    #### Hydrogen
                    h1(HTML(paste0(tags$sup("2"), "H/",tags$sup("1"),"H - Plots"))),
                    h4("Ratio vs Retention Time"),
                    plotOutput('h_plot_rt', inline = F, width = 1200, height = 800),
                    h4("Ratio vs Total Ion Current"),
                    plotOutput('h_plot_tic', inline = F, width = 1200, height = 800),
                    h4("Ratio Distribution for Samples"),
                    plotOutput('h_plot_box', inline = F, width = 1200, height = 800),
                    h4("Ratio Density Plot"),
                    plotOutput('h_plot_density', inline = F, width = 1200, height = 800),
                    h4("Mean Ratio vs TIC"),
                    plotOutput('h_plot_mean', inline = F, width = 1200, height = 800),
                    
                    h1(HTML(paste0(tags$sup("2"), "H/",tags$sup("1"),"H - Tables"))),
                    h4("File Summary Table"),
                    DT::dataTableOutput('h_table_sample'),
                    h4("Group Summary Table"),
                    DT::dataTableOutput('h_table_group'),
                    
                    
                    
                    
                    #DT::dataTableOutput('gplot')
                     #   })
                    )
                    )
           ),
           navbarMenu("More",
                      tabPanel("Table",
                               DT::dataTableOutput("table")
                      ),
                      tabPanel("About",
                               fluidRow(
                                 column(6,
                                        
                                 ),
                                 column(3,
                                        img(class="img-polaroid",
                                            src=paste0("http://upload.wikimedia.org/",
                                                       "wikipedia/commons/9/92/",
                                                       "1919_Ford_Model_T_Highboy_Coupe.jpg")),
                                        tags$small(
                                          "Source: Photographed at the Bay State Antique ",
                                          "Automobile Club's July 10, 2005 show at the ",
                                          "Endicott Estate in Dedham, MA by ",
                                          a(href="http://commons.wikimedia.org/wiki/User:Sfoskett",
                                            "User:Sfoskett")
                                        )
                                 )
                               )
                      )
           )
)
