#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
#library(knitr)
library(survival)
library(survminer)


panel_size_mela_Mb=0.844437
panel_size_pallocca_Mb=0.844437

z<-read.table("ALL.IR-filtered.noCNV.breakmulti.clean.ir_refilter.cosm_id.annovar.hg19_multianno.clean.damaging.concise_tsv.all_tmb.header_added.clinical_data", sep="\t", quote="", header = T, na.strings=c("","NA") )
for(c in colnames(z)){#scalo i tmb
  if(startsWith(c,"tmb")){
    z[,c]<-z[,c]/panel_size_mela_Mb
  }
}
#converto le date
for(c in colnames(z)){
  if(endsWith(c,"_date")){
    z[,c]<-as.Date(z[,c])
  }
}
z$surv_date <- as.Date(ifelse(is.na(as.vector(z$death_date)), as.vector(z$last_followup_date), as.vector(z$death_date)),origin="1970-01-01")
z$surv_time <- z$surv_date - z$therapy1_start_date
z$pfs_time <- z$surv_date - z$therapy1_progression_date
z$death=as.logical(z$death)
z$death_by_melanome=as.logical(z$death_by_melanome)
z$overall_surv<-Surv(as.numeric(z$surv_time),z$death)
z$progfree_surv<-Surv(as.numeric(z$pfs_time),z$death)

z_pallocca<-read.table("pallocca_clinical_data.panel_tmb.header_added", sep="\t", quote="", header = T, na.strings=c("","NA") )
z_pallocca$overall_surv<-Surv(as.numeric(z_pallocca$surv_time),z_pallocca$death)#non ho il dato morto o vivo, assumo tutti morti
z_pallocca$progfree_surv<-Surv(as.numeric(z_pallocca$PFS),z_pallocca$death)#non ho il dato morto o vivo, assumo tutti morti
for(c in colnames(z_pallocca)){#scalo i tmb
  if(startsWith(c,"tmb")){
    z_pallocca[,c]<-z_pallocca[,c]/panel_size_pallocca_Mb
  }
}
z_pallocca$death=as.logical(z_pallocca$death)
z_Snyder<-z_pallocca[z_pallocca$study=="Snyder",]
z_VanAllen<-z_pallocca[z_pallocca$study=="VanAllen",]
z_Hugo<-z_pallocca[z_pallocca$study=="Hugo",]
z_Riaz<-z_pallocca[z_pallocca$study=="Riaz",]

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   #titlePanel("Old Faithful Geyser Data"),
   plotOutput("Plot", height = "800px"),
      # Sidebar with a slider input for number of bins 
   fluidRow(
     column(6,
            sliderInput("mutation_n",
                        "Mutation load (Mb):",
                        min = 1,
                        max = 40,
                        value = 10,
                        step=1,
                        width = "100%",
                        animate = animationOptions(interval = 3000, loop = FALSE, playButton = NULL, pauseButton = NULL)
            )
     ),
     column(2,
            selectInput("vaf","min vaf (%):",
                        c("0.05"="05",
                          "0.07"="07",
                          "0.10"="10",
                          "0.15"="15")
            ),
            selectInput("damaging","damaging:",
                        c("damaging"="damaging",
                          "all"="all")
            )
     ),
     column(2,
            selectInput("therapy1__therapy_type","therapy_type:",
                 c("Immuno"="Immuno",
                   "Target"="Target",
                   "Any"="any"
                 )
            ),
            selectInput("surv_type","survival:",
                        c("Overall"="overall",
                          "Progression free"="progfree"
                        )
            )
     ),
     column(2,
            selectInput("data_source","project:",
                        c("ACC MELA pilot"="acc_mela_pilot",
                          "Pallocca meta"="pallocca",
                          "Snyder"="Snyder",
                          "VanAllen"="VanAllen",
                          "Hugo"="Hugo",
                          "Riaz"="Riaz")
            ) 
     )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  output$Plot <- renderPlot({
    #project selection
    z_sel<-z
    if(input$data_source=="pallocca"){
      z_sel<-z_pallocca
    }
    if(input$data_source=="Snyder"){
      z_sel<-z_Snyder
    }
    if(input$data_source=="VanAllen"){
      z_sel<-z_VanAllen
    }
    if(input$data_source=="Hugo"){
      z_sel<-z_Hugo
    }
    if(input$data_source=="Riaz"){
      z_sel<-z_Riaz
    }
    
    #mutation selection
    tmb_col=paste("tmb",input$vaf,input$damaging,sep="_")
    print(tmb_col)
    z_sel$tmb_high_low<-as.factor(ifelse(z_sel[,tmb_col] >= input$mutation_n,"high","low"))
    
    
    #sample selection
    z_sel2<-z_sel
    if(input$therapy1__therapy_type!="any"){
      z_sel2<-z_sel[z_sel$therapy1__therapy_type==input$therapy1__therapy_type,]
    }
    
    z_sel2$surv_sel<-z_sel2$overall_surv
    if(input$surv_type=="progfree"){
      z_sel2$surv_sel<-z_sel2$progfree_surv
    }
    
    km <- survfit(surv_sel ~ tmb_high_low, data = z_sel2, conf.type = "log-log")
    ggsurvplot(
      km,                      # survfit object with calculated statistics.
      data = z_sel2,                # data used to fit survival curves. 
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      conf.int = TRUE,         # show confidence intervals for 
      # point estimaes of survival curves.
      xlim = c(0,2000),        # present narrower X axis, but not affect
      # survival estimates.
      break.time.by = 365,     # break X axis in time intervals by 500.
      ggtheme = theme_minimal(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T, # colour risk table text annotations.
      risk.table.y.text = FALSE # show bars instead of names in text annotations
      # in legend of risk table
    )
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
