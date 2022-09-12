library(shiny)
library(shinydashboard)
library(rintrojs)
library(shinyhelper)
library(magrittr)
library(nnet)
library(MASS)
library(genalg)
library(survminer)
library(survival)
library(randomForestSRC)
library(party)

load("allmodel.RData")
header <- dashboardHeader(title = "3S Score")

sidebar <- dashboardSidebar(
  width=300,
  introjsUI(),
  introBox(
    fluidRow(
      column(12, 
             align = "middle",
             flowLayout(actionButton("button1", "Click for instructions")
             )
      )
      
    ),
    #actionButton("button1", "Click for instructions"),
    data.step = 1,
    data.intro = "Input these 2 clinical information below.",
    fluidRow(
      column(12, 
             align = "middle",
             flowLayout(helpText(h4("  Please enter the patient's clinical information"))
             )
      )
      
    ),
    selectInput("grade","grade", 
                choices = list("I or II" = "oneortwo", "III" = "three","IV" = "four"), selected = 1),
    selectInput("primary","primary", 
                choices = list("yes" = "yes", "no" = "no"), selected = 1)
  ),
  fluidRow(
    column(12, 
           align = "middle",
           flowLayout(helpText(h4("Please upload expression - data below."))%>%
                        helper(icon = "question-circle",
                               colour = "white",
                               type = "markdown",
                               content = "mymarkdown")
           )
    )
    
  ),
  introBox(
    data.step = 2,
    data.intro = "Please input the 69 gene expression data of your patients.",
    fileInput("file1", "Upload CSV File",
              multiple = TRUE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv"))
        
      
    
  ),
  introBox(
    data.step = 3,
    data.intro = "You can click here to see our example data.",
    fluidRow(
      column(12, 
             align = "middle",
             flowLayout(downloadButton('downloadData', 'Click to download example data', class = "butt"),
                        tags$head(tags$style(".butt{background-color:#000000;}"))
             )
      )

    )
    #downloadButton('downloadData', 'Click to download sampledata', class = "butt"),
    #tags$style(type='text/css', "#run_report { width:50%; margin-left: 100px;}"),
    #tags$head(tags$style(".butt{background-color:#000000;}"))
  )
  
  
  
)

body <- dashboardBody(
  HTML('<meta name="viewport" content="width=1024">'),
  
  sidebarLayout(
    sidebarPanel(
      width = 6,
      hr(),
      h3("Below are the results about 3S score:"),
      hr(),
      introBox(
        data.step = 4,
        data.intro = "You can get points of 2 clinical information here.",
      
        h3("The points of Grade:"),
        textOutput("point_grade"),
        
        
        
        h3("The points of Primary:"),
        textOutput("point_primary")
      ),
      
      introBox(
        data.step = 5,
        data.intro = "You can get 3S score related result here.",
        h3("The Risk score value:"),
        textOutput("riskscore"),
        h3("The points of Risk score:"),
        textOutput("pointsss_risk")
        
      ),
      introBox(
        data.step = 6,
        data.intro = "You can get total points and predicted mortality of this patient.",
        h3("Total points:"),
        textOutput("point_total"),
        
        h3("Mortality of"),
        textOutput("rate_year"),
        tags$head(tags$style("#rate_year{color: red;
                                 font-size: 20px;
                                 font-style: italic;
                                 }"
        )
        )
        
      )
    ),
    mainPanel(
      width = 6,
      introBox(
        data.step = 7,
        data.intro = "You can get the nomogram base on 3S score here.",
        h3("The patient-estimated nomogram based on the 3S score"),
        img(src = "nomo123.jpg",height=300)
        
      ),
      
      introBox(
        data.step = 8,
        data.intro = "Our contact is here.",
        hr(),
        h4("If you have any problems with using the web-app, or suggestions for improvement, please contact Weikaixin Kong at 1510307407@pku.edu.cn"),
        hr(),
        h4("Department of Molecular and Cellular Pharmacology"),
        h4("School of Pharmaceutical Sciences"),
        h4("Peking University Health Science Center"),
        h4("100191 Beijing"),
        h4("China")
      )
      
    )
    #width=500,
  )
)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output,session) {
  
  
  observe_helpers(help_dir = "helpfiles")
  
  
  datasetInput <- reactive({
    df <- read.csv("example.csv",
                   header = T,
                   sep = ",")
    df
  })
  genelist <- reactive({
    df2 <- read.csv("lassogenelist.csv",
                   header = T,
                   sep = ",")
    df2
  })
  
  
  
  pairlist <- reactive({
    df2 <- read.csv("lassopairlist.csv",
                    header = T,
                    sep = ",")
    df2
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("example data", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE,col.names = F)
    }
  )
  
  
  
  
  output$point_grade=renderText(
    {
      if(input$grade=="oneortwo")
      {
        "39"
      }
      else if(input$grade=="three")
      {
        "69"
      }
      else if(input$grade=="four")
      {
        "100"
      }
      
    }
  )
  
  output$point_primary=renderText(
    {
      
      if(input$primary=="yes")
      {
        "0"
      }
      else
      {
        "39"
      }
      
      
      
    }
  )
  
  
  
  datariskscore <- reactive({
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = T,
                   sep = ",",
                   row.names = 1)
    df$newcol=df[,1]
    
    rt=df[genelist()[,1],]
    tcgaPair=data.frame()
    sampleNum=ncol(rt)
    for(i in 1:(nrow(rt)-1)){
      for(j in (i+1):nrow(rt)){
        pair=ifelse(rt[i,]>rt[j,], 1, 0)
        pairRatio=sum(pair)/sampleNum
        if((pairRatio>-1) & (pairRatio<2)){
          rownames(pair)=paste0(rownames(rt)[i],"|",rownames(rt)[j])
          tcgaPair=rbind(tcgaPair, pair)
        }
      }
    }

    tcgaOut=rbind(ID=colnames(tcgaPair), tcgaPair)
    tcgaOut=tcgaOut[pairlist()[,1],]
    tcgaOut=t(tcgaOut)


    for(i in 1:ncol(tcgaOut))
    {
      tcgaOut[,i]=as.numeric(as.character(tcgaOut[,i]))
    }


    tcgaOut=as.data.frame(tcgaOut)
    for(i in 1:ncol(tcgaOut))
    {
      tcgaOut[,i]=as.numeric(as.character(tcgaOut[,i]))
    }


    m5= predict(object = rsfmodel5,newdata=tcgaOut)$predicted[1]
    m13= predict(object = rsfmodel13,newdata=tcgaOut)$predicted[1]
    m29= predict(object = rsfmodel29,newdata=tcgaOut)$predicted[1]
    m39= predict(object = coxmodel39,newdata=tcgaOut)[1]

    ensembleinput=as.data.frame(cbind(m5,m13,m29,m39))
    datariskscore=predict(object = rsfmodel_ensemble,newdata=ensembleinput)$predicted[1]-27
    as.numeric(datariskscore)
  })
  

  output$riskscore=renderText(
    {

      as.character(datariskscore())
      
      
    }
    
    
    
  )
  
  
  
  
  
  output$pointsss_risk=renderText(
    {
        scoreofrisk=((92-12)/3.5)*(as.numeric(datariskscore())/100-0)+12
        scoreofrisk=round(scoreofrisk,digits = 2)
        as.character(scoreofrisk)
      
      
    }
    
    
    
  )
  
  output$point_total=renderText(
    {
     
        sumvalue=0
        scoreofrisk=((92-12)/3.5)*(as.numeric(datariskscore())/100-0)+12

        
        sumvalue=sumvalue+scoreofrisk
        
        ######################################
        if(input$grade=="oneortwo")
        {
          sumvalue=sumvalue+39
        }
        else if(input$grade=="three")
        {
          sumvalue=sumvalue+69
        }
        else if(input$grade=="four")
        {
          sumvalue=sumvalue+100
        }

        #########################
        if(input$primary=="no")
        {
          sumvalue=sumvalue+39
        }

        ##############################
        
        sumvalue=round(sumvalue,digits = 2)
        as.character(sumvalue)
      }
      
      
    
  )
  
  #############we will calculate the final prob
  output$rate_year=renderText(
    {
        
      sumvalue=0
      scoreofrisk=((92-12)/3.5)*(as.numeric(datariskscore())/100-0)+12
      
      
      sumvalue=sumvalue+scoreofrisk
      
      ######################################
      if(input$grade=="oneortwo")
      {
        sumvalue=sumvalue+39
      }
      else if(input$grade=="three")
      {
        sumvalue=sumvalue+69
      }
      else if(input$grade=="four")
      {
        sumvalue=sumvalue+100
      }
      
      #########################
      if(input$primary=="no")
      {
        sumvalue=sumvalue+39
      }
      
      
      ##############################
      

          
          if(sumvalue>40&sumvalue<=60)
          {
            onevalue=((0.0626-0.0409)/20)*(sumvalue-40)+0.0409
            
            twovalue=((0.1816-0.1213)/20)*(sumvalue-40)+0.1213
            
            threevalue=((0.2642-0.1796)/20)*(sumvalue-40)+0.1796
          }
          else if(sumvalue>60&sumvalue<=80)
          {
            onevalue=((0.0954-0.0626)/20)*(sumvalue-60)+0.0626
            
            twovalue=((0.2671-0.1816)/20)*(sumvalue-60)+0.1816
            
            threevalue=((0.3785-0.2642)/20)*(sumvalue-60)+0.2642
          }
          else if(sumvalue>80&sumvalue<=100)
          {
            onevalue=((0.1440-0.0954)/20)*(sumvalue-80)+0.0954
            
            twovalue=((0.3823-0.2671)/20)*(sumvalue-80)+0.2671
            
            threevalue=((0.5216-0.3785)/20)*(sumvalue-80)+0.3785
          }
          else if(sumvalue>100&sumvalue<=120)
          {
            onevalue=((0.2142-0.1440)/20)*(sumvalue-100)+0.1440
            
            twovalue=((0.5261-0.3823)/20)*(sumvalue-100)+0.3823
            
            threevalue=((0.6811-0.5216)/20)*(sumvalue-100)+0.5216
          }
          else if(sumvalue>120&sumvalue<=140)
          {
            onevalue=((0.3117-0.2142)/20)*(sumvalue-120)+0.2142
            
            twovalue=((0.6858-0.5261)/20)*(sumvalue-120)+0.5261
            
            threevalue=((0.8300-0.6811)/20)*(sumvalue-120)+0.6811
          }
          else if(sumvalue>140&sumvalue<=160)
          {
            onevalue=((0.4396-0.3117)/20)*(sumvalue-140)+0.3117
            
            twovalue=((0.8338-0.6858)/20)*(sumvalue-140)+0.6858
            
            threevalue=((0.9359-0.8300)/20)*(sumvalue-140)+0.8300
          }
          else if(sumvalue>160&sumvalue<=180)
          {
            onevalue=((0.5925-0.4396)/20)*(sumvalue-160)+0.4396
            
            twovalue=((0.9381-0.8338)/20)*(sumvalue-160)+0.8338
            
            threevalue=((0.9858-0.9359)/20)*(sumvalue-160)+0.9359
          }
          else if(sumvalue>180&sumvalue<=200)
          {
            onevalue=((0.7514-0.5925)/20)*(sumvalue-180)+0.5925
            
            twovalue=((0.9866-0.9381)/20)*(sumvalue-180)+0.9381
            
            threevalue=((0.9986-0.9858)/20)*(sumvalue-180)+0.9858
          }
          else if(sumvalue>200&sumvalue<=220)
          {
            onevalue=((0.8844-0.7514)/20)*(sumvalue-200)+0.7514
            
            twovalue=((0.9988-0.9866)/20)*(sumvalue-200)+0.9866
            
            threevalue=((1.0000-0.9986)/20)*(sumvalue-200)+0.9986
          }
          else if(sumvalue>220&sumvalue<=240)
          {
            onevalue=((0.9647-0.8844)/20)*(sumvalue-220)+0.8844
            
            twovalue=((1.0000-0.9988)/20)*(sumvalue-220)+0.9988
            
            threevalue=1
          }
          else if(sumvalue>240)
          {
            onevalue=1
            twovalue=1
            threevalue=1
          }
          a=round(onevalue*100,2)
          b=round(twovalue*100,2)
          c=round(threevalue*100,2)
          paste0("1-year:",as.character(a),"%;  ","3-year:",as.character(b),"%;  ",
                 "5-year:",as.character(c),"%;")

        
        
        
      }
      
      
      
    
  )
  
  
  
  
  observeEvent(input$button1,
               introjs(session, options = list("nextLabel"="Next",
                                               "prevLabel"="Prev",
                                               "skipLabel"="End tour")
               )
  )
  
  
  
}

shinyApp(ui, server)
