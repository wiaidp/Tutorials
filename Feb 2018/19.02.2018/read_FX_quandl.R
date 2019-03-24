#install.packages("shiny")
#install.packages("Quandl")
#install.packages("sendmailR")

library(shiny)
library(Quandl)
library(sendmailR)


sender <- "martingjesus@gmail.com"
r <- c("martingjesus@gmail.com", "a.sidler@edrofine.ch", "MBenshemesh@gmail.com")
sender <- "wlmr@zhaw.ch"
r <- c("wlmr@zhaw.ch")
subject <- "EMA Signal change"
pass <- ""

QuandlAPIKey <- "cxvP9tfN13a9tvMmsXag"

Quandl.api_key(QuandlAPIKey)

# This code will provide a backtest in Excel of a strategy based in a exponential moving average crossover. 
# The inputs we have to provide are:
# FromDate: The starting date of the backtest.
# ToDate: The finishing date of the backtest.
# Instrument: The instrument that will be used to backtest the strategy.
# EMALength: The length used for the calculation of the exponential moving average.
# TradeSize: The # of units/shares/contracts we are purchasing of the instrument each time we trade.

# The strategy logic is the well-known crossover. We buy when the price crosses over the exponential moving average and we
# sell when the price crosses under the exponential moving average. In this case, the strategy is of the type stop&reverse or
# always-in. The signal that triggers a sell short or buy also works for a buy to cover or sell.

# The output of this function is a data frame object containing the backtest information and a csv file that will be created
# in the current working directory.


# Reading data and changing the date format
#Data <- read.csv(file= paste(getwd(),"/data/EURUSD.txt", sep = ""))
Data <- Quandl("ECB/EURUSD")
Data$Date <- rev(Data$Date)
Data$Value <- rev(Data$Value)
Data$Date <- as.Date(Data$Date, format = "yyyy/mm/dd")

ui <- fluidPage(
      title = "EMA Backtest",
        
      fluidRow(column(4,
        
          sliderInput(inputId = "EMALength", label = "EMA Length", value = 14, 
                            min = 1, max = 500),
          numericInput(inputId = "Alpha", label = "Alpha", value = 0.2, min = NA, max = NA, step = NA,
                       width = NULL)
        ),
        
        column(4,
                dateInput(inputId = "FromDate", label = "From Date", value = "2008/01/01", min = NULL, max = NULL,
                          format = "yyyy/mm/dd", startview = "month", weekstart = 0,
                          language = "en", width = NULL),
                dateInput(inputId = "ToDate", label = "To Date", value = "2015/01/01", min = NULL, max = NULL,
                          format = "yyyy/mm/dd", startview = "month", weekstart = 0,
                          language = "en", width = NULL)
        ),
        
        column(4,
          
         numericInput(inputId = "TradeSize", label = "Trade size", value = 100000, min = NA, max = NA, step = NA,
                            width = NULL),
         checkboxInput(inputId = "TradeArrows", label = "Trade Arrows", value = FALSE, width = NULL),

        numericInput(inputId = "ArrowWidth", label = "Arrow width", value = 4, min = NA, max = NA, step = NA,
                     width = NULL)
        )
      ),

        plotOutput("plot")
      )

   


server <- function(input,output, session){
  
  if (F)
  {
    observe({
      
      # Invalidate this observer every hour (360000 milliseconds)
      invalidateLater(360000, session)
      
      Data <- Quandl("ECB/EURUSD")
      Data$Date <- rev(Data$Date)
      Data$Value <- rev(Data$Value)
      Data$Date <- as.Date(Data$Date, format = "yyyy/mm/dd")
      
      body <- paste("There has been a change in the signal of the EMA algorithm, now the position
    is", toString(Data$SignalChange[0]))
      
       sendmail(
         from = sender,
         to = r,
         subject= subject,
         body = body,
         smtp = list(host.name = "mail.mydomain.com", port = 465, user.name = sender, passwd = alas9781, ssl = TRUE),
         authenticate = TRUE,
         send = TRUE
       )
  
    })
  }
  output$plot <- renderPlot({
    
    # Subsetting rows within the dates inputted and removing unnecessary columns in the dataset
    Data <- subset(Data, Data$Date >= as.Date(input$FromDate)
                   & Data$Date <= as.Date(input$ToDate))
    
    
    # Weighting multiplier
    W = input$Alpha
    
    # Calculating a Simple moving average
    SMACalc = NA
    for (i in 1:length(Data$Value)){
      if (i>= input$EMALength)
        SMACalc[i] <- round(mean(c(Data$Value[i-input$EMALength+1:i])), 5) else
          SMACalc[i] <- NA
        
    }
    
    # Calculating a Exponential moving average
    EMACalc = NA
    for (i in 1:length(Data$Value)){
      if (i> input$EMALength)
        EMACalc[i] <- round((Data$Value[i] - EMACalc[i-1])*W + EMACalc[i-1],5) else
          if (i == input$EMALength)
            EMACalc[i] <- SMACalc[i] else
              EMACalc[i] <- NA
            
    }  
    
    # Inspecting Data
    head(Data)
    
    # Resetting row names
    rownames(Data) <- seq(length=nrow(Data))
    
    # Calculating the EMA Values and inserting the in the Data dataframe
    EMA <- EMACalc
    Data$EMA <- EMA
    
   
    # Checking closing prices vs. EMA values to see if we should be long (price above EMA) 
    # or short (price below EMA) 
    Data$PriceAboveEMA <- NA
    for (i in 1:nrow(Data)){
      if (is.na(Data$EMA[i]) == FALSE) {
        if (Data$Value[i] > Data$EMA[i])
          Data$PriceAboveEMA[i] <- 1
        else
          Data$PriceAboveEMA[i] <- -1
      }
    }
    
    # Calculating daily PL 
    Data$DailyPL <- NA
    for (i in 2:nrow(Data)){
      if (is.na(Data$PriceAboveEMA[i-1]) == FALSE) {
        Data$DailyPL[i] <- round(input$TradeSize*(Data$Value[i]-Data$Value[i-1])*Data$PriceAboveEMA[i-1],2)
      }
    }
    
    # Calculating the cumulative equity gained (or loss)
    Data$CumPL <- NA
    Data$CumPL <- cumsum(c(rep(0, input$EMALength+1), Data$DailyPL[(input$EMALength+2):nrow(Data)]))
    
    # Signal change. This column marks the actual trades with the number of lots included.
    Data$SignalChange <- 0
    Data$SignalChange[2:length(Data$PriceAboveEMA)] <- Data$PriceAboveEMA[2:length(Data$PriceAboveEMA)] - Data$PriceAboveEMA[1:(length(Data$PriceAboveEMA)-1)]
    
    # The function  plots the EMA alongside the price series with markers for 
    # trades
    par(mfrow = c(2,1))
    plot(y = Data$CumPL/1000, x = as.Date(Data$Date), type = "l", main = "Equity", xlab = "Daily time series",
         ylab = "$ (k)", col = "darkgreen", lwd = 2)
    plot(y = Data$Value, x = as.Date(Data$Date), type = "l", main = "Price series + EMA",
         xlab = "Daily time series", ylab = "Price", col = "black", lwd = 2)
    lines(y = Data$EMA[!is.na(Data$EMA)], x = as.Date(Data$Date[!is.na(Data$EMA)]), type = "l", col = "darkgray", lwd = 2)
    
    if (input$TradeArrows == TRUE){
    arrows(x0 = Data$Date[Data$SignalChange > 0], y0 = (Data$Value[Data$SignalChange > 0]-0.04), x1 = Data$Date[Data$SignalChange > 0], 
           y1 = (Data$Value[Data$SignalChange > 0]-0.01), lwd = input$ArrowWidth, col = "blue")
    arrows(x0 = Data$Date[Data$SignalChange < 0], y0 = (Data$Value[Data$SignalChange < 0]+0.04), x1 = Data$Date[Data$SignalChange < 0], 
           y1 = (Data$Value[Data$SignalChange < 0]+0.01), lwd = input$ArrowWidth, col = "red")
    }
    
})
}

shinyApp(ui = ui, server = server)
