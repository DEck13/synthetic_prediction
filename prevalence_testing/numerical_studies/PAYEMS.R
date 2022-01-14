## 'PAYEMS' is the Federal Reserve of St. Louis' name for the monthly employment figure,
## (i.e. the number of persons employed in the US, in thousands)
## a time series that goes back to the 1950s.

### In this file, we produce two simple kinds output:

# (1) a set of shock times
# (2) a large matrix of data indexed by month, where the index
# includes the shock times as a proper subset

#significant digits
options(scipen = 7)

# load packages (not all of these may be necessary)
library("data.table")
library("dplyr")
library("tseries")
library("quantmod")
library('Rsolnp')
library('msos')
library('tikzDevice')
library('xtable')
# load packages
require('forecast')

## Time series under study: # persons on nonfarm payrolls
getSymbols("PAYEMS", src = 'FRED')

###          COVARIATES            ###
###           START                ###

## https://fred.stlouisfed.org/series/W825RC1: transfer receipts
getSymbols("W825RC1", src = 'FRED')

## https://fred.stlouisfed.org/series/RPI: real personal income
getSymbols("RPI", src = 'FRED')

## https://fred.stlouisfed.org/series/PCE: Personal Consumption Expenditures
getSymbols("PCE", src = 'FRED')

## https://fred.stlouisfed.org/series/INDPRO: Industrial Production
getSymbols("INDPRO", src = 'FRED') 

## https://fred.stlouisfed.org/series/CPIAUCSL: Consumer Price Index
getSymbols("CPIAUCSL", src = 'FRED')

## https://fred.stlouisfed.org/series/FEDFUNDS: Federal Funds Index
getSymbols("FEDFUNDS", src = 'FRED')

## https://fred.stlouisfed.org/series/LNS12000031: Black Employment Count
getSymbols("LNS12000031", src = 'FRED') ######## NOTE: IF THIS is not super important as variable, let's drop it

###          COVARIATES          ###
###             END              ###



### CONSTRUCTION OF THE DONOR POOL ###
###           START                ###

#https://www.federalreservehistory.org/essays/oil-shock-of-1973-74
#https://www.federalreservehistory.org/essays/recession-of-1981-82
#credit-control program initiated in March 1980 by the Carter administration
#https://www.history.com/news/us-economic-recessions-timeline
#https://fred.stlouisfed.org/series/FEDFUNDS#0

# As a rule, when the shock occurs mid-month, we take that to be the shock-time, even though
# the shock effect is distributed across that month as well as the following month(s)
shock_time_vec <- c('1957-04-01', ## Flu hits US
                    '1958-08-01', ## Fed Funds rate increases from .68 to 1.53 in one month
                    '1973-10-01', ## OAPEC oil embargo begins
                    '1980-03-01', ## program announced on March 14th by Carter
                    '2001-09-01', ## 9/11 attacks occur 1/3 of way into month
                    '2020-03-01') ## COVID shutdown - included in the time series of interest

                    # These are all T* points not T*+1.

### CONSTRUCTION OF THE DONOR POOL ###
###           END                  ###


## Now, we take the outcome variable and covariates and smash them into a 
## large matrix, including lags of the covariates.
datasets <- list(W825RC1, RPI, PCE, INDPRO, CPIAUCSL, FEDFUNDS, LNS12000031)
df <- PAYEMS
for (i in 1:length(datasets)) {df <- merge(df, datasets[[i]])}

# Hit every column with the differenced-log transformation
difflog_df <- data.frame(diff(as.matrix(log(df))))

# We create lags of the covariates
# https://stackoverflow.com/questions/38119225/debugging-function-to-create-multiple-lags-for-multiple-columns-dplyr

difflog_df.lag <- shift(difflog_df, n=1:2, give.names = T)  ##column indexes of columns to be lagged as "[,startcol:endcol]", "n=1:3" specifies the number of lags (lag1, lag2 and lag3 in this case)

# We we combine and original series and the lags
merged <- bind_cols(difflog_df, difflog_df.lag)

#Now add the row names
row.names(merged) <- row.names(difflog_df)

# We do not need any rows prior to 1954-07-01
merged <- merged[row.names(merged) >= '1954-07-01', ]

#Now, due to the release date for each of these monthly series,
#we unfortunately cannot use the same-month data points for some of 
#these columns.  We now drop them...

merged <- subset(merged, select = -c(W825RC1, 
                                    RPI, 
                                    PCE, 
                                    INDPRO, 
                                    CPIAUCSL,
                                    LNS12000031))


#Finally, we have missing data in the late 1950s, so we are faced with a choice:

#(1) We can drop all rows with NA entries, which will remove the 1950s donors.
complete_cases_merged <- merged[complete.cases(merged),]
paste('This is a dataset of dimension', dim(complete_cases_merged)[1], 'by', dim(complete_cases_merged)[2])

#(2) We can drop all columns with NA values, which will drop some covariates, but keep all 5 donors.
no_NA_cols_merged <- merged[ , colSums(is.na(merged)) == 0]
paste('This is a dataset of dimension', dim(no_NA_cols_merged)[1], 'by', dim(no_NA_cols_merged)[2])


# And we're done!  Either of the two datasets above will work.  It's a question of whether we 
# want 5 donors with fewer covariates, or fewer donors with more covariates.
