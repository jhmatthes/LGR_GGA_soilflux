# Library of functions to streamline processing LGR GGA output into fluxes.

# format.LGR.output requires two input variables:
#   1. LGR.files = vector of LGR file names to be processed
#   2. LGR.match.files = index of LGR file names in LGR.files that match single day
# Output of this function is a list (LGR.dat) with LGR times, LGR dates, and LGR data (cols 1:14)

format.LGR.output <- function(LGR.files, LGR.match.files){
  # Aggregate LGR files if there is more than one file per day.
  if(length(match.LGR.files)>1){
    for(day.file in 1:length(match.LGR.files)){ 
      # Read LGR data and clean out garbage columns and rows.
      LGR.tmp  <- read.csv(LGR.files[match.LGR.files[day.file]],header=TRUE,skip=1,stringsAsFactors = FALSE)
      LGR.tmp  <- LGR.tmp[,1:14] #only use first 14 columns
      if(sum(LGR.tmp[,1]=="-----BEGIN PGP MESSAGE-----")>0){
        LGR.tmp  <- LGR.tmp[1:(which(LGR.tmp[,1]=="-----BEGIN PGP MESSAGE-----")-1),]
      }
      
      # Format messy column 1 into LGR date & time.
      LGR.date.tmp   <- strsplit(sapply(LGR.tmp$Time, as.character)," ")
      LGR.date <- LGR.time <- vector()
      for(i in 1:length(LGR.date.tmp)){
        LGR.date[i] <- LGR.date.tmp[[i]][3]
        LGR.time[i] <- LGR.date.tmp[[i]][4]
      }
      LGR.tmp.times <- times(LGR.time) 
      
      # Aggregate each file into a big daily file.
      if(day.file == 1){
        LGR.data  = LGR.tmp
        LGR.times = LGR.tmp.times
      } else {
        LGR.data  = rbind(LGR.data,LGR.tmp)
        LGR.times = c(LGR.times,LGR.tmp.times)
      }
    }
  } else { # If there is only one file for this day:
    # Read LGR data and clean out garbage columns and rows.
    LGR.data  <- read.csv(LGR.files[match.LGR.files],
                          header=TRUE,skip=1,stringsAsFactors = FALSE)
    LGR.data  <- LGR.data[,1:14] #only use first 14 columns
    if(sum(LGR.data[,1]=="-----BEGIN PGP MESSAGE-----")>0){
      LGR.data  <- LGR.data[1:(which(LGR.data[,1]=="-----BEGIN PGP MESSAGE-----")-1),]
    }
    
    # Format messy column 1 into LGR date & time.
    LGR.date.tmp   <- strsplit(sapply(LGR.data$Time, as.character)," ")
    LGR.date <- LGR.time <- vector()
    for(i in 1:length(LGR.date.tmp)){
      LGR.date[i] <- LGR.date.tmp[[i]][3]
      LGR.time[i] <- LGR.date.tmp[[i]][4]
    }
    LGR.times <- times(LGR.time) 
  }
  
  # Format output
  LGR.dat <- list()
  LGR.dat$time <- LGR.times
  LGR.dat$date <- LGR.date
  LGR.dat$data <- LGR.data
  
  return(LGR.dat)
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
