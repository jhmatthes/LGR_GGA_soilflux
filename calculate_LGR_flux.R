# Process LGR output files into fluxes. 
# JHM (jaclyn.hatala.matthes@gmail.com), 06/22/16
# Three datasets are required to run this processing script:
#     1. Seperate directory (folder) with raw LGR files (.txt files) with 
#         dates in titles like gga_2016-06-22_f0000.txt
#     2. Text (.csv) file with all jar replicates listed, soil weight, jar size, 
#         and flux measurement start times for each timepoint (time_0, time_1, time_2 etc.)
#     3. Text (.csv) file that assigns dates for each timepoint (time_0 is 2016-06-20, etc.)

#remember to do this:
#install.packages("chron")

#load libraries
library(chron)
source("LGR_GGA_functionlib.R")

#set path to data 
data.path <- "data/"

# List files with the output from the LGR analyzer
LGR.files    <- list.files(paste(data.path,"LGR_data/",sep=""),full.names=TRUE)

# Load file with dates that correspond to timepoints
date.time <- read.csv(paste(data.path,"time_dates.csv",sep=""),header=TRUE, stringsAsFactors = FALSE)

# Load file with the experimental replicate start times, jar size, soil mass
all.data = readLines(paste(data.path,"jar_data.csv",sep=""))
skip.second  = all.data[-2]
jar.dat  = data.frame(read.csv(textConnection(skip.second), stringsAsFactors = FALSE, header = TRUE))
jar.dat$id   = paste(jar.dat$site,jar.dat$rep,sep="")
jar.time.col = grep("time",colnames(jar.dat))

# Experiment constants
ts = 5 #timestep for LGR analyzer = 5 seconds in this case
jar.area       = 0.004118 # Cross section area of Mason jar in m2 - DOUBLE CHECK THIS
total.reps     = nrow(jar.dat) #number of experimental replicates (total # jars)
buffer.start   = 12 # LGR timepoint to start on (after 12*ts = 60 seconds)
buffer.end     = 24 # LGR timepoint to end on (was 34 --> changed to 2 minute window (24))
fit.timepoints = 12 # Number of LGR timepoints to include for flux slope calculation
plot.slope     = 1  # Make plots for CO2 concentration vs time for each jar, 0/1 = no/yes

# Calculate total mol in system: n = (PV)/(RT) = (101325*Vol_system)/(8.31441*298.15)
# Vol_system = LGR chamber + tubing + jar = 0.315 + 0.001 + 8 or 16 oz*(0.0295735 liters/oz)
vol.system = 0.315 + 0.001 + jar.dat$jar_size*0.0295735
n_total    = (101325*vol.system)/(8.31441*298.15) # n = (P * vol_system) / (R * T)

# Always constants
ppm.to.mol = 10^-6 
C.mol.to.g = 12 #g/mol
sec.to.hour = 360 #s/hr

# Set up storage for things we're going to calculate
dates <- times <- reps <- CH4.sl <- CH4.r2 <- CO2.r2 <- CO2.sl <- flux.timepoint <- flux.dates <- flux.id <- flux.site <- flux.treat <- vector()
CH4.mol.rate <- CH4.flux.mass <- CH4.flux.area <- CO2.mol.rate <- CO2.flux.mass <- CO2.flux.area <- vector()

# Loop over each date and calculate the fluxes by matching the LGR files with the time files.
for(day in 1:(nrow(date.time)-sum(date.time$date==""))){
  
  # Find LGR data files that match each date (i.e., for time_0, time_1, etc.)
  match.LGR.files <- grep(date.time$date[day],LGR.files)
  
  # Aggregate and format LGR data
  LGR.data <- format.LGR.output(LGR.files, match.LGR.files) 
  
  # Check if there are replicates missing on this day. If so, only use reps with times.
  if(sum(jar.dat[jar.time.col[day]]=="")>0){
    jar.dat = data.frame(jar.dat[which(jar.dat[jar.time.col[day]]!=""),],
                         row.names = 1:length(which(jar.dat[jar.time.col[day]]!="")))
  } else { 
    # Re-load file, in case previous clause clipped jar.dat for a different day
    all.data = readLines(paste(data.path,"jar_data.csv",sep=""))
    skip.second  = all.data[-2]
    jar.dat  = data.frame(read.csv(textConnection(skip.second), stringsAsFactors = FALSE, header = TRUE))
    jar.dat$id   = paste(jar.dat$site,jar.dat$rep,sep="")
  }
  
  # Loop over jar replicates and calculate fluxes.
  for(jar in 1:nrow(jar.dat)){
    
    # Find replicate start time and match to LGR time file.
    rep.time <- times(jar.dat[jar,jar.time.col[day]])
    match.times <- which(abs(LGR.data$time-rep.time)==min(abs(LGR.data$time-rep.time)))
    
    # Define the LGR flux time period, get LGR data for that period.
    flux.period <- (match.times+buffer.start):(match.times+buffer.end) #find location of flux values 
    CO2.conc <- LGR.data$data$X.CO2.d_ppm[flux.period]   #grab the matching CO2 concentrations
    CH4.conc <- LGR.data$data$X.CH4.d_ppm[flux.period]   #grab the matching CH4 concentrations
    dat.time <- LGR.data$time[flux.period]               #grab the matching times
    
    rep.index <- (day-1)*total.reps+jar #index for storage: (day number-1) * jars each day + nth jar
    seconds   <- seq(1,(fit.timepoints+1)*ts,by=ts) #seconds for lm fit
    
    # Calculate flux by fitting a line
    lm.CH4    <- lm(CH4.conc ~ seconds)
    CH4.sl[rep.index]   <- summary(lm.CH4)$coefficients[2]
    CH4.r2[rep.index]   <- summary(lm.CH4)$r.squared
    
    # At same temp and pressure: V_CO2/V_T = n_CO2/n_total, n_CO2 = (ppm*10^-6) * n_total
    # The calculation is the same for the slope (ppm/s): (slope*10^-6) * n_total = mol/s 
    CH4.mol.rate[rep.index]  = CH4.sl[rep.index]*ppm.to.mol*n_total[jar] #mol/s
    CH4.flux.mass[rep.index] = (CH4.sl[rep.index]*ppm.to.mol*n_total[jar]*C.mol.to.g*sec.to.hour)/
      jar.dat$soil_weight[jar] #g CO2 / (g soil * hour)
    CH4.flux.area[rep.index] = (CH4.sl[rep.index]*ppm.to.mol*n_total[jar])/jar.area #umol/(m^2 * s)
    
    # Calculate CO2 flux by fitting a line
    lm.CO2    <- lm(CO2.conc ~ seconds)
    CO2.sl[rep.index]   <- summary(lm.CO2)$coefficients[2]
    CO2.r2[rep.index]   <- summary(lm.CO2)$r.squared
    
    # At same temp and pressure: V_CO2/V_T = n_CO2/n_total, n_CO2 = (ppm*10^-6) * n_total
    # The calculation is the same for the slope (ppm/s): (slope*10^-6) * n_total = mol/s 
    CO2.mol.rate[rep.index]  = CO2.sl[rep.index]*ppm.to.mol*n_total[jar] #mol/s
    CO2.flux.mass[rep.index] = (CO2.sl[rep.index]*ppm.to.mol*n_total[jar]*C.mol.to.g*sec.to.hour)/
      jar.dat$soil_weight[jar] #g CO2 / (g soil * hour)
    CO2.flux.area[rep.index] = (CO2.sl[rep.index]*ppm.to.mol*n_total[jar])/jar.area #umol/(m^2 * s)
    
    # Make plots of CO2.conc vs time to visually inspect
    if(plot.slope == 1){
      plot(seconds,CO2.conc,main=paste("Rep num: ",jar.dat$id[jar],sep=""))
    }
    
    # Replicate/date bookkeeping
    flux.dates[rep.index] = date.time$date[day]
    flux.timepoint[rep.index] = date.time$time[day]
    flux.id[rep.index]    = jar.dat$id[jar]
    flux.site[rep.index]  = jar.dat$site[jar]
    flux.treat[rep.index] = jar.dat$exp.group[jar]
    
  } #end jar loop
} #end measuring day loop

flux.dat <- data.frame(id=flux.id,date=flux.dates,timepoint=flux.timepoint,
                       site=flux.site,treat=flux.treat,group=paste(flux.site,flux.treat,sep=""),
                       CO2.flux.mass, CO2.flux.area, CO2.mol.rate, CO2.r2,
                       CH4.flux.mass, CH4.flux.area, CH4.mol.rate, CH4.r2)

#flux.dat[flux.dat$timepoint=="time_1",]

write.csv(flux.dat,file="LGR_flux_output.csv",row.names =FALSE)
