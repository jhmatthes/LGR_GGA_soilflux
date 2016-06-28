# Process LGR output files into fluxes. 
# JHM (jaclyn.hatala.matthes@gmail.com), 06/22/16
# Three datasets are required to run this processing script:
#     1. Seperate directory (folder) with raw LGR files (.txt files) with 
#         dates in titles like gga_2016-06-22_f0000.txt
#     2. Text (.csv) file with all jar replicates listed, soil weight, jar size, 
#         and flux measurement start times for each timepoint (time_0, time_1, time_2 etc.)
#     3. Text (.csv) file that assigns dates for each timepoint (time_0 is 2016-06-20, etc.)
# This function uses format_LGR_output.R, which cleans and aggregates LGR data for each day. 

calculate_LGR_flux <- function(data.path,init){
 
  # Always constants
  ppm.to.mol = 10^-6 
  C.mol.to.g = 12 #g/mol
  sec.to.hour = 360 #s/hr
  
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
  
  total.reps     = nrow(jar.dat) #number of experimental replicates (total # jars)
  
  # Set up storage for things to calculate
  dates <- times <- reps <- CH4.sl <- CH4.r2 <- CO2.r2 <- CO2.sl <- flux.timepoint <- flux.dates <- flux.id <- flux.site <- flux.treat <- vector()
  CH4.mol.rate <- CH4.flux.mass <- CH4.flux.area <- CO2.mol.rate <- CO2.flux.mass <- CO2.flux.area <- vector()
  
  # Loop over each date and calculate the fluxes by matching the LGR files with the time files.
  for(day in 1:(nrow(date.time)-sum(date.time$date==""))){
    
    # Print message to keep track of progree
    print(paste("Working on: day ",day," of ",(nrow(date.time)-sum(date.time$date=="")),sep=""))
    
    # Find LGR data files that match each date (i.e., for time_0, time_1, etc.)
    match.LGR.files <- grep(date.time$date[day],LGR.files)
    
    # Aggregate and format LGR data
    LGR.data <- format_LGR_output(LGR.files, match.LGR.files) 
    
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
      flux.period <- (match.times+init$flux.start):(match.times+init$flux.end) #find location of flux values 
      CO2.conc <- LGR.data$data$X.CO2.d_ppm[flux.period]   #grab the matching CO2 concentrations
      CH4.conc <- LGR.data$data$X.CH4.d_ppm[flux.period]   #grab the matching CH4 concentrations
      dat.time <- LGR.data$time[flux.period]               #grab the matching times
      
      rep.index <- (day-1)*total.reps+jar #index for storage: (day number-1) * jars each day + nth jar
      seconds   <- seq(1,((init$flux.end - init$flux.start)+1)*init$lgr.ts,by=init$lgr.ts) #seconds for lm fit
      
      # Calculate flux by fitting a line
      lm.CH4    <- lm(CH4.conc ~ seconds)
      CH4.sl[rep.index]   <- summary(lm.CH4)$coefficients[2]
      CH4.r2[rep.index]   <- summary(lm.CH4)$r.squared
      
      # At same temp and pressure: V_CO2/V_T = n_CO2/n_total, n_CO2 = (ppm*10^-6) * n_total
      # The calculation is the same for the slope (ppm/s): (slope*10^-6) * n_total = mol/s 
      CH4.mol.rate[rep.index]  = CH4.sl[rep.index]*ppm.to.mol*init$n.total[jar] #mol/s
      
      if(init$method == "jar"){ # Calculate flux per mass soil
        CH4.flux.mass[rep.index] = (CH4.sl[rep.index]*ppm.to.mol*init$n.total[jar]*C.mol.to.g*sec.to.hour)/
          jar.dat$soil_weight[jar] #g CO2 / (g soil * hour)
        CH4.flux.area[rep.index] = NA
      } else if (init$method == "chamber"){ # Calculate flux per surface area
        CH4.flux.mass[rep.index] = NA
        CH4.flux.area[rep.index] = (CH4.sl[rep.index]*ppm.to.mol*init$n.total[jar])/init$jar.area #umol/(m^2 * s)
      }
      
      # Calculate CO2 flux by fitting a line
      lm.CO2    <- lm(CO2.conc ~ seconds)
      CO2.sl[rep.index]   <- summary(lm.CO2)$coefficients[2]
      CO2.r2[rep.index]   <- summary(lm.CO2)$r.squared
      
      # At same temp and pressure: V_CO2/V_T = n_CO2/n_total, n_CO2 = (ppm*10^-6) * n_total
      # The calculation is the same for the slope (ppm/s): (slope*10^-6) * n_total = mol/s 
      CO2.mol.rate[rep.index]  = CO2.sl[rep.index]*ppm.to.mol*init$n.total[jar] #mol/s
      
      if(init$method == "jar"){ # Calculate flux per mass soil
      CO2.flux.mass[rep.index] = (CO2.sl[rep.index]*ppm.to.mol*init$n.total[jar]*C.mol.to.g*sec.to.hour)/
        jar.dat$soil_weight[jar] #g CO2 / (g soil * hour)
      CO2.flux.area[rep.index] = NA
      } else if(init$method == "chamber"){
        CO2.flux.mass[rep.index] = NA
        CO2.flux.area[rep.index] = (CO2.sl[rep.index]*ppm.to.mol*init$n.total[jar])/init$jar.area #umol/(m^2 * s)
      }

      # Make plots of CO2.conc vs time to visually inspect
      if(init$plot.slope == 1){
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
  
  flux.dat <- flux.dat[!is.na(flux.dat$id),] # Remove any NA rows if some days are missing reps
  flux.dat$CO2.flux.mass <- flux.dat$CO2.flux.mass*1000 #scale to mg-C / (g-soil * h)
  flux.dat$CH4.flux.mass <- flux.dat$CH4.flux.mass*10^6 #scale to ug-C / (g-soil * h)
  
  # Save calculated fluxes
  write.csv(flux.dat,file="data/output/LGR_flux_output.csv",row.names =FALSE)
  
  return(flux.dat)
}
# # Add second header for units info to flux.dat 
# comment(flux.dat$id)    = "site * location * forest * treatment * rep"
# comment(flux.dat$date)  = "YYYY-MM-DD"
# comment(flux.dat$timepoint) = "experiment measurement time"
# comment(flux.dat$site)  = "site * location * forest"
# comment(flux.dat$treat) = "experimental treatment: A = dry, B = med, C = wet"
# comment(flux.dat$group) = "site * location * forest * treatment"
# comment(flux.dat$CO2.flux.mass) = "mg-C / (g-soil * h)"
# comment(flux.dat$CO2.flux.area) = "umol / (m^2 * s)"
# comment(flux.dat$CO2.mol.rate)  = "mol / s"
# comment(flux.dat$CO2.r2)  = "fit of flux slope"
# comment(flux.dat$CH4.flux.mass) = "ug-C / (g-soil * h)"
# comment(flux.dat$CH4.flux.area) = "umol / (m^2 * s)"
# comment(flux.dat$CH4.mol.rate)  = "mol / s"
# comment(flux.dat$CH4.r2)  = "fit of flux slope"
# 
# write.csv3(flux.dat,file="data/output/LGR_flux_output_test.csv")
