# This file outlines the workflow for the Lab Soil Flux Moisture Manipulation Experiment
# Three datasets are required to run this processing script:
#     1. Seperate directory (folder) with raw LGR files (.txt files) with 
#         dates in titles like gga_2016-06-22_f0000.txt
#     2. Text (.csv) file with all replicates listed (rep.data), soil weight, jar size, 
#         and flux measurement start times for each timepoint (time_0, time_1, time_2 etc.)
#     3. Text (.csv) file that assigns dates for each timepoint (date.data: time_0 is 2016-06-20, etc.)
# Summer 2016, Jaclyn Hatala Matthes

# Load libraries
library(chron)
library(ggplot2)
library(plyr)

# Load R function source files (everything from project path R/)
lapply(list.files(path="R/",pattern = "[.]R$", recursive = TRUE, full.names=TRUE), source)

# Set input data 
data.path = "data/input/" 
rep.data  = "jar_data.csv"
date.data = "time_dates.csv"

# Set experiment constants: init = list of intialization parameters
init <- list ()
init$method = "jar" #valid values = "chamber" (field) or "jar (lab)
init$lgr.ts = 5     #timestep for LGR analyzer = 5 seconds in this case
init$surf.area  = 0.004118 # x-section area of measurement chambers (Mason jar in m2)
init$flux.start = 12 # LGR timepoint to start on (after 12*ts = 60 seconds)
init$flux.end   = 24 # LGR timepoint to end on (was 34 --> changed to 2 minute window (24))
init$plot.slope = 1  # Make plots for CO2 concentration vs time for each jar, 0/1 = no/yes

# Load file with the experimental replicate data.
# Here, we need this for the volume for each replicate, either chamber/jar (called rep_volume)
all.data = readLines(paste(data.path,"jar_data.csv",sep=""))
skip.second  = all.data[-2]
jar.dat  = data.frame(read.csv(textConnection(skip.second), stringsAsFactors = FALSE, header = TRUE))

# Calculate total mol in system: n = (PV)/(RT) = (101325*Vol_system)/(8.31441*298.15)
# Vol_system = LGR chamber + tubing + jar = 0.315 + 0.001 + 8 or 16 oz*(0.0295735 liters/oz)
oz.to.liters    = 0.0295735
chamber.vol     = jar.dat$rep_volume*oz.to.liters # Volume of chambers/jars
vol.system      = 0.315 + 0.001 + chamber.vol
init$n.total    = (101325*vol.system)/(8.31441*298.15) # n = (P * vol_system) / (R * T)

# Run LGR processing function to calculate fluxes
# Also requires format_LGR_output.R to be loaded
flux.dat <- calculate_LGR_flux(data.path,rep.data,date.data,init)


# Integrate total CO2 respired across the experiment by replicate
flux.dat.wide <- reshape(flux.dat, idvar = "id", timevar = "timepoint", direction = "wide")
timepoint.cols <- grep("CO2.flux.mass",colnames(flux.dat.wide))
flux.dat.new  <- flux.dat.wide[,c(1,15,16,17,timepoint.cols)]
colnames(flux.dat.new)[1:4] <- c("id", "site", "treat", "group")
timepoint.new <- grep("time_",colnames(flux.dat.new))
meas.timestep <- 24 #hours
mass.CO2.total <- vector()
for(n in 1:nrow(jar.dat)){
  mass.CO2.total[n] <- (meas.timestep/2)*
    (sum(flux.dat.new[n,grep("time_1",colnames(flux.dat.new)):timepoint.new[length(timepoint.new)]]) + 
       sum(flux.dat.new[n,(grep("time_1",colnames(flux.dat.new))+1):timepoint.new[length(timepoint.new)-1]]))
}

flux.dat.new$CO2.total <- mass.CO2.total 

### Plotting stuff, in progress - MOVE TO SEPARATE FILE
# Plot fluxes by group

pdf(file="figures/CO2flux_sitegroups_byday.pdf")
# Calculate Day 0, group summaries
CO2.day0 <- summarySE(flux.dat[flux.dat$timepoint=="time_0",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 0, standard error of the mean
ggplot(CO2.day0, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() + geom_point() + labs(title = "Day 0, 6/20/16") + 
  ylab("CO2 Flux [mg-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16)) 

# Calculate Day 1, group summaries
CO2.day1 <- summarySE(flux.dat[flux.dat$timepoint=="time_1",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 1, standard error of the mean
ggplot(CO2.day1, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 1, 6/22/16") + 
  ylab("CO2 Flux [mg-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))

# Calculate Day 2, group summaries
CO2.day2 <- summarySE(flux.dat[flux.dat$timepoint=="time_2",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 2, standard error of the mean
ggplot(CO2.day2, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 2, 6/23/16") + 
  ylab("CO2 Flux [mg-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))

# Calculate Day 3, group summaries
CO2.day3 <- summarySE(flux.dat[flux.dat$timepoint=="time_3",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 3, standard error of the mean
ggplot(CO2.day3, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 3, 6/24/16") + 
  ylab("CO2 Flux [mg-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))

# Calculate Day 4, group summaries
CO2.day4 <- summarySE(flux.dat[flux.dat$timepoint=="time_4",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 4, standard error of the mean
ggplot(CO2.day4, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 4, 6/25/16") + 
  ylab("CO2 Flux [mg-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))
dev.off()

# Time series plots by site
pdf(file="figures/site_timeseries.pdf")
sites <- unique(flux.dat$site)
for(s in 1:length(sites)){
  g = ggplot(subset(flux.dat,site %in% sites[s]), aes(x = date, y=CO2.flux.mass)) + geom_boxplot()
  print(g + aes(colour = factor(treat))  + scale_y_continuous(limits = c(-1, 16)) + labs(title = sites[s]))
}
dev.off()

# Calculate Total CO2 flux, group summaries
CO2.total <- summarySE(flux.dat.new, 
                      measurevar="CO2.total", groupvars=c("site","treat"))

# Plot Day 0, standard error of the mean
ggplot(CO2.total, aes(x=site, y=CO2.total, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.total-se, ymax=CO2.total+se), width=.1) +
  geom_line() + geom_point() + labs(title = expression('Total CO'[2]*' Respired') )+ 
  ylab(expression('CO'[2]*' Flux [mg-C / g-soil]')) + 
  scale_y_continuous(limits = c(0, 700)) 

# CH4 flux summary
pdf(file="figures/CH4flux_sitegroups_byday.pdf")
# Calculate Day 0, group summaries - NONE HAVE SIGNIFICANT FLUX (all R2 < 0.5)
# Calculate Day 1, group summaries
CH4.day1 <- summarySE(flux.dat[flux.dat$timepoint=="time_1"&flux.dat$CH4.r2>0.5,], 
                      measurevar="CH4.flux.mass", groupvars=c("site","treat"))

# Plot Day 1, standard error of the mean
ggplot(CH4.day1, aes(x=site, y=CH4.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CH4.flux.mass-se, ymax=CH4.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 1, 6/22/16") + 
  ylab("CH4 Flux [ug-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-0.3, 0.01))

# Calculate Day 2, group summaries
CH4.day2 <- summarySE(flux.dat[flux.dat$timepoint=="time_2"&flux.dat$CH4.r2>0.5,], 
                      measurevar="CH4.flux.mass", groupvars=c("site","treat"))

# Plot Day 1, standard error of the mean
ggplot(CH4.day2, aes(x=site, y=CH4.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CH4.flux.mass-se, ymax=CH4.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 2, 6/23/16") + 
  ylab("CH4 Flux [ug-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-0.3, 0.01))

# Calculate Day 3, group summaries
CH4.day3 <- summarySE(flux.dat[flux.dat$timepoint=="time_3"&flux.dat$CH4.r2>0.5,], 
                      measurevar="CH4.flux.mass", groupvars=c("site","treat"))

# Plot Day 1, standard error of the mean
ggplot(CH4.day3, aes(x=site, y=CH4.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CH4.flux.mass-se, ymax=CH4.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 3, 6/24/16") + 
  ylab("CH4 Flux [ug-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-0.3, 0.01))

# Calculate Day 1, group summaries
CH4.day4 <- summarySE(flux.dat[flux.dat$timepoint=="time_4"&flux.dat$CH4.r2>0.5,], 
                      measurevar="CH4.flux.mass", groupvars=c("site","treat"))

# Plot Day 1, standard error of the mean
ggplot(CH4.day4, aes(x=site, y=CH4.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CH4.flux.mass-se, ymax=CH4.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 4, 6/22/16") + 
  ylab("CH4 Flux [ug-C / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-0.3, 0.01))

dev.off()

