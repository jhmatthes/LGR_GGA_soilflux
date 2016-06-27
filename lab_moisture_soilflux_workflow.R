# This file outlines the workflow for the Lab Soil Flux Moisture Manipulation Experiment
# Summer 2016, Jaclyn Hatala Matthes

# Load libraries
library(chron)
library(ggplot2)
library(plyr)

# Load source files
lapply(list.files(path="R/",pattern = "[.]R$", recursive = TRUE, full.names=TRUE), source)

# 

# Set path to data 
data.path <- "data/"

# Set experiment constants: init = list of intialization parameters
init <- list ()
init$ts = 5 #timestep for LGR analyzer = 5 seconds in this case
init$jar.area  = 0.004118 # Cross section area of Mason jar in m2 - DOUBLE CHECK THIS
init$buffer.start   = 12 # LGR timepoint to start on (after 12*ts = 60 seconds)
init$buffer.end     = 24 # LGR timepoint to end on (was 34 --> changed to 2 minute window (24))
init$fit.timepoints = 12 # Number of LGR timepoints to include for flux slope calculation
init$plot.slope     = 1  # Make plots for CO2 concentration vs time for each jar, 0/1 = no/yes

# Calculate total mol in system: n = (PV)/(RT) = (101325*Vol_system)/(8.31441*298.15)
# Vol_system = LGR chamber + tubing + jar = 0.315 + 0.001 + 8 or 16 oz*(0.0295735 liters/oz)
oz.to.liters    = 0.0295735
init$vol.system = 0.315 + 0.001 + jar.dat$jar_size*oz.to.liters
init$n_total    = (101325*vol.system)/(8.31441*298.15) # n = (P * vol_system) / (R * T)

flux.dat <- read.csv("data/output/LGR_flux_output.csv",header=TRUE)

# Integrate total CO2 respired across the experiment by replicate
flux.dat.wide <- reshape(flux.dat, idvar = "id", timevar = "timepoint", direction = "wide")
timepoint.cols <- grep("CO2.flux.mass",colnames(flux.dat.wide))
flux.dat.wide <- flux.dat.wide[,c(1,timepoint.cols)]
meas.timestep <- 1
mass.CO2.total <- vector()
for(n in 1:total.reps){
  mass.CO2.total[n] <- (meas.timestep/2)*
    (sum(flux.dat.wide[n,grep("time_1",colnames(flux.dat.wide))]:timepoint.cols[length(timepoint.cols)]) + 
       sum(flux.dat.wide[n,grep("time_1",colnames(flux.dat.wide))]:timepoint.cols[length(timepoint.cols)-1]))
}

### Plotting stuff, in progress
# Plot fluxes by group

pdf(file="figures/sitegroups_byday.pdf")
# Calculate Day 0, group summaries
CO2.day0 <- summarySE(flux.dat[flux.dat$timepoint=="time_0",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 0, standard error of the mean
ggplot(CO2.day0, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() + geom_point() + labs(title = "Day 0, 6/20/16") + 
  ylab("CO2 Flux [mg-CO2 / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16)) 

# Calculate Day 1, group summaries
CO2.day1 <- summarySE(flux.dat[flux.dat$timepoint=="time_1",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 1, standard error of the mean
ggplot(CO2.day1, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 1, 6/22/16") + 
  ylab("CO2 Flux [mg-CO2 / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))

# Calculate Day 2, group summaries
CO2.day2 <- summarySE(flux.dat[flux.dat$timepoint=="time_2",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 2, standard error of the mean
ggplot(CO2.day2, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 2, 6/23/16") + 
  ylab("CO2 Flux [mg-CO2 / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))

# Calculate Day 3, group summaries
CO2.day3 <- summarySE(flux.dat[flux.dat$timepoint=="time_3",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 3, standard error of the mean
ggplot(CO2.day3, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 3, 6/24/16") + 
  ylab("CO2 Flux [mg-CO2 / (g-soil hr)]") + 
  scale_y_continuous(limits = c(-1, 16))

# Calculate Day 4, group summaries
CO2.day4 <- summarySE(flux.dat[flux.dat$timepoint=="time_4",], 
                      measurevar="CO2.flux.mass", groupvars=c("site","treat"))

# Plot Day 4, standard error of the mean
ggplot(CO2.day4, aes(x=site, y=CO2.flux.mass, colour=treat)) + 
  geom_errorbar(aes(ymin=CO2.flux.mass-se, ymax=CO2.flux.mass+se), width=.1) +
  geom_line() +
  geom_point() + labs(title = "Day 4, 6/25/16") + 
  ylab("CO2 Flux [mg-CO2 / (g-soil hr)]") + 
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

