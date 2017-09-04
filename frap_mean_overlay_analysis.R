rm(list = ls(all = TRUE))

# Functions for Overlay Analysis
## The following script is designed to analyse Fluorescence Recovery after Photobleaching experiments in "batch" mode. This means, that you as the user can define the path on the file system of your computer to the measurement data in .csv-format. Preferably, you stored your data in seperate folders for each individual experiment. The script initially loops through the folders of the file path you give to it and loads all of your data in a data frame called "frap". Then you need to define the rows in your data you want to evaluate, which row represents the background and which the correction ROI of observational photobleaching.
## In the variable "postIntervall" you need to tell the script the index of your measurement data, where your post bleaching intervall startet and ended. 
## In the end all of the evaluated data are stored in .csv-file (time, mean, standard deviation, standard error)

#### USER INPUT 

# Please add the file path to the FOLDER, that contains your measurements.
path = "/Users/tomkache/Documents/Studium/Biologie/SS 2017 B. Sc. Biology FSU Jena/Bachelorarbeit/FRAP-Experimente/S100A11-EGFP/S100A11-EGFP_woBLM_29_07_17/"

# What are your csv-files called? 
file_name = "/test.csv"

# In which columns of your csv-data are the data of the measurement (bleached ROI), the background and the fading (ROI, that contains data about observational photobleaching)?
colMeas = 2
colBack = 4
colFade = 5

# What interval of the prebleaching phase would you like to use for correction? (Please type in the row numbers from:to)
preInterVal = 20:100

# At which rownumber does the postbleaching phase start? 
postInterval = 103

# Add the name of your measurement. The averaged FRAP data and the results of the curve fitting will be stored at the location on your coputer, that you provided in the "path" variable.
measurement_name = "S100A11_woBLM_Ncl"


#### END OF USER INPUT


# Function definitions


frap.normalize = function(intensity, pre_phase = 20:100) {
  
  min <- intensity[which.min(intensity)]
  max <- mean(intensity[pre_phase])
  
  scaled <- for(i in 1:length(intensity)){
    intensity[i] <- (intensity[i] - min)/(max - min)
  }
  scaled = intensity
  return(scaled)
}

frap.se = function(sd, sampleSize){
  
  return((sd/sqrt(sampleSize)))
  
}

frap.kinetics = function(t, tau, mf, corr, tau2, B, y0){
  
  y = corr + (mf - corr) * (1 - exp(-t / tau)) * (y0 + B * exp(-t / tau2))
  return(y)

}

frap.BBCorr = function(measurement, background, fading){
  
  measurement = sweep(measurement, 1, background)
  fading = sweep(fading, 1, background)
  measurent = sweep(measurement, 1, fading)
  
  return(measurement)
}

# Loading the required packages

require(readr)
require(ggplot2)
require(drc)

####### START


# Batch Import of the Measurement Data


# The next step is to bring together all the measurement rows. 
## The time axis is skipped and only imported initially (see above) as the time intervals are all the same.

folder = 1

frap = read_csv(paste0(path,
                       folder,
                       file_name),
                skip = 1)


folder = folder + 1

while(folder < 30){
  path_loop = paste0(path, 
                folder, 
                file_name)
  
  dat <- read_csv(path_loop,
                 skip = 1,
                 col_types = cols(`Axis [s]` = col_skip()))
  frap = cbind(frap, dat)
  
  folder = folder + 1
  
}

# Seperate the ROIs of the measurement, the background and the observational photobleaching

frap_time = frap[[1]]

i = colMeas
frap_meas = data.frame(frap[[i]])

while(i <= length(frap)){
  
  dat = frap[[i]]
  frap_meas = cbind(frap_meas, dat)
  
  i = i+4
  
}

i = colBack
frap_back = data.frame(frap[[i]])

while(i <= length(frap)){
  
  dat = frap[[i]]
  frap_back = cbind(frap_back, dat)
  
  i = i+4
  
}



i = colFade
frap_fade = data.frame(frap[[i]])
while(i <= length(frap)){
  
  dat = frap[[i]]
  frap_fade = cbind(frap_fade, dat)
  
  i = i+4
  
}




# Correcting for Background & Observational Photobleaching
## Subtracting Background from Measurement and Fading ROIs (cb = corrected for background)

frap_meas_cb = frap_meas - frap_back
frap_fade_cb = frap_fade - frap_back

## Removing Observational Photobleaching by division (measurement / fading)

frap_meas_cor = frap_meas_cb / frap_fade_cb


## Normalising
frap_meas_norm = as.data.frame(frap.normalize(frap_meas_cor[[1]], pre_phase = preInterVal))
i = 2
while(i <= length(frap_meas_cor)){
  
  dat = frap.normalize(frap_meas_cor[[i]], pre_phase = preInterVal)
  frap_meas_norm = cbind(frap_meas_norm, dat)
  i = i + 1
  
}


# Calculating the mean, standard deviaton and standard error to obtain the "average" curve
## We will loop through each row of the normalized data frame and calculate these values

frap_meas_mean = rowMeans(frap_meas_norm)
frap_meas_sd = apply(frap_meas_norm,
                    1,
                    sd)
sample_size = length(frap_meas_norm)

frap_meas_se = apply(as.matrix(frap_meas_sd),
                    1:2,
                    function(x) x/sqrt(sample_size))


# Shifting the Curve: Time begins at Bleaching Event

#frap_time = frap_time - frap_time[which.min(frap_meas_mean)]

# Getting all the data in one data frame (time, mean of measurement, sd, se)

frap_mean = data.frame(frap_time,
                       frap_meas_mean,
                       frap_meas_sd,
                       frap_meas_se)


# Overlayed measurements

i = 2
plot(frap_time,
     frap_meas_norm[[1]],
     type = "l",
     main = paste("Overlay[FRAP(t)]:",
                  measurement_name),
     xlab = "Zeit [s]",
     ylab = "FRAP(t) [normalisiert]")

while(i < length(frap_meas_norm)){
  
  lines(frap_time,
        frap_meas_norm[[i]],
        add = TRUE)
  i = i+1
  
}



# Plotting the Data with Standard Error
ggplot(frap_mean,
       aes(x = frap_time,
           y = frap_meas_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = frap_meas_mean - 3*frap_meas_se,
                  ymax = frap_meas_mean + 3*frap_meas_se),
              alpha = 0.3)+ 
  ylab("FRAP(t) [normalisiert]") +
  xlab("Zeit [s]") +
  ggtitle(paste("Av[FRAP(t)]:",
                measurement_name,
                "(+ SE)"))

plot(frap_mean$frap_time,
     frap_mean$frap_meas_mean,
     pch = 1,
     cex = 0.8)
lines(frap_mean$frap_time,
      frap_mean$frap_meas_mean + frap_mean$frap_meas_se,
      lty = 2,
      add = TRUE)
lines(frap_mean$frap_time,
      frap_mean$frap_meas_mean - frap_mean$frap_meas_se,
      lty = 2,
      add = TRUE)



# Curve Fitting to the Data with determined weights
## Getting the data in to shape (t should start at 0)

fit_data = data.frame(frap_mean$frap_time[(postInterval):length(frap_mean$frap_time)] - frap_mean$frap_time[postInterval[1]], frap_mean$frap_meas_mean[postInterval:length(frap_mean$frap_time)], frap_mean$frap_meas_sd[(postInterval):length(frap_mean$frap_time)])


fit = nls(fit_data[[2]] ~ A*(1-exp(-( fit_data$frap_mean.frap_time..postInterval..length.frap_mean.frap_time..... / tau))),
          data = fit_data,
          start = list(A = 1, tau = 1),
          weights = (fit_data[[3]]))

fit_data$fitted = predict(fit, newdata = fit_data[[1]])



## See Additionial 1



ggplot(fit_data, aes(x = fit_data[[1]], y = fit_data[[2]])) + 
  geom_point() +
  geom_line(aes(x = fit_data[[1]], y = fit_data[[4]])) +
  ggtitle(paste("Fit[FRAP(t)]:",
                measurement_name)) +
  xlab("Zeit [s]") + 
  ylab("FRAP(t)")

ggplot(fit_data, aes(x = residuals(fit))) + 
  geom_line(stat = "density") +
  ggtitle(paste("Residuen[FRAP(t)]:",
                measurement_name)) +
  xlab("FRAP(t)") + 
  ylab("Dichte")
  
  




# Saving the obtained data in user-defined file path

fit_output = as.vector(summary(fit)[10])


write.csv(frap_mean,
          file = paste0(path,
                        measurement_name,
                        "_meanData.csv"))



write.csv(fit_output,
          file = paste0(path,
                        measurement_name,
                        "_fitData.csv"))









#Additional 1:
## Trying the Ellenberg model for fitting to diffusion coefficients (I*(1-(w^2*(w^2+4*pi*D*t)^(-1)))^(1/2))

# omega = 2
# 
# fit_alt = nls(fit_data$frap_mean.frap_meas_mean.postInterval.length.frap_mean.frap_time.. ~ I * (1 - (w^2 * (w^2 + 4 * pi * D * fit_data$frap_mean.frap_time..postInterval..length.frap_mean.frap_time.....)^(-1))^(1/2)), data = fit_data, start = list(I = 0.9, D = 2, w = 2))


## Using the drc package by Christian Ritz et al.
# fit = drm((frap_mean$frap_meas_mean) ~ frap_mean$frap_time,
#           fct = AR.3(names = c("cor","Mf", "tau")),
#           start = c(-25000, 0.85, 1),
#           subset = postInterval:length(frap_mean$frap_time)
#           #,weights = frap_mean$frap_meas_sd)
# )




