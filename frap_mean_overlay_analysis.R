rm(list = ls(all = TRUE))

# Functions for Overlay Analysis
## The following script is designed to analyse Fluorescence Recovery after Photobleaching experiments in "batch" mode. This means, that you as the user can define the path on the file system of your computer to the measurement data in .csv-format. Preferably, you stored your data in seperate folders for each individual experiment. The script initially loops through the folders of the file path you give to it and loads all of your data in a data frame called "frap". Then you need to define the rows in your data you want to evaluate, which row represents the background and which the correction ROI of observational photobleaching.
## In the variable "postIntervall" you need to tell the script the index of your measurement data, where your post bleaching intervall startet and ended. 
## In the end all of the evaluated data are stored in .csv-file (time, mean, standard deviation, standard error)






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

measurement_name = "S100A11_wBLM_Ncl"


require(readr)
require(ggplot2)
require(drc)

####### START
path = "/Users/tomkache/Documents/Studium/Biologie/SS 2017 B. Sc. Biology FSU Jena/Bachelorarbeit/FRAP-Experimente/S100A11-EGFP/S100A11-EGFP_wBLM_29_07_17/"
file_name = "/test.csv"

colMeas = 2
colBack = 4
colFade = 5

postIntveral = 103:202

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

# Seperate the ROIs of the nucleus, the cytoplasm, the background and the observational photobleaching

frap_time = frap[[1]]

i = colMeas
frap_ncl = data.frame(frap[[i]])

while(i <= length(frap)){
  
  dat = frap[[i]]
  frap_ncl = cbind(frap_ncl, dat)
  
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

frap_ncl_cb = frap_ncl - frap_back
frap_fade_cb = frap_fade - frap_back

## Removing Observational Photobleaching by division (measurement / fading)

frap_ncl_cor = frap_ncl_cb / frap_fade_cb


## Normalising
frap_ncl_norm = as.data.frame(frap.normalize(frap_ncl_cor[[1]]))
i = 2
while(i <= length(frap_ncl_cor)){
  
  dat = frap.normalize(frap_ncl_cor[[i]])
  frap_ncl_norm = cbind(frap_ncl_norm, dat)
  i = i + 1
  
}


# Calculating the mean, standard deviaton and standard error to obtain the "average" curve
## We will loop through each row of the normalized data frame and calculate these values

frap_ncl_mean = rowMeans(frap_ncl_norm)
frap_ncl_sd = apply(frap_ncl_norm,
                    1,
                    sd)
sample_size = length(frap_ncl_norm)

frap_ncl_se = apply(as.matrix(frap_ncl_sd),
                    1:2,
                    function(x) x/sqrt(sample_size))


# Shifting the Curve: Time begins at Bleaching Event

#frap_time = frap_time - frap_time[which.min(frap_ncl_mean)]

# Getting all the data in one data frame (time, mean of measurement, sd, se)

frap_mean = data.frame(frap_time,
                       frap_ncl_mean,
                       frap_ncl_sd,
                       frap_ncl_se)


# Overlayed measurements

i = 2
plot(frap_time,
     frap_ncl_norm[[1]],
     type = "l",
     main = paste("Overlay[FRAP(t)]:",
                  measurement_name),
     xlab = "Zeit [s]",
     ylab = "FRAP(t) [normalisiert]")

while(i < length(frap_ncl_norm)){
  
  lines(frap_time,
        frap_ncl_norm[[i]],
        add = TRUE)
  i = i+1
  
}



# Plotting the Data with Standard Error
ggplot(frap_mean,
       aes(x = frap_time,
           y = frap_ncl_mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = frap_ncl_mean - 3*frap_ncl_se,
                  ymax = frap_ncl_mean + 3*frap_ncl_se),
              alpha = 0.3)+ ylab("FRAP(t) [normalisiert]") +
  xlab("Zeit [s]") +
  ggtitle(paste("Av[FRAP(t)]:",
                measurement_name,
                "(+ SE)"))

plot(frap_mean$frap_time,
     frap_mean$frap_ncl_mean,
     pch = 1,
     cex = 0.8)
lines(frap_mean$frap_time,
      frap_mean$frap_ncl_mean + frap_mean$frap_ncl_se,
      lty = 2,
      add = TRUE)
lines(frap_mean$frap_time,
      frap_mean$frap_ncl_mean - frap_mean$frap_ncl_se,
      lty = 2,
      add = TRUE)



# Curve Fitting to the Data with determined weights

#nls(egfp_mean$egfp_ncl_mean ~ mf*(1-exp(-egfp_mean$egfp_time/tau)),start = list(mf = 0.9, tau = 1), algorithm = "port")


# fit = drm((frap_mean$frap_ncl_mean) ~ frap_mean$frap_time,
#           fct = AR.3(names = c("cor","Mf", "tau")),
#           start = c(-25000, 0.85, 1),
#           subset = postIntveral
#           #,weights = frap_mean$frap_ncl_sd)
# )




plot(frap_mean$frap_time,
     frap_mean$frap_ncl_mean,
     main = paste("Av[FRAP(t)]:",
                  measurement_name,
                  "(fitted)"),
     xlab = "Zeit [s]",
     ylab = "FRAP(t) [normalisiert]")

plot(fit,
     add = TRUE)


fit_output = as.vector(summary(fit)[3])


write.csv(frap_mean,
          file = paste0(path,
                        measurement_name,
                        "_meanData.csv"))



write.csv(fit_output,
          file = paste0(path,
                        measurement_name,
                        "_fitData.csv"))







#test = nls(frap_mean$frap_ncl_mean[104:202] ~ A*(1-exp(-frap_mean$frap_time[104:202] / tau)), start = list(A = 0.8, tau = 1))
