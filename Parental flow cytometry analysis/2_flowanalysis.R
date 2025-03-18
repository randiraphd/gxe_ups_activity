#Adapted by Randi Avery in 2024 from code by Mahlon Collins:
  #Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
  #Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570


#flow (.fcs) file names need to be in the following format:
#date_initials_plate_environment_strain_reporter_replicate.fcs
#controls should be named as follows:
#date_initials_plate_negctrlenvironment_strain_reporter_replicate.fcs (no fluorescent controls actually have no reporter, but the reporter name is included to know which day data was recorded - one reporter/day analyzed)


#base_dir=getwd()
base_dir <- "~/myproteasome"
proj  <- "/reporter_characterization/environments/combined"
proj_dir <- paste0(base_dir, proj)
setwd(proj_dir)


## load all required packages
source("~/myproteasome/reporter_characterization/load_flow_packages.R") 
## -----
## required for merging flowsets into a single flowframe
source(file = "https://raw.githubusercontent.com/mac230/flow_scripts/master/set2frame.R")
library(flowCore)
library(ggplot2)
library(stringr)
library(reshape2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)


#colors for plots
col_by <- "#2166ACFF" 
col_rm <- "#BF3232FF"
mymaroon = '#862334'
mygold = '#FFD966'

needed_dirs <- c("/fcs", "/results", "/tables",
                 "/dataframes", "/sessions",
                 "/dataframes/gated", "/dataframes/ungated", "/fcs/neg_ctrls")

dir_maker <- function(x){
    ifelse(!dir.exists(paths = paste0(proj_dir, x)),
           dir.create(path = paste0(proj_dir, x)),
           paste0("dir ", paste0(proj_dir, x), " exists."))
}
sapply(X = needed_dirs, FUN = dir_maker)

work_dir       <- paste0(proj_dir, "/fcs")
results_dir    <- paste0(proj_dir, "/results")
tables_dir     <- paste0(proj_dir, "/tables")
sessions_dir   <- paste0(proj_dir, "/sessions")
frame_dir      <- paste0(proj_dir, "/dataframes")
gated_dir      <- paste0(frame_dir, "/gated/")
ungated_dir    <- paste0(frame_dir, "/ungated/")
negctrls_dir    <- paste0(base_dir, "/fcs/neg_ctrls/") #not used

out_log        <- paste0(results_dir, "/2.6.25output_log")
write(x = '2.6.25 Start Here',
      file = out_log,
      append = T) 

print('put .fcs files in the newly made fcs folders')

## -----
## read in all the fcs files in a directory 

all_flow <- read.flowSet(files = NULL,
                         path = work_dir,
                         pattern = ".*.fcs",
                         alter.names = T,
                         min.limit = 2)

 
## get the number of cells per file; write output
cells_per_file <- fsApply(x = all_flow,
                          FUN = function(x) {
                              nrow(as.data.frame(exprs(x))) 
                          })

## write information string to the output log
read_string <- paste0("1. On ",
                      Sys.time(),
                      " files were read in via 'read.flowSet'.\n\nCells per file:\n")

write(x = read_string,
      file = out_log,
      append = T)
#this file saves in the 'results' folder

## for nicely aligned output in the resulting file
write(capture.output(cells_per_file),
    file = out_log,
    append = T,
    sep = "\n")

all_filtered <- all_flow

## -----
## now convert flowframes to dataframes and merge to a single frame
## start w/ ungated cells
ungated_frames <- fsApply(x = all_filtered,
                          FUN = function(x) {
                              
                              ## extract file name; we'll
                              ## add this to the dataframe
                              file_name <- rep(x = x@description$GUID,
                                               times = nrow(as.data.frame(exprs(x))))
                              
                              ## extract the time the sample was run;
                              ## then convert to a numeric value
                              ## BTIM = "beginning time"
                              file_time <- rep(x = x@description$`$BTIM`,
                                               times = nrow(as.data.frame(exprs(x))))
                              time_conv <- HmsToSec(file_time)
                              
                              ## combine flow data and file name
                              ## into a single data frame 
                              cbind(as.data.frame(exprs(x)),
                                    file_name,
                                    file_time,
                                    time_conv)
})

ungated_final <- do.call("rbind", ungated_frames)


## -----
## log fluorescence values and compute TFT ratio
ungated_final$log_GFP   <- log10(x = ungated_final$GFP.A) 
ungated_final$log_RFP   <- log10(x = ungated_final$mCherry.A)


ungated_final$TFT_ratio <- -log(x = ungated_final$mCherry.A / ungated_final$GFP.A,
                                base = 2)


## -----
## now gate the cells to capture the haploid cell population.
## we identify haploids as a sharp peak in the lower end of
## the fsc density plot.  I take 10% above and below the
## max density value.
ungated_final$file_string <- as.character(ungated_final$file_name)

## ungated = all cells
ungated_final$ungated <- rep(T, nrow(ungated_final)) #adds a column that says 'TRUE'


#run the following code (building 'densities' df) and make sure the fsc gating makes sense. There is an upper limit to reading in the data, so some samples produce a peak that is not accurate - based on the upper limit value being used.
#if there are more samples that what RStudio shows, separate data, and view all fsc plots
  #savedata = ungated_final
  #ungated_final = savedata
  #ungated_final$plate = str_split_i((ungated_final$file_string), '_', 3) 
  #ungated_final = ungated_final[ungated_final$plate == "plate02", ]

#6/10/24 - running the 'densities' code below showed that many Low Nitrogen samples produced double peaks.
#as outlined in the Methods, the lower peaks were were used for gating, so those samples will be gated separately from the rest of the data, and combined later


ungated_final$split = str_split_i((ungated_final$file_string), '_', 4)
#ungated_finalsave = ungated_final
ungated_finallowN = ungated_final[ungated_final$split == "LowN", ] #used for lower peak
ungated_finalrest = ungated_final[ungated_final$split != "LowN", ] #used for rest of the data


#used for gating LowN samples - using lower peak
densitieslowN <- lapply(X = unique(ungated_finallowN$file_string), 
                    FUN = function(x) {

                      #to run one at a time
                      #x=X[1]

                      ## file name for dataframe creation
                      f_name <- x
                      #f_name = "20240126_RRA_plate02_BTZ_RM_UFD_008.fcs" #example of one file

                      ## per file dataframe for calculating density
                      dat  <- ungated_finallowN[ungated_finallowN$file_name == f_name, ]

                      ## get the density for FSC
                      dens <- density(dat$FSC.A)


                      #makes an fsc plot of all the samples.

                      ## 10% +/- the max FSC density peak
                      #dens_max <- dens$x[which.max(dens$y)]
                      checkmax = dens$x[which.max(dens$y)]
                      if(checkmax > 110000){
                        dfx = dens$x
                        dfy = dens$y
                        df = cbind.data.frame(dfx, dfy)
                        cutdf = df[df$dfx < 110000, ] #picked this value visually based on the plot
                        newmax = max(cutdf$dfy)
                        secondmax = cutdf$dfx[cutdf$dfy == newmax]
                        newmaxindex = which(dens$x == secondmax)
                        dens_max = dens$x[newmaxindex]} else {
                          dens_max = dens$x[which.max(dens$y)]
                        }

                      dens_up <- (0.1 * dens_max) + dens_max
                      dens_down <- dens_max - (0.1 * dens_max)
                      data.frame(f_name, dens_up, dens_down, stringsAsFactors = F)
                      #print(dens) #gives summary
                      #uncomment following 3 lines to make plots. Then comment, and run for loop again to make df
                      #plot(dens, main = f_name)
                      #abline(v=dens_up)
                      #abline(v=dens_down)
                    })

densitieslowN <- do.call("rbind", densitieslowN)


## gated = all cells w/in 10% +/- FSC max density
densities <- lapply(X = unique(ungated_finalrest$file_string),
                    FUN = function(x) {
                        
                      #to run one at a time
                      #x=X[1]
                      
                        ## file name for dataframe creation
                        f_name <- x
                        #f_name = "20240126_RRA_plate02_BTZ_RM_UFD_008.fcs"

                        ## per file dataframe for calculating density
                        dat  <- ungated_finalrest[ungated_finalrest$file_name == f_name, ]
                        
                        ## get the density for FSC
                        dens <- density(dat$FSC.A)
                        
                        
                        ## 10% +/- the max FSC density peak
                        #dens_max <- dens$x[which.max(dens$y)]
                        checkmax = dens$x[which.max(dens$y)]
                        if(checkmax > 225000){  #put limit so maxed out value isn't picked for peak
                          dfx = dens$x
                          dfy = dens$y
                          df = cbind.data.frame(dfx, dfy)
                          cutdf = df[df$dfx < 225000, ]
                          newmax = max(cutdf$dfy)
                          secondmax = cutdf$dfx[cutdf$dfy == newmax]
                          newmaxindex = which(dens$x == secondmax)
                          dens_max = dens$x[newmaxindex]} else {
                            dens_max = dens$x[which.max(dens$y)]
                          }
                          
                        dens_up <- (0.1 * dens_max) + dens_max
                        dens_down <- dens_max - (0.1 * dens_max)
                        data.frame(f_name, dens_up, dens_down, stringsAsFactors = F)
                        #print(dens) #gives summary
                        #uncomment following 3 lines to make plots. Then comment, and run for loop again to make df
                        #plot(dens, main = f_name)
                        #abline(v=dens_up)
                        #abline(v=dens_down)
                    })


## bind each file's density estimate into a single dataframe
densities <- do.call("rbind", densities)

densities = rbind(densities, densitieslowN) #combine with the low nitrogen samples

#adds a column in ungated_final saying if that datapoint is within the FSC peak
ungated_final$gated <- sapply(X = 1:nrow(ungated_final), #future idea: add a line to print file name since this command takes very long, especially with 800 samples!
                              FUN = function(x) {
                                  
                                  ## get file name for current row
                                  f_name <- ungated_final[x, "file_string"]

                                  ## get corresponding row of density estimates
                                  dens_row <- densities[densities$f_name == f_name, ]
                                  
                                  ## 10% +/- the max FSC density peak
                                  dens_up   <- dens_row$dens_up
                                  dens_down <- dens_row$dens_down
                                  
                                  ## create a gate and subset to keep cells in range
                                  ifelse(ungated_final[x, "FSC.A"] > dens_down &
                                         ungated_final[x, "FSC.A"] < dens_up,
                                         yes = T, no = F)

                     })


ungated_per_file <- tapply(X = ungated_final$FSC.A,
                           INDEX = ungated_final$file_name,
                           FUN = length)

#shows how many cells are saved in each file
gated_per_file <- tapply(X = ungated_final$FSC.A[ungated_final$gated == T],
                         INDEX = ungated_final$file_name[ungated_final$gated == T],
                         FUN = length)

fsc_gated_counts <- data.frame(kept = gated_per_file,
                               excluded = ungated_per_file - gated_per_file,
                               percent = gated_per_file / ungated_per_file)

gate_string <- paste0("\n3. On ", Sys.time(),
                      " Filtered cells to grab 10% +/- the central FSC peak\n\n",
                      "The following counts were obtained:\n")

write(x = gate_string,
      file = out_log,
      append = T)

cat(capture.output(fsc_gated_counts),
    file = out_log,
    append = T,
    sep = "\n")


## -----
## write out summary statistics for gated and ungated files
## extract only numeric parameters for summary statistics 
f_params <- colnames(ungated_final)
f_params <- f_params[sapply(X = f_params,
           FUN = function(x) {
               is.numeric(ungated_final[, x])
               })]

## -----
## write out summary statistics for all parameters 
## make sure all columns are printed together 
options(width = 200)
#only prints values for first 100 samples
for(p in 1:length(f_params)) {

    gated_final    <- ungated_final[ungated_final$gated == T, ]

    ungated_mean   <- tapply(X = ungated_final[, f_params[p]],
                             INDEX = ungated_final$file_name,
                             FUN = mean)
    gated_mean     <- tapply(X = gated_final[, f_params[p]],
                             INDEX = gated_final$file_name,
                             FUN = mean)
    mean_frac      <- round(x = gated_mean / ungated_mean,
                            digits = 2)
    ungated_median <- tapply(X = ungated_final[, f_params[p]],
                             INDEX = ungated_final$file_name,
                             FUN = median)
    gated_median   <- tapply(X = gated_final[, f_params[p]],
                             INDEX = gated_final$file_name,
                             FUN = median)
    median_frac      <- round(x = gated_median / ungated_median,
                              digits = 2)
    ungated_sd     <- tapply(X = ungated_final[, f_params[p]],
                             INDEX = ungated_final$file_name,
                             FUN = sd)
    gated_sd       <- tapply(X = gated_final[, f_params[p]],
                             INDEX = gated_final$file_name,
                             FUN = sd)
        
    ungated_cv     <- ungated_mean / ungated_sd

    gated_cv       <- gated_mean / gated_sd

    sum_stats      <- data.frame(ungated_mean, gated_mean, mean_frac,
                                 ungated_median, gated_median, median_frac,
                                 ungated_sd, gated_sd,
                                 ungated_cv, gated_cv)

    sum_string     <- paste0("\nOn ", Sys.time(), " Obtained summary statistics for parameter: ",
                             f_params[p], "\n\n")

    write(x = sum_string,
      file = out_log,
      append = T)

    cat(capture.output(sum_stats),
        file = out_log,
        append = T,
        sep = "\n")
    
}


## -----
## extract strain, reporter, and replicate as factors
ungated_final$strain <- sapply(X = 1:length(ungated_final$file_string),
                               FUN = function(x){
                                   strsplit(x = ungated_final$file_string[x], split = "_")[[1]][5]
                               })
ungated_final$strain_factor <- as.factor(ungated_final$strain)

## now write the unique values we obtain for strain in ea. dataset
write(x = sprintf("%s", c("\nUnique values for 'strain' variable:\n",
                          unique(ungated_final$strain))),
      file = out_log,
      append = T)

## now for reporter
ungated_final$reporter <- sapply(X = 1:length(ungated_final$file_string),
                                 FUN = function(x){
                                     strsplit(x = ungated_final$file_string[x], split = "_")[[1]][6] #reporter named in 6th element
                                 })
ungated_final$reporter_factor <- as.factor(ungated_final$reporter)

## log the results
write(x = sprintf("%s", c("\nUnique values for 'reporter' variable for ungated set:\n",
                          unique(ungated_final$reporter))),
      file = out_log,
      append = T)

## now for replicates
ungated_final$replicate <- sapply(X = 1:length(ungated_final$file_string),
                                  FUN = function(x){
                                      as.numeric( 
                                          gsub(pattern = ".fcs",
                                               replacement = "",
                                               x = strsplit(x = ungated_final$file_string[x],
                                                            split = "_")[[1]][7]))
                                  })
## for later use in plotting individual replicates
ungated_final$replicate_factor <- as.factor(ungated_final$replicate)


## log the results
write(x = sprintf("%s", c("\nUnique values for 'replicate' variable for ungated set:\n",
                          unique(ungated_final$replicate))),
      file = out_log,
      append = T)

ungated_final$environment = str_split_i((ungated_final$file_string), '_', 4)

write(x = sprintf("%s", c("\nUnique values for 'environment' variable for ungated set:\n",
                          unique(ungated_final$environment))),
      file = out_log,
      append = T)


## -----
## get plate information and create a variable that 
## expresses time relative to the first sample for ea. plate

ungated_final$plate = str_split_i((ungated_final$file_string), '_', 3)

## log relative time values per file
time_string <- paste0("\nOn ", Sys.time(), "Extracted time values for each file/plate\n")

write(x = time_string,
      file = out_log,
      append = T) 

time_frame <- lapply(X = unique(ungated_final$file_string),
                     FUN = function(x) {

                         sub_set <- ungated_final$file_string == x
                         f_params <- c("file_string", "strain", "plate",
                                     "file_time", "time_conv", 'reporter')
                         dat <- ungated_final[sub_set, f_params]
                         dat_out <- data.frame(file = unique(dat$file_string),
                                               strain = unique(dat$strain),
                                               plate = unique(dat$plate),
                                               file_time = unique(dat$file_time),
                                               time_conv = unique(dat$time_conv),
                                               reporter = unique(dat$reporter))
                                               })

time_frame <- do.call("rbind", time_frame)

## convert time to a relative value based on the first sample of ea. plate

#only works if there are two plates - comment out indicated plates if only using one plate
time_mins = lapply(X = unique(time_frame$reporter),
                      FUN = function(x){
                        sub_set <- time_frame$reporter == x 
                        f_params <- c("file","time_conv", 'plate')
                        dat = time_frame[sub_set, f_params]
                        plate_one_min <- min(dat$time_conv[dat$plate == "plate01"])
                        plate_two_min <- min(dat$time_conv[dat$plate == "plate02"])#comment out
                        plate1df = dat[dat$plate == 'plate01',]
                        plate2df = dat[dat$plate == 'plate02',]#comment out
                        plate1df$min_time = plate_one_min
                        plate2df$min_time = plate_two_min #comment out
                        newdf = rbind(plate1df, plate2df) #comment out
                        #newdf = plate1df #comment out if two plates
                      })
time_mins <- do.call("rbind", time_mins)
newtimeframe = merge(time_frame,time_mins, by = c('file','plate','time_conv')) #combine the df with the min values
newtimeframe$time_rel = newtimeframe$time_conv - newtimeframe$min_time
time_frame = newtimeframe


cat(capture.output(time_frame),
    file = out_log,
    append = T,
    sep = "\n")


## -----
## save the final dataframe for future use:
save(ungated_final,
     file = paste0(frame_dir, "/6.10.24ungated_final.R"))

write(x = paste0("\nOn ", Sys.time(), " saved final dataframe as: ",
                 paste0(frame_dir, "/", "2.6.25ungated_final.R"),
                 "\n"),
                 file = out_log,
                 append = T)

save(work_dir, results_dir, tables_dir, sessions_dir,
     frame_dir, gated_dir, ungated_dir, out_log,
     file = paste0(frame_dir, "/dir_structure.R"))

write(x = paste0("\nOn ", Sys.time(), " saved dir structure as: ",
                 paste0(frame_dir, "/dir_structure.R"),
                 "\n"),
                 file = out_log,
                 append = T)

## -----
## re-load the data by uncommenting the following:
## load(file = paste0(frame_dir, "/", "ungated_final.R"))
## load(file = paste0(frame_dir, "/dir_structure.R"))

## -----
## extract the median from each biological replicate and use this
## value to build a dataframe w/ n. replicate observations per strain
## per reporter.  This dataframe is what we'll use for stats and for
## creating stripcharts, boxplots, and heatmaps
## 'aggregate' creates a new dataframe from x by applying FUN to
## all unique combinations of the factors supplied to the 'by'
## argument - in this case, grab the mean of numeric data and
## keep everything else a factor

ungated_medians <- aggregate.data.frame(x = ungated_final,
                                        by = list(ungated_final$strain_factor,
                                                  ungated_final$replicate_factor,
                                                  ungated_final$reporter_factor,
                                                  ungated_final$environment, 
                                                  ungated_final$plate,
                                                  ungated_final$reporter),
                                        FUN = function(x) {
                                            ifelse(is.numeric(x), median(x), as.character(x))
                                        },
                                        ## simplify results to vector 
                                        simplify = T)


#adding time_rel by merging these df's
time_frame$file_string = time_frame$file
newdf = merge(time_frame, ungated_medians, by = c('file_string','plate','time_conv', 'strain', 'reporter', 'file_time'))
ungated_medians = newdf

## -----
## 'aggregate' seems to strip the levels from factors, so add
## these back using the values present in the original dataframe
ungated_medians$strain_factor    <- as.factor(ungated_medians$strain)
ungated_medians$replicate_factor <- as.factor(ungated_medians$replicate)
ungated_medians$reporter_factor  <- as.factor(ungated_medians$reporter)
ungated_medians$environment_factor <- as.factor(ungated_medians$environment)
ungated_medians$gating           <- rep("ungated", nrow(ungated_medians))

## -----
## adjust for the effect of time on the TFT ratio,
## then convert it to a Z-score: 
#essentially fits data to a smoothened line, makes mean be zero
ungated_loess <- loess(formula = TFT_ratio ~ time_rel,
                       data = ungated_medians)
#bump back to same mean
ungated_medians$TFT_loess <- ungated_loess$residuals + mean(ungated_medians$TFT_ratio)

## scale to 'wild-type' median with sd = 1. 
by_median <- median(ungated_medians$TFT_loess[ungated_medians$strain == "BY"]) 
ungated_medians$TFT_scaled <- scale(x = ungated_medians$TFT_loess,
                                    center = by_median,
                                    scale = T)

## -----

## drop the deletion strain level from the dataset
ungated_medians$strain_final <- droplevels(x = ungated_medians$strain_factor)
## levels(ungated_medians$strain_final)

## -----
## have to extract medians from the gated data separately
gated_final <- ungated_final[ungated_final$gated == T, ]

gated_medians <- aggregate.data.frame(x = gated_final,
                                       by = list(gated_final$strain_factor,
                                                 gated_final$replicate_factor,
                                                 gated_final$reporter_factor,
                                                 gated_final$environment,
                                                 gated_final$plate, #added to keep plates separate - SC is on plate 1 and 2
                                                 gated_final$reporter), 
                                      FUN = function(x) {
                                           ifelse(is.numeric(x), median(x), as.character(x))
                                       },
                                       ## simplify results to vector 
                                       simplify = T)

## 'aggregate' seems to strip the levels from factors, so add
## these back using the values present in the original dataframe
gated_medians$strain_factor    <- as.factor(gated_medians$strain)
gated_medians$replicate_factor <- as.factor(gated_medians$replicate)
gated_medians$reporter_factor  <- as.factor(gated_medians$reporter)
gated_medians$environment_factor <- as.factor(gated_medians$environment)
gated_medians$gating           <- rep("gated", nrow(gated_medians))


#combine to get time data
newdf2 = merge(time_frame, gated_medians, by = c('file_string','plate','time_conv', 'strain', 'reporter', 'file_time'))
gated_medians = newdf2


## -----
## adjust for the effect of time on the TFT ratio,
## then convert it to a Z-score: 
#TFT_loess is used in the final analysis
gated_loess <- loess(formula = TFT_ratio ~ time_rel,
                     data = gated_medians)
gated_medians$TFT_loess <- gated_loess$residuals + mean(gated_medians$TFT_ratio)

## scale to 'wild-type' median with sd = 1
by_median <- median(gated_medians$TFT_loess[gated_medians$strain == "BY"])
gated_medians$TFT_scaled <- scale(x = gated_medians$TFT_loess,
                                  center = by_median,
                                  scale = T)

gated_medians$strain_final = gated_medians$strain

####------

#now analyze negative controls to see if any samples don't pass the GFP fluorescence of the neg controls. If they don't pass, throw those out.
#neg ctrls are strains with no TFTs
negSC = gated_medians[gated_medians$environment == 'negctrlSC', ]
avgneg = mean(negSC$log_GFP)

#plot neg ctrls
ggplot(negSC, aes(x=reporter,  y=log_GFP, fill = strain)) +
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "Reporter", y = "log_GFP") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6, alpha = 0.4)+
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18)) +
  ggtitle("neg controls") +
  scale_fill_manual(values = c(col_by, col_rm))

#just want to see which samples might be under any neg ctrl value
maxneg = max(negSC$log_GFP)
#samples in this dataframe (that aren't negative controls) were thrown out
belowneg = gated_medians[gated_medians$log_GFP < maxneg, ]

#making x axis
belowneg$grp2 = paste0(belowneg$reporter, ' in ', belowneg$environment)

#this plot is what I used to show I threw out 4xUb and UFD in Low N
ggplot(belowneg, aes(x=grp2,  y=log_GFP, fill = strain)) +
  geom_dotplot(binaxis='y', stackdir='center') + 
  labs( x = "grp", y = "log_GFP") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("neg controls")

#only removed samples that didn't pass neg ctrl threshold for gated medians. 
gated_medians$grp2 = paste0(gated_medians$environment, '_', gated_medians$reporter)
filtered_gated = gated_medians[gated_medians$grp2 != "LowN_UFD",]
filtered_gated = filtered_gated[filtered_gated$grp2 != "LowN_4xUb",]
gated_medians = filtered_gated

#save data
all_data <- list()
all_data[[1]] <- ungated_medians
all_data[[2]] <- gated_medians
names(all_data) <- c("ungated_medians", "gated_medians")
save(all_data,
     file = paste0(frame_dir, "/2.7.25gating_all_data.R"))


#includes both sets of SC samples, when SC was also ran on plate02
write.table(all_data$gated_medians,
            file = "filtered_gated_medians2.7.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')
write.table(all_data$ungated_medians,
            file = "ungated_medians2.7.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


# #As indicated above in the gating code, Low Nitrogen produced two FSC peaks.
# #In order to decide to go with the lower FSC peak, I duplicated all Low Nitrogen files, then ran above code on one set of Low N files
# #gating the lower peak and the other set of Low N files gating the higher peak. Then I compared the -log2(RFP/GFP) to see if there
# #was a significant difference. 
# 
# #plotted time each sample was taken
# ggplot(all_data$gated_mediansLowN, aes(x=time_rel,  y=log_GFP, color = environment)) + 
#   geom_point(aes(shape=reporter, alpha = strain)) +
#   #scale_shape_manual(values=c(0, 16))+
#   geom_text(label=all_data$gated_mediansLowN$replicate_factor, vjust=-1.1)+
#   scale_x_continuous(breaks=seq(0,3000,500))+
#   theme(axis.text.x = element_text(size =10))
# #the LowN and LowNduplicated sets have the exact same timepoints, which they should
# 
# #now plotting the -log2(RFP/GFP) of the different reporters for the set gating the low peak and the set gating the high peak
# gated_mediansLowN$grp3 = paste0(gated_mediansLowN$reporter, '_', gated_medianLowNs$environment)
# ggplot(gated_mediansLowN, aes(x=grp3,  y=TFT_ratio, fill = strain)) +
#   geom_boxplot(outlier.size = 0.3, position = 'identity', alpha = 0.7) + 
#   labs( x = "Reporter and Environment", y = "-log2(RFP/GFP)") + #-log2(RFP/GFP)
#   geom_dotplot(binaxis='y', dotsize= 5 , alpha = 0.4, binwidth = 0.01, stackdir='center')+ 
#   #geom_jitter(color= 'pink', size=2, alpha=0.9) +
#   #geom_point(position = position_jitter(w = 0.1, h = 0)) +
#   theme(axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 10)) +
#   ggtitle("Different LowN Gating. LowNdup = smaller fsc peak") + 
#   scale_fill_manual(values = c(col_by, col_rm))
# 
# #t-tests to see if gating made a significant difference between -log2(RFP/GFP) of gating of low vs high peak.
# #Instead of a for loop, I just typed in the diff reporters and ran the following lines.
# #dependent variable was -log2(RFP/GFP). Performed t-test separately for BY and RM strains for each reporter.
# mydf = gated_mediansLowN[gated_mediansLowN$reporter == 'Thr', ] #re-ran this line for each reporter
# mydfby = mydf[mydf$strain == 'BY', ]
# mydfrm = mydf[mydf$strain == 'RM', ]
# t.test(TFT_ratio ~ environment, data = mydfby)
# t.test(TFT_ratio ~ environment, data = mydfrm)
# 
# #the above analysis, plus the distribution of densities for the fsc peaks, showed that the lower peak was a more accurate fsc value, and what I used moving forward for analysis.
