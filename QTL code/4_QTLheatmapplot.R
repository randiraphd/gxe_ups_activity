#Adapted by Randi Avery in 2024 from code by Mahlon Collins:
  #Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
  #Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570

#Analyzes and describes the QTL data. Organizes data for further plots made in "5_analyzingQTL_Plots.R"
#Makes heatmap for Figure 3A and Supplementary Figure 4A

base_dir <- "~/myproteasome/"
proj  <- "2024.07.05_allHL_GxE/"
proj_dir <- paste0(base_dir, proj)
peak_dir <- paste0(base_dir, proj, "peaks/")
setwd(proj_dir)


#read in QTL data from previous script: ~/myproteasome/paper_final_QTLs_24.07.12.R
qtls = read.table(file = "QTLs_averaged_keptsinglereps24.07.12.txt", header = T)

#change names of reporters and environments to how I want them in the plot:
qtls[qtls == "rpn4_degron_redo"] = "Rpn4"
qtls[qtls == "rpn4_degron"] = "Rpn4"
qtls[qtls == "4x_Ub"] = "4xUb"
qtls[qtls == "Asn_N-end"] = "Asn"
qtls[qtls == "Phe_N-end"] = "Phe"
qtls[qtls == "Thr_N-end"] = "Thr"
qtls[qtls == "Low_Glucose"] = "Low G" #no underscore is ok here
qtls[qtls == "Low_Nitrogen"] = "Low N"#no underscore is ok here
qtls[qtls == "Bortezomib"] = "BTZ"


#### Heatmap

#the heatmap should only have QTLs that are present in both replicates!
qtls.hm = qtls[qtls$number_of_reps == 2,] #should be 416 QTLs


#load packages and functions
library("lattice")
library("latticeExtra")
library("RColorBrewer")
library("grid")
source(paste0(proj_dir, "QTL_scripts/gTest.R"))
source(paste0(proj_dir, "QTL_scripts/x_qtl_seq_functions_170831.R"))
source(paste0(proj_dir, "QTL_scripts/mp_JB_170901.R"))
source(paste0(proj_dir, "QTL_scripts/peaksFromVector.R"))

sepBetweenChr <- 0
trimFromEnd   <- 15e3
obsMin        <- 10
LoessSpan     <- 0.1
## same as in Albert 2014
AF_thres       <- 0.09653124
## multipool LOD threshold, our usual value
multi_thres    <- 4.5


## -----
## <<heatmap_setup>>
## build a dummy dataframe to drop QTLs into
## ~120 bins if you use 100 kb windows for the map

## read in chr lengths
## need this for building the heatmap
chr_lengths <- read.table(paste0(proj_dir, "QTL_scripts/sacCer3ChromLengths.txt"), header = F)

## don't need mitochondrial chromosome 
chr_lengths <- chr_lengths[1:16, ]

## add the gcoords and numeric (not roman) chromosomes:
chr_lengths$chr <- 1:16
chr_lengths$gcoords <- getGcoords(chr = 1:16,
                                  pos = chr_lengths$chr,
                                  spacing = sepBetweenChr)
chr_lengths <- chr_lengths[, c(3, 1, 2, 4)]
names(chr_lengths) <- c("chr", "chr_r", "length", "gcoords")


#Combine reporter and environment. Can do both directions to sort by either environment or reporter
qtls.hm$rep_env = paste0(qtls.hm$reporter, ' in ', qtls.hm$environment) #this was actually already in the table as 'rep'
qtls.hm$env_rep = paste0(qtls.hm$environment, ' | ', qtls.hm$reporter)
#save df's in order to replace 'reporter' column
qtls_saved = qtls.hm
#qtls.hm = qtls_saved

####then choose if i want reporter or environment first####
#sorted by reporter
qtls.hm$reporter.hm = qtls.hm$rep_env
#added 'unique' on 10/18/24
reporters <- unique(qtls.hm$rep_env)
#sorted by environment
#qtls.hm$reporter.hm = qtls.hm$env_rep
#reporters <- qtls.hm$env_rep


##df has to have these columns in this order:
# chr_bins <- data.frame(reporter.hm = rep(reporters[k], length(bins)),
#                        chr = chr_indices,
#                        LOD = rep(0, length(bins)),
#                        delta_AF = rep(0, length(bins)),
#                        left_Index = rep(0, length(bins)),
#                        max_Index = rep(0, length(bins)),
#                        right_Index = rep(0, length(bins)),
#                        bin = bins)

keeps = c('reporter.hm', 'chr', 'LOD', 'delta_AF', 'left_Index', 'max_Index', 'right_Index')
qtls.hm = qtls.hm[keeps]

#add rows of 'known loci' so i can see what bins they go to
#knownloci = read.table(file = "knownlociQTLformat.txt", header = TRUE, sep = '\t')
#knownloci$reporter.hm = paste0(knownloci$reporter.hm, ' in ', knownloci$reporter.hm) #so then later when I split the sample name, it is just there twice
#qtls.hm2 = rbind(qtls.hm, knownloci)
#reporters <- unique(qtls.hm2$reporter.hm)

## <<build_heatmap_dataframe>>
#k <- 1

final <- list()
#here, 'reporters' is the 47 samples
for (k in 1:length(reporters)) {    
  
  ## make a list of 16 containing the chr and bin as a dataframe
  ## have to round up the chr lengths so the bins are all the same
  ## size (1e5 bp)
  
  ## we need to make an extra bin for ea. chromosome, since the
  ## lengths don't end in 1e5 multiples.  'round_up' allows you
  ## to do this using a multiple of your choice.
  ## So, chrI = 230218 bp and round_up(230218, 1e5) = 3e5
  round_up <- function(from, to) {
    ceiling(from / to) * to
  }
  
  round_down <- function(from, to) {
    floor(from / to) * to
  }    
  
  ## assign each bin to the correct chromosome
  chr_indices <- unlist(sapply(X = 1:16, FUN = function(x) {
    rep(x,
        times = length(seq(from = 1e5,
                           to = round_up(chr_lengths$length[x],
                                         to = 1e5),
                           by = 1e5)))
  }))
  
  ## create chromosome bins - this will serve as a conditioning factor
  bins <- seq(from = 1e5,
              to = length(chr_indices) * 1e5,
              by = 1e5)
  
  ## when we assign peaks to bins, we need to add the length of
  ## all preceding chromosomes plus the peak position.  we'll use
  ## the bin lengths, not the chromosome lengths for this purpose.
  ## so, a QTL on chrII at base 150 = 3e5 (length of chr I bins) + 150
  bin_sums <- c(0, unlist(sapply(X = 2:16, FUN = function(x) {
    max(bins[chr_indices == x - 1])
  })))
  
  ## assign the peaks for ea. reporter to a dataframe 
  out <- qtls.hm[qtls.hm$reporter.hm == reporters[k], ]
  
  ## assign each peak a bin for the heatmap 
  for (b in 1:nrow(out)) {
    out$bin[b] <- round_down(from = out$max_Index[b] + bin_sums[out$chr[b]],
                             to = 1e5)
  }
  
  ## dummy dataframe that contains all bins    
  chr_bins <- data.frame(reporter.hm = rep(reporters[k], length(bins)),
                         chr = chr_indices,
                         LOD = rep(0, length(bins)),
                         delta_AF = rep(0, length(bins)),
                         left_Index = rep(0, length(bins)),
                         max_Index = rep(0, length(bins)),
                         right_Index = rep(0, length(bins)),
                         bin = bins)
  
  ## the bins above aren't exact, so assign them to one
  ## of the bins we defined using 'which.min'
  ## assign peaks to bins via 'which.min'
  for (r in 1:nrow(out)) {
    ind <- out$chr[r]
    chrs <- chr_bins$bin[chr_bins$chr == ind & chr_bins$bin > out$bin[r]]
    out$bin[r] <- chrs[which.min(abs(out$bin[r] - chrs))]
  }
  
  ## if there's a peak, replace the current row with the peak row
  ## if not, just leave everything at 0
  for (n in 1:nrow(chr_bins)) {
    test <- min(abs(chr_bins[n, "bin"] - out$bin))
    ind  <- which.min(abs(chr_bins[n, "bin"] - out$bin))
    chr_bins[n, ] <- if(test == 0) out[ind, ] else chr_bins[n, ]
  }
  
  ## assign the data frame w/ each peak and peak bin position to a list
  ## that we'll collapse after we've run through all the reporters 
  final[[k]] <- chr_bins
}

## collapse output to a single frame
test <- do.call("rbind", final)


#order rows
#by reporter
test$reporter.hm <- factor(test$reporter.hm, levels=c("Rpn4 in BTZ", "Rpn4 in AZC", "Rpn4 in 4NQO", "Rpn4 in LiAc", "Rpn4 in YNB", "Rpn4 in Low N", "Rpn4 in Low G", "Rpn4 in SC", 
                                                      "4xUb in BTZ", "4xUb in AZC", "4xUb in 4NQO", "4xUb in LiAc", "4xUb in YNB", "4xUb in Low G", "4xUb in SC", 
                                                      "UFD in BTZ", "UFD in AZC", "UFD in 4NQO", "UFD in LiAc", "UFD in YNB", "UFD in Low N", "UFD in Low G", "UFD in SC", 
                                                      "Thr in BTZ", "Thr in AZC", "Thr in 4NQO", "Thr in LiAc", "Thr in YNB", "Thr in Low N", "Thr in Low G", "Thr in SC", 
                                                      "Phe in BTZ", "Phe in AZC", "Phe in 4NQO", "Phe in LiAc", "Phe in YNB", "Phe in Low N", "Phe in Low G", "Phe in SC", 
                                                      "Asn in BTZ", "Asn in AZC", "Asn in 4NQO", "Asn in LiAc", "Asn in YNB", "Asn in Low N", "Asn in Low G", 
                                                      "Asn in SC"), order=TRUE)

#test$reporter.hm <- factor(test$reporter.hm, levels = unique(test$reporter.hm))


#by environment
# test$reporter.hm <- factor(test$reporter.hm, levels=c("BTZ | Rpn4", "BTZ | 4xUb", "BTZ | UFD", "BTZ | Thr", "BTZ | Phe", "BTZ | Asn", 
#                                                       "AZC | Rpn4", "AZC | 4xUb", "AZC | UFD", "AZC | Thr", "AZC | Phe", "AZC | Asn", 
#                                                       "4NQO | Rpn4", "4NQO | 4xUb", "4NQO | UFD", "4NQO | Thr", "4NQO | Phe", "4NQO | Asn", 
#                                                       "LiAc | Rpn4", "LiAc | 4xUb", "LiAc | UFD", "LiAc | Thr", "LiAc | Phe", "LiAc | Asn", 
#                                                       "YNB | Rpn4", "YNB | 4xUb", "YNB | UFD", "YNB | Thr", "YNB | Phe", "YNB | Asn",
#                                                       "Low N | Rpn4", "Low N | UFD", "Low N | Thr", "Low N | Phe", "Low N | Asn", 
#                                                       "Low G | Rpn4", "Low G | 4xUb", "Low G | UFD", "Low G | Thr", "Low G | Phe", "Low G | Asn", 
#                                                       "SC | Rpn4", "SC | 4xUb", "SC | UFD", "SC | Thr", "SC | Phe", 
#                                                       "SC | Asn"), order=TRUE)



## position where the labels for the x axis (chr) go
chr_labels <- sapply(X = 1:16, FUN = function(x) {
  mean(test$bin[test$chr == x])
})

## position where the dividing lines on the x axis go
chr_cutoffs <- sapply(X = 1:15, FUN = function(x) {
  max(test$bin[test$chr == x]) + 5e4
})


## this object demarcates the bins on the colorkey
## 0.075 is the smallest delta_AF value to display for the QTLs
#by looking at the df, the min and max is just within -0.7 and 0.7
delta_AF_ats = c(-0.7, -0.50, -0.25, -0.075, 0.075, 0.25, 0.50, 0.7)

## "RdBu" is colorlbind friendly:
delta_AF_cols <- c(brewer.pal(10, "RdBu")[10:7],
                   "white", 
                   brewer.pal(10, "RdBu")[4:1])


#using df with bins to see how many unique regions/bins the 416 qtls map to
savetest = test
binnedqtls = test
binnedqtls$empty = paste0(binnedqtls$left_Index, '_', binnedqtls$max_Index, '_', binnedqtls$right_Index)
binnedqtls = binnedqtls[binnedqtls$empty != '0_0_0',]
#binnedqtls <- binnedqtls[!duplicated(binnedqtls), ]
binnedqtls$uniquebins = paste0(binnedqtls$chr, '_', binnedqtls$bin)
length(unique(binnedqtls$uniquebins))
#78 unique bins for the 416 QTLs!

#how many bins overall?
test$uniquebins = paste0(test$chr, '_', test$bin)
length(unique(test$uniquebins))
#128

#separate reporter and environment

binnedqtls$rep = str_split_i(binnedqtls$reporter.hm, ' in ', 1)
binnedqtls$env = str_split_i(binnedqtls$reporter.hm, ' in ', 2)

write.table(binnedqtls,
            file = "416QTLsinBins10.17.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

write.table(binnedqtls,
            file = "416QTLsinKnownLociBins10.31.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')
#further analyzing these in ~/myproteasome/paper_GxEQTLbins24.10.20.R


binnedqtls = read.table(file = "416QTLsinBins10.17.24.txt", header = T, sep = '\t')

#how many QTLs go into each bin?
totqtlsinbins = binnedqtls %>% count(uniquebins)
colnames(totqtlsinbins) <- c('uniquebins', 'allQTLs')

ggplot(data = totqtlsinbins, aes(y=n, x=reorder(uniquebins, -n))) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text=element_text(size=8),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



{
  pdf(file = paste0(proj_dir, "2024.10.31_delta_AF_heatmap.pdf"),
      height = 10, width = 15)
  print( 
    levelplot(delta_AF ~ bin * reporter.hm,
              data = test,
              xlab = "Chromosome",
              ylab = "Reporter and Environment", #rra this didn't show up
              col.regions = delta_AF_cols,
              at = delta_AF_ats,
              colorkey = list(col = delta_AF_cols,
                              at = 0:7,
                              labels = as.character(c(-0.7, -0.50, -0.25, -0.075, 
                                                      0.075, 0.25, 0.50, 0.7)),
                              ## reduce the height of the colorkey a bit:
                              height = 0.5,
                              title = "delta AF"),
              scales = list(x = list(at = chr_labels,
                                     tck = c(1, 0),
                                     alternating = 1,
                                     labels = as.roman(1:16)),
                            y = list(alternating = 1,
                                     tck = c(1, 0))),
              par.settings = list(layout.widths = list(left.padding = 3,
                                                       right.padding = 3.5),
                                  par.ylab.text = list(cex = .75, #rra changed from 1.24
                                                       col = "white")),
              panel = function(...){
                panel.levelplot(...)
                panel.abline(v = chr_cutoffs,
                             lty = 1, col = gray(0.4))
                panel.abline(h = c(8.5, 15.5,23.5,31.5,39.5), #horizontal lines for ordered by reporter. Start with 0.5, and each row is 1.0
                #panel.abline(h = c(6.5, 12.5,18.5,24.5,30.5,35.5,41.5), #horizontal lines for ordered by environment
                             lty = 1, col = gray(0.4))
              })
  )

  dev.off()
}