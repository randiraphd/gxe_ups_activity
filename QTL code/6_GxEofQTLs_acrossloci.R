#Code by Randi Avery for 
  #Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity. bioRxiv, 2024.11.21.624644. https://doi.org/10.1101/2024.11.21.624644

#binning genome to find where GxE QTLs fit


base_dir <- "~/myproteasome/"
proj  <- "2024.07.05_allHL_GxE/"
proj_dir <- paste0(base_dir, proj)
peak_dir <- paste0(base_dir, proj, "peaks/")
setwd(proj_dir)

mymaroon = '#862334'
mygold = '#FFD966'
col_by <- "#2166ACFF" 
col_rm <- "#BF3232FF"

library(ggplot2)
library(dplyr)
library(MASS)
library(tidyr)
library(stringr)
library(janitor)

#load in overlaps
overlaps_merge_filtered = read.table(file = "overlappingQTLs_forGxE_24.07.24_filtered.txt", header = TRUE, sep = '\t')

#separate out P/A vs sign change

signchangeGxE = overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir == 'FALSE', ]
paGxE = overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir == 'Absent', ]

#average sign change values
signchangeGxE$reporter.hm = paste0(signchangeGxE$rep1, '_and_', signchangeGxE$environment_2)
signchangeGxE$chr = signchangeGxE$chrm
#Repurposing previous code: lod and dAF are not used here. We just need the same df format for the code, and the position of the QTL
signchangeGxE$LOD = (as.numeric(signchangeGxE$rep_1_LOD) + as.numeric(signchangeGxE$rep_2_LOD))/2
signchangeGxE$delta_AF = (as.numeric(signchangeGxE$rep_1_delta_AF) + as.numeric(signchangeGxE$rep_2_delta_AF))/2
signchangeGxE$left_Index = (as.numeric(signchangeGxE$rep_1_left_Index) + as.numeric(signchangeGxE$rep_2_left_Index))/2
signchangeGxE$max_Index = (as.numeric(signchangeGxE$rep_1_max_Index) + as.numeric(signchangeGxE$rep_2_max_Index))/2
signchangeGxE$right_Index = (as.numeric(signchangeGxE$rep_1_right_Index) + as.numeric(signchangeGxE$rep_2_right_Index))/2

keeps = c("reporter.hm", "chr",         "LOD" ,        "delta_AF" ,   "left_Index",  "max_Index",   "right_Index")
avgsignchGxE = signchangeGxE[keeps]
#use this for plotting

paGxE$reporter.hm = paste0(paGxE$repA2, '_and_', paGxE$repB2)
paGxE$chr = paGxE$chrm
paGxE$LOD = as.numeric(paGxE$rep_1_LOD)
paGxE$delta_AF = as.numeric(paGxE$rep_1_delta_AF)
paGxE$left_Index = as.numeric(paGxE$rep_1_left_Index)
paGxE$max_Index = as.numeric(paGxE$rep_1_max_Index)
paGxE$right_Index = as.numeric(paGxE$rep_1_right_Index)

paGxE.hm = paGxE[keeps]
#use this for plotting


#Lines 63-200 are adapted from code from Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
#and Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570

#paste from heatmap code
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


reporters <- unique(avgsignchGxE$reporter.hm)

##df has to have these columns in this order:
# chr_bins <- data.frame(reporter.hm = rep(reporters[k], length(bins)),
#                        chr = chr_indices,
#                        LOD = rep(0, length(bins)),
#                        delta_AF = rep(0, length(bins)),
#                        left_Index = rep(0, length(bins)),
#                        max_Index = rep(0, length(bins)),
#                        right_Index = rep(0, length(bins)),
#                        bin = bins)

## <<build_heatmap_dataframe>>
#k <- 1

final <- list()
#change for each hm running on:
hm = paGxE.hm
reporters <- unique(hm$reporter.hm)
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
  out <- hm[hm$reporter.hm == reporters[k], ]
  
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
gxetest <- do.call("rbind", final)

#using df with bins to see how many unique regions/bins the 416 qtls map to

binnedsignchangeGxE = gxetest
#run again for PA
binnedPAGxE = gxetest


binnedsignchangeGxE$empty = paste0(binnedsignchangeGxE$left_Index, '_', binnedsignchangeGxE$max_Index, '_', binnedsignchangeGxE$right_Index)
binnedsignchangeGxE = binnedsignchangeGxE[binnedsignchangeGxE$empty != '0_0_0',]
binnedsignchangeGxE$uniquebins = paste0(binnedsignchangeGxE$chr, '_', binnedsignchangeGxE$bin)
length(unique(binnedsignchangeGxE$uniquebins))
#only 8 bins! But the QTLs on 14 are split into 2 bins, even tho the peaks are only 56,625 bp apart

binnedPAGxE$empty = paste0(binnedPAGxE$left_Index, '_', binnedPAGxE$max_Index, '_', binnedPAGxE$right_Index)
binnedPAGxE = binnedPAGxE[binnedPAGxE$empty != '0_0_0',]
binnedPAGxE$uniquebins = paste0(binnedPAGxE$chr, '_', binnedPAGxE$bin)
length(unique(binnedPAGxE$uniquebins))
#70 unique bins for all 254 P/A GxE!


#how many QTLs go into each bin?
totbinnedPAGxE = binnedPAGxE %>% count(uniquebins)
totbinnedsignchangeGxE  = binnedsignchangeGxE %>% count(uniquebins)
#want to combine as well
totbinnedPAGxE$GxE = 'PA'
totbinnedsignchangeGxE$GxE = 'SignChange'
totbinnedGxE = rbind(totbinnedPAGxE, totbinnedsignchangeGxE)

ggplot(data = totbinnedGxE, aes(fill = GxE, y=n, x=reorder(uniquebins, n, sum))) +
  geom_bar(position="stack", stat="identity") +
  theme(axis.text=element_text(size=8),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_x_discrete(limits = rev) #reverses order of x-axis


#see if category (ubiquitin independent vs dependent | starvation/not) pile up at certain loci. Do together for sign change and P/A GxE
binnedPAGxE$GxE = 'PA'
binnedPAGxE$reporter = str_split_i((binnedPAGxE$reporter.hm), '_', 1)
binnedPAGxE$environment = str_split_i((binnedPAGxE$reporter.hm), '_', -1)
binnedsignchangeGxE$GxE = 'signchange'
binnedsignchangeGxE$reporter = str_split_i((binnedsignchangeGxE$reporter.hm), '_', 1)
binnedsignchangeGxE$environment = str_split_i((binnedsignchangeGxE$reporter.hm), '_', -1)
binnedGxE = rbind(binnedPAGxE, binnedsignchangeGxE)

#another way to get there, from above, but doesn't keep all the info
totbinnedGxE = binnedGxE %>% count(uniquebins)

binnedGxE$ubsystem = ifelse(binnedGxE$reporter == '4xUb' |binnedGxE$reporter == 'Rpn4', 'Independent', 'Dependent')
binnedGxE$starve = ifelse(binnedGxE$environment == 'Low N' |binnedGxE$environment == 'Low G'|binnedGxE$environment == 'YNB' , 'Starvation', 'Not Starvation')
#each row would be one count, so adding a column with 1, so they pile up
binnedGxE$n = 1

write.table(binnedGxE,
            file = "271GxEpairsinBins24.10.22.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

binnedGxE = read.table(file = "271GxEpairsinBins24.10.22.txt", header = TRUE, sep = '\t')
known_binnedGxE = read.table(file = "knownloci_binnedGxEcomparisons24.10.31.txt", header = TRUE, sep = '\t')

#naming BUL2
known_binnedGxE[known_binnedGxE == "13_9e+06"] = "BUL2"
known_binnedGxE[known_binnedGxE == "PA"] = "Presence / Absence"
known_binnedGxE[known_binnedGxE == "signchange"] = "Sign Change"

#manually putting in labels
#want labels of 'known' loci over bars
#actually only need one label per bar, and then need position
known_binnedGxE$name = NA
known_binnedGxE['22', "name"] = "UBR1"
known_binnedGxE['41', "name"] = "UBC6"
known_binnedGxE['5', "name"] = "RPT6"
known_binnedGxE['233', "name"] = "NTA1"
known_binnedGxE['38', "name"] = "MKT1"
known_binnedGxE['28', "name"] = "IRA2"
known_binnedGxE['131', "name"] = "HAP1"
known_binnedGxE['99', "name"] = "DOA10"
known_binnedGxE['139', "name"] = "BUL2"

#got position from other DF with totals
known_binnedGxE$pos = NA
known_binnedGxE['22', "pos"] = 7
known_binnedGxE['41', "pos"] = 1
known_binnedGxE['5', "pos"] = 10
known_binnedGxE['233', "pos"] = 1
known_binnedGxE['38', "pos"] = 8
known_binnedGxE['28', "pos"] = 15
known_binnedGxE['131', "pos"] = 11
known_binnedGxE['99', "pos"] = 5
known_binnedGxE['139', "pos"] = 14


#Figure 5A
ggplot(data = known_binnedGxE, aes(fill = GxE, y=n, x=reorder(uniquebins, n, sum))) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c(mymaroon, mygold))+
  theme(text=element_text(size=30),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank()) +
  labs( x = "Genomic bins (n = 71)", y = 'Number of GxE loci') +
  scale_x_discrete(limits = rev)+ #reverses order of x-axis
  #geom_text(aes(x=reorder(uniquebins, n, sum), y=pos, label=name, fontface = 'italic'),  
  #          angle = 90, hjust = -0.1) + #, hjust = 0.1, size = 4.5, vjust = 0.1)+
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(0,16,5)) #reposition bottom of x axis
#saved as 8.5x17 inches pdf



##Testing INTEGRATING PARENTAL GXE WITH QTL GXE##
#scatter plot - 40 points
#one axis, p-value for parental GxE
#other axis, number of GxE pairs of QTLs

#amount of GxE in QTLs (DONT use proportion - using amount basically takes into account number of QTLs as well)
stats_done = read.table(file = "GxEstatsSample24.07.31.txt", header = TRUE, sep = '\t', check.names=FALSE)
stats_done$ubsystem = ifelse(stats_done$rep == '4xUb' |stats_done$rep == 'Rpn4', 'Independent', 'Dependent')
stats_done$starve = ifelse(stats_done$environment == 'Low N' |stats_done$environment == 'Low G'|stats_done$environment == 'YNB' , 'Starvation', 'Not Starvation')
stats_done_conclude = stats_done[stats_done$environment != 'SC', ]
#removing sample that I don't have data for parents
stats_done_conclude = stats_done_conclude[stats_done_conclude$Sample != 'UFD_in_Low N', ]

#GxE strength in parents
results_table_loess = read.table(file = "~/myproteasome/reporter_characterization/environments/combined/results_lmer_40tests_10.17.24_TFTtime_loess.txt", header = TRUE, sep = '\t', check.names=FALSE)
results_table_loess[results_table_loess == "LowG"] = "Low G"
results_table_loess[results_table_loess == "LowN"] = "Low N"
results_table_loess$Sample = paste0(results_table_loess$reporter, "_in_", results_table_loess$environment2)

#conclusions concatenated 
concludescat = merge(stats_done_conclude, results_table_loess, by = 'Sample')
concludescat$proportion = concludescat$`Total GxE` / concludescat$`Total Overlapping Pairs`
concludescat$`Total Overlapping Pairs`

#pvalue, amount of GxE, color/shape ubsystem/starve
ggplot(concludescat, aes(x = `Total GxE`, y = -log10(interaction_Pval))) +
  geom_point(aes(color = ubsystem, shape = starve))+
  xlab('Total pairs of QTLs with GxE') +
  ylab('p-value of GxE interaction in parents (-log10)')
#not really interesting

#correlate number of QTLs and amount of GxE at QTLs in that sample
#taking out SC
stats_done_noSC = stats_done[stats_done$environment != 'SC', ]
totqtls$Sample = paste0(totqtls$reporter, "_in_", totqtls$environment)
gxe = merge(stats_done_noSC, totqtls, by = 'Sample')

ggplot(gxe, aes(y = `Total GxE`, x = n)) +
  geom_point(aes(color = ubsystem.x, shape = starve))+
  ylab('Total pairs of QTLs with GxE') +
  xlab('Total QTLs')
  
cor.test(gxe$n, gxe$`Total GxE`, method = 'pearson')


#plot showing number of QTLs in each bin, and then GxE pairs in those same bins, to see if they're correlated

#all QTLs
binnedqtls = read.table(file = "416QTLsinBins10.17.24.txt", header = T, sep = '\t')
#how many QTLs go into each bin?
totqtlsinbins = binnedqtls %>% count(uniquebins)
colnames(totqtlsinbins) <- c('uniquebins', 'allQTLs')

binnedGxE = read.table(file = "271GxEpairsinBins24.10.22.txt", header = TRUE, sep = '\t')
totGxEqtlsinbins = binnedGxE %>% count(uniquebins)
colnames(totGxEqtlsinbins) <- c('uniquebins', 'GxEQTLpairs')

#does more QTLs mean more GxE? use bins
#i'm out of names
df2 = merge(totGxEqtlsinbins, totqtlsinbins, by = 'uniquebins', all = TRUE) 
#there are some bins that have QTLs, but no GxE, so put 0's there
df2[is.na(df2)] <- 0

##point out top bins
#don't have the gene names in there - just add manually
df2$name = NA
df2["69", "name"] = "RPT6"
df2["29", "name"] = "MKT1"
df2["18", "name"] = "HAP1"
df2["32", "name"] = "IRA2"


#Fig 5B
ggplot(df2, aes(y = GxEQTLpairs, x = allQTLs, stroke = NA)) +
  geom_point(alpha = 0.4, position=position_jitter(height=.2, width=.2), size = 2.75)+
  ylab('Total pairs of QTLs with GxE') +
  xlab('Total QTLs') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_line(color="grey", linewidth = 0.3),
        axis.line = element_line(colour = "grey"))+ 
  ggtitle('Numbers of QTLs and GxE pairs in each bin that isnt empty') +
  geom_text(aes(allQTLs, GxEQTLpairs, label=name, fontface = 'italic'), position = position_dodge(width = 1),
            vjust = -0.6) +
  scale_y_continuous(breaks = seq(0,16,5)) 


cor.test(df2$allQTLs, df2$GxEQTLpairs, method = 'spearman')


#thought it looked like there was a limit to GxE, so i cut off the top 4 points, but the correlation was actually weaker
#test = df[df$allQTLs < 26,]
#cor.test(test$allQTLs, test$GxEQTLpairs, method = 'pearson')


#Want more info on the known loci - made following table in heatmap script

knownloci_binnedQTLs = read.table(file = "416QTLsinKnownLociBins10.31.24.txt", header = TRUE, sep = '\t')
#now just name bins based off the locus - just opened the file in Excel and looked at the unique bin names
knownloci_binnedQTLs[knownloci_binnedQTLs == "7_4600000"] = "RPT6"
knownloci_binnedQTLs[knownloci_binnedQTLs == "5_3600000"] = "UBC6"
knownloci_binnedQTLs[knownloci_binnedQTLs == "7_5e+06"] = "UBR1"
knownloci_binnedQTLs[knownloci_binnedQTLs == "9_6100000"] = "DOA10_left" #straddles two bins
knownloci_binnedQTLs[knownloci_binnedQTLs == "9_6200000"] = "DOA10_right"
knownloci_binnedQTLs[knownloci_binnedQTLs == "10_6900000"] = "NTA1"
knownloci_binnedQTLs[knownloci_binnedQTLs == "12_8500000"] = "HAP1"
knownloci_binnedQTLs[knownloci_binnedQTLs == "14_10400000"] = "MKT1"
knownloci_binnedQTLs[knownloci_binnedQTLs == "14_10500000"] = "MKT1_maybe"
knownloci_binnedQTLs[knownloci_binnedQTLs == "15_10900000"] = "IRA2"

#since nothing went to DOA10_right, can simplify that name
knownloci_binnedQTLs[knownloci_binnedQTLs == "DOA10_left"] = "DOA10"
knownloci_binnedQTLs[knownloci_binnedQTLs == "DOA10_right"] = "DOA10"

#Remove rows with NAs (those were the dummy QTLs I made for binning)
knownloci_binnedQTLs = na.omit(knownloci_binnedQTLs)
totqtlsinbins_named = knownloci_binnedQTLs %>% count(uniquebins)
#change n to 'all QTLs'
colnames(totqtlsinbins_named) <- c('uniquebins', 'allQTLs')

write.table(knownloci_binnedQTLs,
            file = "knownloci_binnedQTLs24.10.31.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')



#now add GxE, and specifically sign changes 
binnedGxE[binnedGxE == "7_4600000"] = "RPT6"
binnedGxE[binnedGxE == "5_3600000"] = "UBC6"
binnedGxE[binnedGxE == "7_5e+06"] = "UBR1"
binnedGxE[binnedGxE == "9_6100000"] = "DOA10" #straddles two bins
binnedGxE[binnedGxE == "9_6200000"] = "DOA10"
binnedGxE[binnedGxE == "10_6900000"] = "NTA1"
binnedGxE[binnedGxE == "12_8500000"] = "HAP1"
binnedGxE[binnedGxE == "14_10400000"] = "MKT1"
binnedGxE[binnedGxE == "15_10900000"] = "IRA2"
binnedGxE[binnedGxE == "14_10500000"] = "MKT1_maybe"

write.table(binnedGxE,
            file = "knownloci_binnedGxEcomparisons24.10.31.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

totGxEqtlsinbins_named = binnedGxE %>% count(uniquebins)
colnames(totGxEqtlsinbins_named) <- c('uniquebins', 'GxEQTLcomparisons')

#now only want to know where sign changes are
signchangebinned = binnedGxE[binnedGxE$GxE == "signchange", ]
totsignchangebins_named = signchangebinned %>% count(uniquebins)
colnames(totsignchangebins_named) <- c('uniquebins', 'signchangeOnly')

#merge by column - but put zeros in for missing
inknownloci = merge(totqtlsinbins_named, totGxEqtlsinbins_named, all = TRUE)
inknownloci[is.na(inknownloci)] <-0
inknownloci = merge(inknownloci, totsignchangebins_named, all = TRUE)
inknownloci[is.na(inknownloci)] <-0


write.table(inknownloci,
            file = "knownloci_totals24.10.31.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#do a for loop to see if dAF is all one direction for the bins
#just gonna write by hand
knownloci = c("UBR1" , "UBC6" , "RPT6" , "NTA1" , "MKT1_maybe" , "MKT1" , "IRA2" , "HAP1" , "DOA10")
#for (l in knownloci){
  locus = knownloci_binnedQTLs[knownloci_binnedQTLs$uniquebins == "UBR1", ]
#  print(l)
#  print(locus$delta_AF)
#}
  
#for all 416 QTLs
  # UBR1 - diff signs
  # UBC6 - all positive!
  # RPT6 - all positive!
  # NTA1 - one neg out of nine, and it's the furthest QTL to the right in the locus
  # MKT1_maybe - the Rpn4 in AZC is positive
  # MKT1 - all diff
  # IRA2 - all diff, majority negative tho
  # HAP1 - all diff, majority negative tho
  # DOA10 - 2 neg out of 10, and those are at the very left of the bin, while DOA10 is on the right of the bin (checked and nothing to the next bin over, so captured all DOA10 loci)


