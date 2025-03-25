#Edited for use by Randi Avery in 2024 from code by Mahlon Collins:
  #Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
  #Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570


#produces QTL traces and LOD plots - provided in Supplementary File 3, Supplementary Figure 3, and Figure 4B

base_dir <- "~/myproteasome/" 


#############
## USER INPUT
#############
## set the specific directory you'll work 
## in and name the comparison table
## TRAILING SLASH AT END OF DIR
## below, your project, e.g.,
## "2020.08.17_FPFA002_TDH3pr_Arg_N-end_TFT_sorts/"
proj           <- "2024.07.05_allHL_GxE/"
proj_dir       <- paste0(base_dir, proj)
#copy and paste the file names from the previous output into the replicates comparison table
comp_table     <- paste0(proj_dir, "2024.07.05_all_replicates_comparison_table_allGxE.txt")

mpr            <- "_multipoolResults"
rd             <- ".RData"
pop            <- "_pop_0"
#################
## END USER INPUT
#################

needed.dirs <- c("results/", "rdata/", "/peaks")

dir.maker <- function(x){if(!dir.exists(paths = paste0(proj_dir, x)))
                             dir.create(path = paste0(proj_dir, x))}

sapply(X = needed.dirs, FUN = dir.maker) 

results_dir     <- paste0(proj_dir, "results/")
data_dir      <- paste0(proj_dir, "r_output/") 
peaks_dir       <- paste0(proj_dir, "peaks/")
experiment_file <- as.data.frame(read.table(comp_table,
                                            stringsAsFactors=FALSE,
                                            header = TRUE))

## SNPs is a giant table w/ SNP positions 
SNPs <- read.table("~/myproteasome/2024.07.05_allHL_GxE/QTL_scripts/SNPs_Maggie_170809_BY_positions.txt",
                   stringsAsFactors = FALSE,
                   head = FALSE)

for (thisChr in unique(SNPs[,1])){
    SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])
}

SNPs <- rbind(SNPs[SNPs[,1] == "chrI",],
              SNPs[SNPs[,1] == "chrII",],
              SNPs[SNPs[,1] == "chrIII",],
              SNPs[SNPs[,1] == "chrIV",],
              SNPs[SNPs[,1] == "chrV",],
              SNPs[SNPs[,1] == "chrVI",],
              SNPs[SNPs[,1] == "chrVII",],
              SNPs[SNPs[,1] == "chrVIII",],
              SNPs[SNPs[,1] == "chrIX",],
              SNPs[SNPs[,1] == "chrX",],
              SNPs[SNPs[,1] == "chrXI",],
              SNPs[SNPs[,1] == "chrXII",],
              SNPs[SNPs[,1] == "chrXIII",],
              SNPs[SNPs[,1] == "chrXIV",],
              SNPs[SNPs[,1] == "chrXV",],
              SNPs[SNPs[,1] == "chrXVI",])


withMultipool <- TRUE

## common annotations, functions, etc ----------------
## check for Bioconductor and install if not available
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")

bioc_package_installer <- function(x) {
    if (!requireNamespace(x))
        BiocManager::install(x, INSTALL_opts = "--no-lock")
}

bioc_package_installer("VariantAnnotation")

library("VariantAnnotation")
source(paste0(proj_dir, "QTL_scripts/gTest.R"))
source(paste0(proj_dir, "QTL_scripts/x_qtl_seq_functions_170831.R"))
source(paste0(proj_dir, "QTL_scripts/mp_JB_170901.R"))
source(paste0(proj_dir, "QTL_scripts/peaksFromVector.R"))

## data frame with all yeast genes, plus
## chr, pos., strand, and names
geneInfo <- read.table(paste0(proj_dir, "QTL_scripts/ensemblGenes_ensembl83_160307_MOD.txt"),
                       stringsAsFactors = FALSE,
                       sep = "\t",
                       header = TRUE)

## rownames become systematic names
rownames(geneInfo) <- geneInfo[, "geneID"]

## "geneName" is the common name, e.g., 'HOG1'
## for some (many?) rows of 'geneInfo', there
## is no 'geneName', so it's just an empty string
## e.g., head(allNames)
allNames <- geneInfo[, "geneName"]

names(allNames) <- geneInfo[, 1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

sepBetweenChr <- 1e5
trimFromEnd   <- 15e3
obsMin        <- 10
LoessSpan     <- 0.1
## same as in Albert 2014
AFThres       <- 0.09653124
## multipool LOD threshold, determined by null sorts
multiThres    <- 4.5
mpr  <- "_multipoolResults"
rd   <- ".RData"
#p    <- 35

for (p in 1:nrow(experiment_file)) {

    ## load 'multipoolResults.RData' for replicate 1
    ## adds 'multiPeaks' and 'multipoolOutput'
    load(paste0(data_dir,
                experiment_file[p, 1],
                mpr, rd))
    print(experiment_file[p, 1])
    ## load ".RData" for replicate 1
    ## adds 'theseCounts'
    load(paste0(data_dir,
                experiment_file[p, 1], 
                rd))

    multipeaks_one    <- multiPeaks
    multipool_out_one <- multipoolOutput
    these_counts_one  <- theseCounts

    #get peak information again to include allele frequency difference this time for the peaks
    multi_table_one <- do.call(rbind.data.frame, multipeaks_one)
    
    ## get the chromosome number by looping over the
    ## multiPeaks list and repeating the chromosome number
    ## for each QTL in a chromosome.  if there are no QTLs
    ## on a chromosome, do nothing
    multi_chr_number_one <- vector()
    for(l in 1:length(multipeaks_one)){
      if(!is.null(multipeaks_one[[l]]))
        multi_chr_number_one <- c(multi_chr_number_one,
                              rep(l, nrow(multipeaks_one[[l]])))
    }
    
    multi_table_one$chr <- multi_chr_number_one
    
    ## these need to be 'if' clauses to account for
    ## situations where there are no multiPool peaks
    ## re-order the output 
    if(nrow(multi_table_one) >= 1)
      multi_table_one <- multi_table_one[, c(5, 2, 1, 3, 4)]
    
    ## re-name a couple of columns for easier reading
    if(nrow(multi_table_one) >= 1)
      colnames(multi_table_one) <- c("chr", "LOD", "max_Index",
                                 "left_Index", "right_Index")
    
    
    ## load 'multipoolResults.RData' for replicate 2
    ## adds 'multiPeaks' and 'multipoolOutput'
    load(paste0(data_dir,
                experiment_file[p, 2],
                mpr, rd))

    ## load ".RData" for replicate 2
    ## adds 'theseCounts'
    load(paste0(data_dir,
                experiment_file[p, 2],
                rd))

    multipeaks_two    <- multiPeaks
    multipool_out_two <- multipoolOutput
    these_counts_two  <- theseCounts
    
    #get peak information again to include allele frequency difference this time for the peaks
    multi_table_two <- do.call(rbind.data.frame, multipeaks_two)
    
    ## get the chromosome number by looping over the
    ## multiPeaks list and repeating the chromosome number
    ## for each QTL in a chromosome.  if there are no QTLs
    ## on a chromosome, do nothing
    multi_chr_number_two <- vector()
    for(l in 1:length(multipeaks_two)){
      if(!is.null(multipeaks_two[[l]]))
        multi_chr_number_two <- c(multi_chr_number_two,
                                  rep(l, nrow(multipeaks_two[[l]])))
    }
    
    multi_table_two$chr <- multi_chr_number_two
    
    ## these need to be 'if' clauses to account for
    ## situations where there are no multiPool peaks, #rra, but need to also allow when there is only 1 peak
    ## re-order the output 
    if(nrow(multi_table_two) >= 1)
      multi_table_two <- multi_table_two[, c(5, 2, 1, 3, 4)]
    
    ## re-name a couple of columns for easier reading
    if(nrow(multi_table_two) >= 1)
      colnames(multi_table_two) <- c("chr", "LOD", "max_Index",
                                     "left_Index", "right_Index")

        ## compute BY allele frequencies
        ## and sequencing coverage 
        allele_f_calc <- function(sample, pool, ...) {
            n  <- sample[, paste0(pool, "_ref")]
            d1 <- sample[, paste0(pool, "_ref")]
            d2 <- sample[, paste0(pool, "_alt")]
            n / (d1 + d2)
        }

        cover_calc <- function(sample, pool, ...) {
            a <- sample[, paste0(pool, "_ref")]
            r <- sample[, paste0(pool, "_alt")]
            a + r
        }

        these_counts_one$h_BY_AF <- allele_f_calc(these_counts_one, "high")
        these_counts_one$l_BY_AF <- allele_f_calc(these_counts_one, "low")

        these_counts_two$h_BY_AF <- allele_f_calc(these_counts_two, "high")
        these_counts_two$l_BY_AF <- allele_f_calc(these_counts_two, "low")

        these_counts_one$h_cover <- cover_calc(these_counts_one, "high")
        these_counts_one$l_cover <- cover_calc(these_counts_one, "low")

        these_counts_two$h_cover <- cover_calc(these_counts_two, "high")
        these_counts_two$l_cover <- cover_calc(these_counts_two, "low")

        gcoords <- getGcoords(these_counts_two$chr,
                              these_counts_two$pos,
                              sepBetweenChr)

        ## get plot coordinates for the chromosome line dividers
        chrCutoffs <- sapply(unique(theseCounts$chr),
                             function(x){gcoords[theseCounts$chr == x][1] - sepBetweenChr/2})

        names(chrCutoffs) <- unique(theseCounts$chr)

        ## position of the text on the x axis for the chromosomes 
        chrLabels <- sapply(1:(length(chrCutoffs)-1),
                           function(i)(chrCutoffs[i] + chrCutoffs[i+1])/2)

        ## add half the length of chrXVI
        chrLabels <- c(chrLabels, chrCutoffs[16] + sepBetweenChr + 948066/2)

        names(chrLabels)[16] = "chrXVI"

        ## make the ylimit of the plot the max of either:
        ## 1. the significance threshold
        ## 2. the highest peak of the first replicate
        ## 3. the highest peak of the second replicate
        ylimMax <- max(c(multiThres,
                         sapply(c(multipool_out_one,
                                  multipool_out_two),
                                function(x) {
                                    max(x[[2]][, 2])}))) + 1

        rep_one_col  <- "#3A3A3A"     ## dark gray (back)
        rep_two_col  <- "#8A8A8A"     ## light gray (front)
        sig_line_col <- "#981F1FCC"   ## red

        roll_high_one <- rollLoessByChrWithWeights(data.frame(these_counts_one[, "chr"],
                                                              these_counts_one[, "h_BY_AF"],
                                                              gcoords,
                                                              median(these_counts_one[, "h_cover"]),
                                                              stringsAsFactors = FALSE),
                                                   LoessSpan)

        roll_low_one <- rollLoessByChrWithWeights(data.frame(these_counts_one[, "chr"],
                                                             these_counts_one[, "l_BY_AF"],
                                                             gcoords,
                                                             median(these_counts_one[, "l_cover"]),
                                                             stringsAsFactors = FALSE),
                                                  LoessSpan)

        roll_high_two <- rollLoessByChrWithWeights(data.frame(these_counts_two[, "chr"],
                                                              these_counts_two[, "h_BY_AF"],
                                                              gcoords,
                                                              median(these_counts_two[, "h_cover"]),
                                                              stringsAsFactors = FALSE),
                                                   LoessSpan)

        roll_low_two <- rollLoessByChrWithWeights(data.frame(these_counts_two[, "chr"],
                                                             these_counts_two[, "l_BY_AF"],
                                                             gcoords,
                                                             median(these_counts_two[, "l_cover"]),
                                                             stringsAsFactors = FALSE),
                                                  LoessSpan)

        ## code below returns the height on the y axis
        ## for the '*''s on the allele frequency difference
        ## plot.  it's called for ea. of the 2 replicates
        ## first, initialize 2 empty lists to put 
        ## (roll_high - roll_low) values into
        ## y1_afd = replicate 1
        ## y2_afd = replicate 2
        y1_afd <- lapply(1:16, function(x) {
                             NULL
                         })
        y2_afd <- y1_afd
        y1afd.df = y1_afd #dfs to save allele frequencies before altering them for star placement
        y2afd.df = y1_afd
        
        ## chromosome names for pattern matching 
        chr_indices <- sapply(1:16, function(x) {
                                  paste0("chr",
                                         as.character(as.roman(x)))
                              })

        ## name the elements according to chromosomes 
        names(y1_afd) <- chr_indices
        names(y2_afd) <- chr_indices
        names(y1afd.df) <- chr_indices
        names(y2afd.df) <- chr_indices
        
        for (y in seq_along(multipeaks_one)) {
            if (!is.null(multipeaks_one[[y]]))
                for (q in seq_along(multipeaks_one[[y]][, 1])) {
                    ## has to return an index you can feed to roll_high_one
                    h_roll_chr <- roll_high_one[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                     x = names(roll_high_one))]
                    l_roll_chr <- roll_low_one[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                    x = names(roll_low_one))]
                    t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
                    t_min <- which.min(abs(multipeaks_one[[y]][q, "maxIndex"] - t_snp))
                    y1_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
                    y1afd.df[[y]][q] = y1_afd[[y]][q] #save for later to make heatmap with allele freq diffs
                    #if there should be a star, add or subtract 0.05 so it goes right above/below the peak
                    y1_afd[[y]][q] <- ifelse(y1_afd[[y]][q] > 0,
                                      ifelse(y1_afd[[y]][q] > AFThres,
                                             y1_afd[[y]][q] + 0.05,
                                             AFThres + 0.05),
                                      ifelse(y1_afd[[y]][q] < -AFThres,
                                             y1_afd[[y]][q] - 0.05,
                                             -AFThres - 0.05))
                }
        }
        
        #get allele frequency differences
        y1afd_values <- data.frame(matrix(ncol=2,nrow=0))
        thenames = c('chr', 'AFD')
        colnames(y1afd_values) = thenames
        for(l in 1:length(y1afd.df)){
          if(!is.null(y1afd.df[[l]]))
            for(i in (y1afd.df[[l]])){
              y1afd_values[nrow(y1afd_values) + 1,] = c(l, i)
            }
        }
        
       #combine with peaks
        multi_table_one = cbind(multi_table_one, y1afd_values, rep = 1)
        #get rid of second chr column
        multi_table_one=multi_table_one[c("chr","LOD","max_Index","left_Index","right_Index","AFD","rep")]
        
        
        for (y in seq_along(multipeaks_two)) {
            if (!is.null(multipeaks_two[[y]]))
                for (q in seq_along(multipeaks_two[[y]][, 1])) {
                    ## has to return an index you can feed to roll_high_two
                    h_roll_chr <- roll_high_two[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                     x = names(roll_high_two))]
                    l_roll_chr <- roll_low_two[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                    x = names(roll_low_two))]
                    t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
                    t_min <- which.min(abs(multipeaks_two[[y]][q, "maxIndex"] - t_snp))
                    y2_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
                    y2afd.df[[y]][q] = y2_afd[[y]][q] #save for later to make heatmap with allele freq diffs
                    y2_afd[[y]][q] <- ifelse(y2_afd[[y]][q] > 0,
                                      ifelse(y2_afd[[y]][q] > AFThres,
                                             y2_afd[[y]][q] + 0.05,
                                             AFThres + 0.05),
                                      ifelse(y2_afd[[y]][q] < -AFThres,
                                             y2_afd[[y]][q] - 0.05,
                                             -AFThres - 0.05))
                }
        }

        #get allele frequency differences
        y2afd_values <- data.frame(matrix(ncol=2,nrow=0))
        thenames = c('chr', 'AFD')
        colnames(y2afd_values) = thenames
        for(l in 1:length(y2afd.df)){
          if(!is.null(y2afd.df[[l]]))
            for(i in (y2afd.df[[l]])){
              y2afd_values[nrow(y2afd_values) + 1,] = c(l, i)
            }
        }
        
        #combine with peaks
        multi_table_two = cbind(multi_table_two, y2afd_values, rep = 2)
        #get rid of second chr column
        multi_table_two = multi_table_two[c("chr","LOD","max_Index","left_Index","right_Index","AFD","rep")]
        
        #now combine peaks and AF diffs for reps 1 and 2
        forheatmap = rbind(multi_table_one, multi_table_two)
        #add environment and reporter as a column
        forheatmap$reporter = experiment_file[p, 3]
        forheatmap$environment = experiment_file[p, 6]
        
        
        ## write the table to a .txt file 
        write.table(x = forheatmap,
                    file = paste0(peaks_dir,
                                  experiment_file[p, 6], "_", #environment name
                                  experiment_file[p, 3], "_", #degron name
                                  experiment_file[p, 4], "_", #rep1 name
                                  experiment_file[p, 5], #rep2 name
                                  "_peaks_AFDs_bothrep.txt"),
                    row.names = F,
                    ## otherwise puts quotes in col names 
                    quote = F,
                    col.names = TRUE,
                    sep = '\t')
        
        
        ## make pdfs for plots
        
        pdf(file = paste0(results_dir,
                          experiment_file[p, 6], "_", #environment name
                          experiment_file[p, 3], "_", #degron name
                          experiment_file[p, 4], "_", #rep1 name
                          experiment_file[p, 5], #rep2 name
                          "_AFD_LOD_combined.pdf"),
            width=11,
            height=11)
        par(mfrow = c(2, 1))

        main1 = gsub(pattern = "_",
                     replacement = " ",
                     experiment_file[p, 3]) #degron name
        main2 = experiment_file[p, 6] #environment
        mytitle = paste(main1, "in", main2) 
        
        plot(gcoords,
             roll_high_one - roll_low_one,
             cex = 0.2,
             col = "grey",

             main = mytitle,
             type = "n",
             xaxt = 'n',
             xlab = "Chromosome",
             ylab = "High - Low Population",
             ylim = c(-0.75, 0.75))

        ## allele frequency difference thresholds
        ## these show the 99.99% quantile lines for
        ## allele freq. difference plots 
        abline(h = c(AFThres, -AFThres),
               lty = 2,
               lwd=2,
               col = sig_line_col)

        ## replicate one allele freq. difference lines
        for (j in unique(SNPs[,1])){
            points(gcoords[these_counts_one[,"chr"] == j],
                   (roll_high_one - roll_low_one)[these_counts_one[, "chr"] == j],
                   type = "l",
                   lwd = 2,
                   col = rep_one_col) ## dark grey
        }

        ## replicate two allele freq. difference lines
        for (j in unique(SNPs[, 1])) {
            points(gcoords[these_counts_two[, "chr"] == j],
                   (roll_high_two - roll_low_two)[these_counts_two[,"chr"] == j],
                   type = "l",
                   lwd = 2,
                   col = rep_two_col) ## light grey
        }

        ## chromosome divider lines
        for (j in c(chrCutoffs, 13566086)) {
            abline(v = j,
                   lty = 2,
                   col = "light blue")
        }

        ## 0 allele frequency difference line
        abline(h = 0,
               lty = 2,
               col = "light blue")

        axis(side = 1,
             at = chrLabels,
             labels = as.roman(1:16), tick=F)

        legend("topleft",
               legend = c("replicate 1", "replicate 2", "99.9% quantile"),
               lty = c(1, 1, 2), lwd = 2, seg.len = 3,
               col = c(rep_one_col, rep_two_col, sig_line_col),
               box.lty = 0, bty = "o", bg = "#ffffffaa")

        ## add stars at peaks
        for (j in 1:16) {
            if(!is.null(multipeaks_one[[j]])) {
                for (thisPeak in 1:nrow(multipeaks_one[[j]])){
                    thisPeakPlotPos <- getGcoords(j, 
                                                  multipeaks_one[[j]][thisPeak, "maxIndex"],
                                                  sepBetweenChr)

                    text(x = thisPeakPlotPos,
                         y = y1_afd[[j]][thisPeak],
                         labels = "*",
                         col = rep_one_col,
                         cex = 2.5)
                }
            }
        }

        for (j in 1:16) {
            if(!is.null(multipeaks_two[[j]])) {
                for (thisPeak in 1:nrow(multipeaks_two[[j]])){
                    thisPeakPlotPos <- getGcoords(j, 
                                                  multipeaks_two[[j]][thisPeak, "maxIndex"],
                                                  sepBetweenChr)

                    text(x = thisPeakPlotPos,
                         y = y2_afd[[j]][thisPeak],
                         labels = "*",
                         col = rep_two_col,
                         cex = 2.5)
                }
            }
        } 

        ## <<LOD_plot>>
        ## this produces and empty plot that we fill
        ## via the 'for' loop below 
        plot(gcoords, rep(0, length(gcoords)),
             type = "n",
             xaxt = "n",
             xlab = "Chromosome",
             ylab = "Multipool LOD",
             ylim = c(0, ylimMax))

        ## add the multipool trace and asterisks
        ## at peaks for the first replicate
        for (j in 1:16) {
            points(getGcoords(paste0("chr", as.roman(j)),
                              multipool_out_one[[j]][[2]][,1],
                              sepBetweenChr),
                   multipool_out_one[[j]][[2]][,2],
                   type = "l",
                   lwd = 2,
                   col = rep_one_col)

            ## add stars at peaks
            if(!is.null(multipeaks_one[[j]])) {
                for (thisPeak in 1:nrow(multipeaks_one[[j]])){
                    thisPeakPlotPos <- getGcoords(j,
                                                  multipeaks_one[[j]][thisPeak, "maxIndex"],
                                                  sepBetweenChr)

                    text(x = thisPeakPlotPos,
                         y = multipeaks_one[[j]][thisPeak, "maxValue"] + 1.25,
                         labels = "*",
                         col = rep_one_col,
                         cex = 2.5)
                }
            }
        }

        ## add the multipool trace and asterisks
        ## at peaks for the second replicate
        for (j in 1:16) {
            points(getGcoords(paste0("chr", as.roman(j)),
                              multipool_out_two[[j]][[2]][,1],
                              sepBetweenChr),
                   multipool_out_two[[j]][[2]][,2],
                   type = "l",
                   lwd = 2,
                   col = rep_two_col)

            ## add stars at peaks
            if(!is.null(multipeaks_two[[j]])) {
                for (thisPeak in 1:nrow(multipeaks_two[[j]])){
                    thisPeakPlotPos <- getGcoords(j,
                                                  multipeaks_two[[j]][thisPeak, "maxIndex"],
                                                  sepBetweenChr)

                    ## fix this so the asterisk height is proportional to peak height
                    text(x = thisPeakPlotPos,
                         y = multipeaks_two[[j]][thisPeak, "maxValue"] + 1.25,
                         labels = "*",
                         col = rep_two_col,
                         cex = 2.5)
                }
            }
        }

        ## non-loop, rest of the plot graphical elements
        ## 0 LOD line 
        abline(h = 0,
               lty = 2,
               col = "light blue")

        ## chromosome divider lines
        for (j in c(chrCutoffs, 13566086)) {
            abline(v = j,
                   lty = 2,
                   col = "light blue")
        }

        axis(side = 1,
             at = chrLabels,
             labels = as.roman(1:16),
             tick = F)

        abline(h = multiThres,
               lty = 2,
               col = sig_line_col,
               lwd = 2)

        legend("topleft",
               legend = c("replicate 1", "replicate 2", "significance threshold"),
               lty = c(1, 1, 2), lwd = 2, seg.len = 3,
               col = c(rep_one_col, rep_two_col, sig_line_col),
               box.lty = 0, bty = "o", bg = "#ffffffaa")

        dev.off()
}

