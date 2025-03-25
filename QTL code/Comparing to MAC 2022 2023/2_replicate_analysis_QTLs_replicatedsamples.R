#Adapted by Randi Avery in 2024 from code by Mahlon Collins:
  #Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
  #Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570

#combines biological replicates of QTLs, makes QTL traces and LOD plots
#Produces Supplementary Figure 2A-D

base_dir <- "~/myproteasome/mac_compare/"


#############
## USER INPUT
#############
## set the specific directory you'll work 
## in and name the comparison table
## TRAILING SLASH AT END OF DIR
## below, your project, e.g.,
## "2020.08.17_FPFA002_TDH3pr_Arg_N-end_TFT_sorts/"
proj           <- "24.07.18_maccompare/"
proj_dir       <- paste0(base_dir, proj)
#copy and paste the file names from the previous output into the replicates comparison table
comp_table     <- paste0(proj_dir, "2024.07.28_4reps_comparison_table_maccompare.txt")

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


#edited for 4 replicates
library(scales) #for getting alpha for base R plot

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
    
    multi_table_one$rep_name = experiment_file[p, 6]
    
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
    
    multi_table_two$rep_name = experiment_file[p, 7]
#####PASTED STARTING HERE for reps 3 and 4
    
    ## load 'multipoolResults.RData' for replicate 1
    ## adds 'multiPeaks' and 'multipoolOutput'
    load(paste0(data_dir,
                experiment_file[p, 3],
                mpr, rd))
    print(experiment_file[p, 3])
    ## load ".RData" for replicate 1
    ## adds 'theseCounts'
    load(paste0(data_dir,
                experiment_file[p, 3], 
                rd))
    
    multipeaks_three    <- multiPeaks
    multipool_out_three <- multipoolOutput
    these_counts_three  <- theseCounts
    
    #get peak information again to include allele frequency difference this time for the peaks
    multi_table_three <- do.call(rbind.data.frame, multipeaks_three)
    
    ## get the chromosome number by looping over the
    ## multiPeaks list and repeating the chromosome number
    ## for each QTL in a chromosome.  if there are no QTLs
    ## on a chromosome, do nothing
    multi_chr_number_three <- vector()
    for(l in 1:length(multipeaks_three)){
      if(!is.null(multipeaks_three[[l]]))
        multi_chr_number_three <- c(multi_chr_number_three,
                                  rep(l, nrow(multipeaks_three[[l]])))
    }
    
    multi_table_three$chr <- multi_chr_number_three
    
    ## these need to be 'if' clauses to account for
    ## situations where there are no multiPool peaks
    ## re-order the output 
    if(nrow(multi_table_three) >= 1)
      multi_table_three <- multi_table_three[, c(5, 2, 1, 3, 4)]
    
    ## re-name a couple of columns for easier reading
    if(nrow(multi_table_three) >= 1)
      colnames(multi_table_three) <- c("chr", "LOD", "max_Index",
                                     "left_Index", "right_Index")
    
    multi_table_three$rep_name = experiment_file[p, 8]
    ## load 'multipoolResults.RData' for replicate 2
    ## adds 'multiPeaks' and 'multipoolOutput'
    load(paste0(data_dir,
                experiment_file[p, 4],
                mpr, rd))
    
    ## load ".RData" for replicate 2
    ## adds 'theseCounts'
    load(paste0(data_dir,
                experiment_file[p, 4],
                rd))
    
    multipeaks_four    <- multiPeaks
    multipool_out_four <- multipoolOutput
    these_counts_four  <- theseCounts
    
    #get peak information again to include allele frequency difference this time for the peaks
    multi_table_four <- do.call(rbind.data.frame, multipeaks_four)
    
    ## get the chromosome number by looping over the
    ## multiPeaks list and repeating the chromosome number
    ## for each QTL in a chromosome.  if there are no QTLs
    ## on a chromosome, do nothing
    multi_chr_number_four <- vector()
    for(l in 1:length(multipeaks_four)){
      if(!is.null(multipeaks_four[[l]]))
        multi_chr_number_four <- c(multi_chr_number_four,
                                  rep(l, nrow(multipeaks_four[[l]])))
    }
    
    multi_table_four$chr <- multi_chr_number_four
    
    ## these need to be 'if' clauses to account for
    ## situations where there are no multiPool peaks, #rra, but need to also allow when there is only 1 peak
    ## re-order the output 
    if(nrow(multi_table_four) >= 1)
      multi_table_four <- multi_table_four[, c(5, 2, 1, 3, 4)]
    
    ## re-name a couple of columns for easier reading
    if(nrow(multi_table_four) >= 1)
      colnames(multi_table_four) <- c("chr", "LOD", "max_Index",
                                     "left_Index", "right_Index")
    
    multi_table_four$rep_name = experiment_file[p, 9]
    ####PASTED UNTIL HERE for reps 3 and 4
    
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
        
        these_counts_three$h_BY_AF <- allele_f_calc(these_counts_three, "high")
        these_counts_three$l_BY_AF <- allele_f_calc(these_counts_three, "low")
        
        these_counts_four$h_BY_AF <- allele_f_calc(these_counts_four, "high")
        these_counts_four$l_BY_AF <- allele_f_calc(these_counts_four, "low")
        
        
        these_counts_one$h_cover <- cover_calc(these_counts_one, "high")
        these_counts_one$l_cover <- cover_calc(these_counts_one, "low")

        these_counts_two$h_cover <- cover_calc(these_counts_two, "high")
        these_counts_two$l_cover <- cover_calc(these_counts_two, "low")
        
        these_counts_three$h_cover <- cover_calc(these_counts_three, "high")
        these_counts_three$l_cover <- cover_calc(these_counts_three, "low")
        
        these_counts_four$h_cover <- cover_calc(these_counts_four, "high")
        these_counts_four$l_cover <- cover_calc(these_counts_four, "low")        

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
                                  multipool_out_two,
                                  multipool_out_three,
                                  multipool_out_four),
                                function(x) {
                                    max(x[[2]][, 2])}))) + 1

        rep_one_col  <- "#3A3A3A"     ## dark gray (back)
        rep_two_col  <- "#8A8A8A"     ## light gray (front)
        rep_three_col   = '#004AFF' #blue
        rep_four_col = '#568AFF' #light blue
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
        #start paste here
        
        roll_high_three <- rollLoessByChrWithWeights(data.frame(these_counts_three[, "chr"],
                                                              these_counts_three[, "h_BY_AF"],
                                                              gcoords,
                                                              median(these_counts_three[, "h_cover"]),
                                                              stringsAsFactors = FALSE),
                                                   LoessSpan)
        
        roll_low_three <- rollLoessByChrWithWeights(data.frame(these_counts_three[, "chr"],
                                                             these_counts_three[, "l_BY_AF"],
                                                             gcoords,
                                                             median(these_counts_three[, "l_cover"]),
                                                             stringsAsFactors = FALSE),
                                                  LoessSpan)
        
        roll_high_four <- rollLoessByChrWithWeights(data.frame(these_counts_four[, "chr"],
                                                              these_counts_four[, "h_BY_AF"],
                                                              gcoords,
                                                              median(these_counts_four[, "h_cover"]),
                                                              stringsAsFactors = FALSE),
                                                   LoessSpan)
        
        roll_low_four <- rollLoessByChrWithWeights(data.frame(these_counts_four[, "chr"],
                                                             these_counts_four[, "l_BY_AF"],
                                                             gcoords,
                                                             median(these_counts_four[, "l_cover"]),
                                                             stringsAsFactors = FALSE),
                                                  LoessSpan)

        #end paste here
        
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
        y3_afd <- y1_afd
        y4_afd <- y1_afd
        y1afd.df = y1_afd #dfs to save allele frequencies before altering them for star placement
        y2afd.df = y1_afd
        y3afd.df = y1_afd
        y4afd.df = y1_afd
        
        
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
        
        names(y3_afd) <- chr_indices
        names(y4_afd) <- chr_indices
        names(y3afd.df) <- chr_indices
        names(y4afd.df) <- chr_indices
        
        
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
        multi_table_one=multi_table_one[c("chr","LOD","max_Index","left_Index","right_Index","AFD","rep", "rep_name")]
       

        
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
        multi_table_two = multi_table_two[c("chr","LOD","max_Index","left_Index","right_Index","AFD","rep", "rep_name")]
        
        
        
        
        ####START PASTE HERE       
        
        for (y in seq_along(multipeaks_three)) {
          if (!is.null(multipeaks_three[[y]]))
            for (q in seq_along(multipeaks_three[[y]][, 1])) {
              ## has to return an index you can feed to roll_high_three
              h_roll_chr <- roll_high_three[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                               x = names(roll_high_three))]
              l_roll_chr <- roll_low_three[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                              x = names(roll_low_three))]
              t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
              t_min <- which.min(abs(multipeaks_three[[y]][q, "maxIndex"] - t_snp))
              y3_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
              y3afd.df[[y]][q] = y3_afd[[y]][q] #save for later to make heatmap with allele freq diffs
              #if there should be a star, add or subtract 0.05 so it goes right above/below the peak
              y3_afd[[y]][q] <- ifelse(y3_afd[[y]][q] > 0,
                                       ifelse(y3_afd[[y]][q] > AFThres,
                                              y3_afd[[y]][q] + 0.05,
                                              AFThres + 0.05),
                                       ifelse(y3_afd[[y]][q] < -AFThres,
                                              y3_afd[[y]][q] - 0.05,
                                              -AFThres - 0.05))
            }
        }
        
        #get allele frequency differences
        y3afd_values <- data.frame(matrix(ncol=2,nrow=0))
        thenames = c('chr', 'AFD')
        colnames(y3afd_values) = thenames
        for(l in 1:length(y3afd.df)){
          if(!is.null(y3afd.df[[l]]))
            for(i in (y3afd.df[[l]])){
              y3afd_values[nrow(y3afd_values) + 1,] = c(l, i)
            }
        }
        
        #combine with peaks
        multi_table_three = cbind(multi_table_three, y3afd_values, rep = 1)
        #get rid of second chr column
        multi_table_three=multi_table_three[c("chr","LOD","max_Index","left_Index","right_Index","AFD","rep", "rep_name")]
        
        
        
        for (y in seq_along(multipeaks_four)) {
          if (!is.null(multipeaks_four[[y]]))
            for (q in seq_along(multipeaks_four[[y]][, 1])) {
              ## has to return an index you can feed to roll_high_four
              h_roll_chr <- roll_high_four[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                               x = names(roll_high_four))]
              l_roll_chr <- roll_low_four[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                              x = names(roll_low_four))]
              t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
              t_min <- which.min(abs(multipeaks_four[[y]][q, "maxIndex"] - t_snp))
              y4_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
              y4afd.df[[y]][q] = y4_afd[[y]][q] #save for later to make heatmap with allele freq diffs
              y4_afd[[y]][q] <- ifelse(y4_afd[[y]][q] > 0,
                                       ifelse(y4_afd[[y]][q] > AFThres,
                                              y4_afd[[y]][q] + 0.05,
                                              AFThres + 0.05),
                                       ifelse(y4_afd[[y]][q] < -AFThres,
                                              y4_afd[[y]][q] - 0.05,
                                              -AFThres - 0.05))
            }
        }
        
        #get allele frequency differences
        y4afd_values <- data.frame(matrix(ncol=2,nrow=0))
        thenames = c('chr', 'AFD')
        colnames(y4afd_values) = thenames
        for(l in 1:length(y4afd.df)){
          if(!is.null(y4afd.df[[l]]))
            for(i in (y4afd.df[[l]])){
              y4afd_values[nrow(y4afd_values) + 1,] = c(l, i)
            }
        }
        
        #combine with peaks
        multi_table_four = cbind(multi_table_four, y4afd_values, rep = 2)
        #get rid of second chr column
        multi_table_four = multi_table_four[c("chr","LOD","max_Index","left_Index","right_Index","AFD","rep", "rep_name")]

        
        ####END PASTE HERE
        

        
        #now combine peaks and AF diffs for reps 1 and 2, 3 and 5
        forheatmap = rbind(multi_table_one, multi_table_two, multi_table_three, multi_table_four)
        #add environment and reporter as a column
        forheatmap$reporter = experiment_file[p, 5]
        forheatmap$environment = experiment_file[p, 10]
        
        
        
        ## write the table to a .txt file 
        write.table(x = forheatmap,
                    file = paste0(peaks_dir,
                                  experiment_file[p, 10], "_", #environment name
                                  experiment_file[p, 5], "_", #degron name
                                  experiment_file[p, 6], "_", #rep1 name
                                  experiment_file[p, 7], "_", #rep2 name
                                  experiment_file[p, 8], "_", #rep3 name
                                  experiment_file[p, 9], #rep4 name
                                  "_peaks_AFDs_4reps_12.14.24.txt"),
                    row.names = F,
                    ## otherwise puts quotes in col names 
                    quote = F,
                    col.names = TRUE,
                    sep = '\t')
        
        
        ## make pdfs for plots
        
        pdf(file = paste0(results_dir,
                          experiment_file[p, 10], "_", #environment name
                          experiment_file[p, 5], "_", #degron name
                          experiment_file[p, 6], "_", #rep1 name
                          experiment_file[p, 7], "_", #rep2 name
                          experiment_file[p, 8], "_", #rep3 name
                          experiment_file[p, 9], #rep4 name
                          "_AFD_LOD_combined_12.14.24.pdf"),
            width=11,
            height=11)
        par(mfrow = c(2, 1))

        main1 = gsub(pattern = "_",
                     replacement = " ",
                     experiment_file[p, 5]) #degron name
        main2 = experiment_file[p, 10] #environment
        mytitle = paste(main1, "in", main2) 
        
        plot(gcoords,
             roll_high_one - roll_low_one, #only first replicate?
             cex = 0.2,
             col = "grey",

             main = mytitle,
             type = "n",
             xaxt = 'n',
             xlab = "Chromosome",
             ylab = "∆RM Allele Frequency (High - Low UPS Activity Pool)",
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

        #start paste
        
        ## replicate three allele freq. difference lines
        for (j in unique(SNPs[,1])){
          points(gcoords[these_counts_three[,"chr"] == j],
                 (roll_high_three - roll_low_three)[these_counts_three[, "chr"] == j],
                 type = "l",
                 lwd = 1.5, #made smaller
                 col = alpha(rep_three_col, 0.5)) ## 
        }
        
        ## replicate four allele freq. difference lines
        for (j in unique(SNPs[, 1])) {
          points(gcoords[these_counts_four[, "chr"] == j],
                 (roll_high_four - roll_low_four)[these_counts_four[,"chr"] == j],
                 type = "l",
                 lwd = 1.5,
                 col = alpha(rep_four_col, 0.5)) ##
        }
        
        
        #end paste
        
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
               legend = c("replicate 1", "replicate 2", "replicate 1 M.A. Collins et al. 2022", "replicate 2 M.A. Collins et al. 2022", "99.9% quantile"),
               lty = c(1, 1, 1, 1, 2), lwd = 2, seg.len = 3,
               col = c(rep_one_col, rep_two_col, alpha(rep_three_col, 0.5), alpha(rep_four_col, 0.5), sig_line_col),
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

        #start paste
        
        ## add stars at peaks
        for (j in 1:16) {
          if(!is.null(multipeaks_three[[j]])) {
            for (thisPeak in 1:nrow(multipeaks_three[[j]])){
              thisPeakPlotPos <- getGcoords(j, 
                                            multipeaks_three[[j]][thisPeak, "maxIndex"],
                                            sepBetweenChr)
              
              text(x = thisPeakPlotPos,
                   y = y3_afd[[j]][thisPeak],
                   labels = "*",
                   col = alpha(rep_three_col, 0.5),
                   cex = 2.5)
            }
          }
        }
        
        for (j in 1:16) {
          if(!is.null(multipeaks_four[[j]])) {
            for (thisPeak in 1:nrow(multipeaks_four[[j]])){
              thisPeakPlotPos <- getGcoords(j, 
                                            multipeaks_four[[j]][thisPeak, "maxIndex"],
                                            sepBetweenChr)
              
              text(x = thisPeakPlotPos,
                   y = y4_afd[[j]][thisPeak],
                   labels = "*",
                   col = alpha(rep_four_col, 0.5),
                   cex = 2.5)
            }
          }
        } 
        
        
        #end paste
        
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
        
        #start paste
        
        
        ## add the multipool trace and asterisks
        ## at peaks for the third replicate
        for (j in 1:16) {
          points(getGcoords(paste0("chr", as.roman(j)),
                            multipool_out_three[[j]][[2]][,1],
                            sepBetweenChr),
                 multipool_out_three[[j]][[2]][,2],
                 type = "l",
                 lwd = 1.5, #made smaller
                 col = alpha(rep_three_col, 0.5))
          
          ## add stars at peaks
          if(!is.null(multipeaks_three[[j]])) {
            for (thisPeak in 1:nrow(multipeaks_three[[j]])){
              thisPeakPlotPos <- getGcoords(j,
                                            multipeaks_three[[j]][thisPeak, "maxIndex"],
                                            sepBetweenChr)
              
              text(x = thisPeakPlotPos,
                   y = multipeaks_three[[j]][thisPeak, "maxValue"] + 1.25,
                   labels = "*",
                   col = alpha(rep_three_col),
                   cex = 2.5)
            }
          }
        }
        
        ## add the multipool trace and asterisks
        ## at peaks for the fourth replicate
        for (j in 1:16) {
          points(getGcoords(paste0("chr", as.roman(j)),
                            multipool_out_four[[j]][[2]][,1],
                            sepBetweenChr),
                 multipool_out_four[[j]][[2]][,2],
                 type = "l",
                 lwd = 1.5, #made smaller
                 col = alpha(rep_four_col, 0.5))
          
          ## add stars at peaks
          if(!is.null(multipeaks_four[[j]])) {
            for (thisPeak in 1:nrow(multipeaks_four[[j]])){
              thisPeakPlotPos <- getGcoords(j,
                                            multipeaks_four[[j]][thisPeak, "maxIndex"],
                                            sepBetweenChr)
              
              ## fix this so the asterisk height is proportional to peak height
              text(x = thisPeakPlotPos,
                   y = multipeaks_four[[j]][thisPeak, "maxValue"] + 1.25,
                   labels = "*",
                   col = alpha(rep_four_col, 0.5),
                   cex = 2.5)
            }
          }
        }
        
        
        
        #end paste

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
               legend = c("replicate 1", "replicate 2", "replicate 1 M.A. Collins et al. 2022", "replicate 2 M.A. Collins et al. 2022", "significance threshold"),
               lty = c(1, 1, 1, 1, 2), lwd = 2, seg.len = 3,
               col = c(rep_one_col, rep_two_col, alpha(rep_three_col, 0.5), alpha(rep_four_col, 0.5), sig_line_col),
               box.lty = 0, bty = "o", bg = "#ffffffaa")

        dev.off()
}

