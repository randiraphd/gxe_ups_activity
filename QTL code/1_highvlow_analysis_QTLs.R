## -----
#Adapted by Randi Avery in 2024 from code by Mahlon Collins:
  #Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
  #Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570


# can then call this from external
i = 1

## trailing slash at the ends of all dirs!!!!
system   <- Sys.info()["nodename"]
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
#table with all .vcf files we are analyzing. the final two columns should NOT have '.vcf' on the end
comp_table     <- paste0(proj_dir, "2024.07.05_comparison_tableRRA_allGxE.txt")
#make a name of file you want for the output of the median coverage
medcovtable <- paste0(proj_dir, "coverage_table_24.07.05.txt")
#################
## END USER INPUT
#################
needed.dirs <- c("results/", "alignments/", "r_output/", "peaks/", "QTL_scripts/")

dir.maker <- function(x){if(!dir.exists(paths = paste0(proj_dir, x)))
                             dir.create(path = paste0(proj_dir, x))}

sapply(X = needed.dirs, FUN = dir.maker)

resultsFolder  <- paste0(proj_dir, "results/")
alignmentDir   <- paste0(proj_dir, "alignments/")
#put all .vcf files in this directory if they aren't already there
output_dir     <- paste0(proj_dir, "r_output/")
peaks_dir      <- paste0(proj_dir, "peaks/")
experimentFile <- read.table(comp_table, #reads in file with all the names of the .vcf to be analyzed
                             stringsAsFactors=FALSE,
                             head=TRUE)

## SNPs is a giant table w/ SNP positions
SNPs           <- read.table("~/myproteasome/2024.07.05_allHL_GxE/QTL_scripts/SNPs_Maggie_170809_BY_positions.txt",
                             stringsAsFactors=FALSE,
                             head=FALSE)

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

bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x,
                                                               INSTALL_opts = '--no-lock')}
bioc_package_installer("VariantAnnotation")
                       #, dependencies = TRUE)

library("VariantAnnotation")
source(paste0(proj_dir, "QTL_scripts/gTest.R"))
source(paste0(proj_dir, "QTL_scripts/x_qtl_seq_functions_170831.R"))
source(paste0(proj_dir, "QTL_scripts/mp_JB_170901.R"))
source(paste0(proj_dir, "QTL_scripts/peaksFromVector.R"))

## data frame with all yeast genes, plus
## chr, pos., strand, and names 
geneInfo <- read.table(paste0(proj_dir, "QTL_scripts/ensemblGenes_ensembl83_160307_MOD.txt"),
                       stringsAsFactors=FALSE,
                       sep="\t",
                       header=TRUE)

## rownames become systematic names 
rownames(geneInfo) <- geneInfo[,"geneID"]

## "geneName" is the common name, e.g., 'HOG1'
## for some (many?) rows of 'geneInfo', there
## is no 'geneName', so it's just an empty string
## e.g., head(allNames)
allNames <- geneInfo[, "geneName"]

names(allNames) <- geneInfo[,1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

sepBetweenChr <- 1e5
trimFromEnd   <- 15e3
obsMin        <- 10
LoessSpan     <- 0.1
## same as in Albert 2014
AFThres       <- 0.09653124
## multipool LOD threshold, our usual value 
multiThres    <- 4.5 

## this loop runs over the comparison table
## and generates all the output for the
## analysis, including plots, R objects,
## and peak tables 
for (i in (1:nrow(experimentFile))){

  #to run over just one file:
  #i = 16
  
    highname <- experimentFile[i,7]
    print(highname)
    highname <-substr(highname,1, (nchar(highname)-29)) #removes _S57_L001_sort_filtered_rmdup from the name
    lowname <- experimentFile[i,8]
    lowname <-substr(lowname,1, (nchar(lowname)-29))
    plotName <- paste0(c(highname, lowname),
                       collapse="_")
    mainPlotName <- paste0(c(highname, lowname),
                           collapse=" ")
    
    thisGene <- experimentFile[i, "gene"]


    highFile <- dir(alignmentDir,
                    pattern=paste0(experimentFile[i, "high"],
                                   ".*vcf$"),
                    full.names=TRUE)
    
    lowFile <- dir(alignmentDir,
                   pattern=paste0(experimentFile[i, "low"],
                                  ".*vcf$"),
                   full.names=TRUE)

    highVCF <- readVcf(highFile)
    lowVCF <- readVcf(lowFile)

    highVCFCounts <- t(sapply(info(highVCF)$AD,
                              function(x){c(x[[1]], x[[2]])}))
    
    lowVCFCounts <- t(sapply(info(lowVCF)$AD,
                             function(x){c(x[[1]], x[[2]])}))
    
    rownames(highVCFCounts) <- sapply(names(rowRanges(highVCF)),
                                      function(x){strsplit(x, "_")[[1]][1]})
    
    rownames(lowVCFCounts) <- sapply(names(rowRanges(lowVCF)),
                                     function(x){strsplit(x, "_")[[1]][1]})

    ## these counts were read from the vcf, so they
    ## don't contain counts for sites without reads
    ## need to build a common SNP table
    theseCounts <- data.frame(SNPs,
                              matrix(0,
                                     nrow=nrow(SNPs),
                                     ncol=4))
    rownames(theseCounts) <- paste(theseCounts[,1],
                                   theseCounts[,2],
                                   sep=":")

    colnames(theseCounts) <- c("chr",
                               "pos",
                               "high_ref",
                               "high_alt",
                               "low_ref",
                               "low_alt")
    
    theseCounts[rownames(highVCFCounts), c("high_ref", "high_alt")] <- highVCFCounts
    theseCounts[rownames(lowVCFCounts), c("low_ref", "low_alt")] <- lowVCFCounts

    gcoords= getGcoords(theseCounts$chr, theseCounts$pos, sepBetweenChr)
    names(gcoords) = rownames(theseCounts)

    ## get plot coordinates for the chromosome line dividers
    chrCutoffs <- sapply(unique(theseCounts$chr),
                         function(x){gcoords[theseCounts$chr == x][1] - sepBetweenChr/2})
    
    names(chrCutoffs) <- unique(theseCounts$chr)

    ## position of the text on the x axis for the chromosomes 
    chrLabels = sapply(1:(length(chrCutoffs)-1),
                       function(i)(chrCutoffs[i] + chrCutoffs[i+1])/2)

    ## add half the length of chrXVI
    chrLabels = c(chrLabels, chrCutoffs[16] + sepBetweenChr + 948066/2)
    names(chrLabels)[16] = "chrXVI"


    ## write out the counts to be analyzed
    setwd(output_dir)
    save(theseCounts,
         file=paste(paste0(experimentFile[i, c(7, 8)],
                           collapse="_"),
                    ".RData",
                    sep = ""))

# compute BY allele frequencies
theseCounts$highBYAF     <- theseCounts$high_ref / (theseCounts$high_ref + theseCounts$high_alt)
theseCounts$lowBYAF      <- theseCounts$low_ref  / (theseCounts$low_ref + theseCounts$low_alt)
theseCounts$highCoverage <- theseCounts$high_ref + theseCounts$high_alt
theseCounts$lowCoverage  <- theseCounts$low_ref + theseCounts$low_alt


    #saves median coverage of each file into a table
    medianCov = apply(theseCounts[,c("highCoverage", "lowCoverage")],
                      2,
                      median)
    medianCovTable = as.data.frame(medianCov)
    tablenames = c(highname, lowname)
    medianCovTable$filename = tablenames
    write.table(x = medianCovTable,
           file = medcovtable,
           sep = "\t",
           col.names=F,
           row.names=FALSE,
           append=T,
           quote = FALSE)

    ## at how many SNPs do we have the minimum number of counts? 
    SNPsAtMinObs = length(theseCounts[,"highCoverage"] >= obsMin & theseCounts[,"lowCoverage"] >= obsMin)



#######################
## <<plot_code>>
#######################

    #####################################################
    ## <<raw_allele_freq_plot>>
    ## this shows the allele frequency for the
    ## two populations you're comparing as separate
    ## plots.  it also marks the MAT, CAN1, and LYP1 loci
    #####################################################
    pdf(file = paste(resultsFolder, plotName, "_all_raw.pdf", sep=""),
        width=11,
        height=8)

    for(thisPop in c("high", "low")){
        plot(gcoords,
             rep(0.5, length(gcoords)),
             ylim=c(0,1),
             main = paste(mainPlotName,
                          thisPop,
                          collapse=" "),
             xaxt='n',
             xlab="chromosome",
             ylab="BY allele frequency",
             type="n")

        ## CAN1 is on chr05, 33466
        ## MATALPHA is on chr03, 198671
        ## LYP1 is on chr XIV, 140385
        abline(v = getGcoords("chrIII", 198671, sepBetweenChr),
               lwd = 2,
               col = "#4b2991")
        
        abline(v = getGcoords("chrV", 33466, sepBetweenChr),
               lwd = 2,
               col = "#c43282")
        
        abline(v = getGcoords(geneInfo[thisGene, "chr"],
                              mean(as.numeric(geneInfo[thisGene, c("start", "end")])),
                              sepBetweenChr),
               lwd = 2,
               col = "#f6a97a")

        points(gcoords, theseCounts[,paste0(thisPop, "BYAF", collapse="")],
               cex=.2,
               col="#00000022")
    
        roll = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"],
                                                    theseCounts[,paste0(thisPop, "BYAF", collapse="")],
                                                    gcoords,
                                                    median(theseCounts[,paste0(thisPop, "Coverage", collapse="")]),
                                                    stringsAsFactors=FALSE),
                                         LoessSpan)

        names(roll) = rownames(theseCounts)
    
    for (j in unique(theseCounts[,"chr"])){
        points(gcoords[theseCounts$chr == j],
               roll[theseCounts$chr == j],
               type="l",
               lwd=2)
    }

        ## make sure we get the end of XVI
        for (j in c(chrCutoffs, 13566086)){
            abline(v = j,
                   lty=2,
                   col="light blue")
        }
        
        abline(h = 0.5,
               lty=2,
               col="light blue")

        legend("topleft",
               legend = c("MAT", "CAN1", "LYP1"),
               lty = 1, lwd = 2, text.font = c(3, 3, 3),
               col = c("#4b2991", "#c43282", "#f6a97a"),
               box.lty = 0, bty = "o", bg = "#ffffffaa")

        legend("topright",
               legend = c(paste0("median coverage: ",
                               median(theseCounts[,paste0(thisPop,
                                                          "Coverage",
                                                          collapse="")]))),
               box.lty = 0, bty = "o", bg = "#ffffffaa")

        axis(side=1,
             at = chrLabels,
             labels = as.roman(1:16),
             tick=F)
}

dev.off()



    #########################################################
    ## <<combined_plot>>
    ## this plot shows the allele frequency trace for the
    ## two populations you're comparing.  it also marks three
    ## loci: (1) MAT, (2) CAN1, and (3) LYP1.  
    #########################################################
    pdf(file = paste(resultsFolder, plotName, "_all_combined.pdf", sep=""),
        width=11,
        height=8)

    plot(gcoords,
         rep(0.5, length(gcoords)),
         ylim=c(0,1),
         cex=.2,
         col="grey",
         main = mainPlotName,
         type="n",
         xaxt='n',
         xlab="chromosome",
         ylab="BY allele frequency")

        ## CAN1 is on chr05, 33466
        ## MATALPHA is on chr03, 198671
        ## LYP1 is on chr XIV, 140385
        abline(v = getGcoords("chrIII", 198671, sepBetweenChr),
               lwd = 2,
               col = "#4b2991")
        
        abline(v = getGcoords("chrV", 33466, sepBetweenChr),
               lwd = 2,
               col = "#c43282")
        
        abline(v = getGcoords(geneInfo[thisGene, "chr"],
                              mean(as.numeric(geneInfo[thisGene, c("start", "end")])),
                              sepBetweenChr),
               lwd = 2,
               col = "#f6a97a")

    cols <- cbind(c("#6f125fAA", "#202020CC"),
                  c("#7f126f22", "#70707022"))
    
    rownames(cols) <- c("high", "low")
    colnames(cols) <- c("lines", "points")

    for (thisPop in c("high", "low")){
        points(gcoords, theseCounts[,paste0(thisPop, "BYAF", collapse="")],
               cex=.2,
               col=cols[thisPop, "points"])
    }

    for(thisPop in c("high", "low")){
        roll = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"],
                                                    theseCounts[,paste0(thisPop, "BYAF", collapse="")],
                                                    gcoords, median(theseCounts[,paste0(thisPop, "Coverage", collapse="")]),
                                                    stringsAsFactors=FALSE),
                                         LoessSpan)
        
        for (j in unique(theseCounts[,"chr"])){
            points(gcoords[theseCounts[,"chr"] == j],
                   roll[theseCounts[,"chr"] == j],
                   type = "l",
                   lwd = 2,
                   col = cols[thisPop, "lines"])
    }
}

for (j in c(chrCutoffs, 13566086)){
    abline(v = j,
           lty = 2,
           col = "light blue")
}

    abline(h = 0.5,
           lty = 2,
           col = "light blue")

    ## high vs. low pool legend
    legend("topright",
           legend = c("high pool", "low pool"),
           col = cols, lty = 1, lwd = 2, box.lty = 0,
           bty = "o", bg = "#ffffffaa")

    ## marker loci legend
    legend("topleft",
           legend = c("MAT", "CAN1", "LYP1"),
           lty = 1, lwd = 2, text.font = c(3, 3, 3),
           col = c("#4b2991", "#c43282", "#f6a97a"),
           box.lty = 0,
           bty = "o", bg = "#ffffffaa")

    axis(side=1, at = chrLabels, labels = as.roman(1:16), tick=F)

dev.off()

##################################################
## <<difference_plot>>
## this plot shows the difference in allele
## frequency between the two populations compared.
## the expectation is 0 for a locus w/ no effect
## on the trait you're measuring.  the red ticked
## line shows significant deflections above chance
##################################################
    pdf(file = paste(resultsFolder, plotName, "_difference.pdf", sep=""),
        width=11,
        height=8)

    rollHigh = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"],
                                                    theseCounts[,"highBYAF"],
                                                    gcoords,
                                                    median(theseCounts[,"highCoverage"]),
                                                    stringsAsFactors=FALSE),
                                         LoessSpan)

    rollLow = rollLoessByChrWithWeights(data.frame(theseCounts[,"chr"],
                                                   theseCounts[,"lowBYAF"],
                                                   gcoords,
                                                   median(theseCounts[,"lowCoverage"]),
                                                   stringsAsFactors=FALSE),
                                        LoessSpan)


    plot(gcoords,
         rollHigh - rollLow,
         cex=0.2,
         col="grey",
         main = mainPlotName,
         type="n",
         xaxt='n',
         xlab="chromosome",
         ylab="High - Low Population",
         ylim=c(-1,1))

    points(gcoords,
           theseCounts[,"highBYAF"] - theseCounts[,"lowBYAF"],
           cex=.2,
           col="#00000022")

    ## allele frequency difference thresholds
    ## these show the significance line for
    ## allele freq. difference plots 
    abline(h = AFThres,
           lty = 2,
           lwd=2,
           col = "red")
    
    abline(h = -AFThres,
           lty = 2,
           lwd=2,
           col = "red")

for (j in unique(SNPs[,1])){
    points(gcoords[theseCounts[,"chr"] == j],
           (rollHigh - rollLow)[theseCounts[,"chr"] == j],
           type="l",
           lwd=2,
           col="black")
}

for (j in c(chrCutoffs, 13566086)){
    abline(v = j,
           lty=2,
           col="light blue")
}
    
    abline(h = 0,
           lty=2,
           col="light blue")

    axis(side=1,
         at = chrLabels,
         labels = as.roman(1:16), tick=F)

dev.off()


############
## multipool
############
# make a switch to kill plotting etc here if no multipool desired

# run multipool
# per chromosome!


    #this takes a couple minutes on one replicate
    multipoolOutput <- lapply(unique(theseCounts[,"chr"]),
                              function(j){
                                  doMultiPoolFromWithinR(theseCounts[theseCounts$chr == j,
                                                                     c("pos", "high_ref", "high_alt")],
                                                         theseCounts[theseCounts$chr == j,
                                                                     c("pos", "low_ref", "low_alt")])
                              })


## call peaks from multipool LODs
multiPeaks <- lapply(multipoolOutput, function(x){
    thisLODTrace = x[[2]]
    theseChrPeaks = callPeaks(thisLODTrace[,2], 4.5, 2)
    theseChrPeaks
})

    setwd(output_dir)
    save(multiPeaks,
         multipoolOutput,
         file=paste(paste0(experimentFile[i, c(7, 8)],
                           collapse="_"),
                    "_multipoolResults.RData",
                    sep=""))

    ## write the multipool peaks to a .txt table
    ## this converts the multiPeaks list to a dataframe
    ## but doesn't provide chromosome number 
    multi_table <- do.call(rbind.data.frame, multiPeaks)

    ## get the chromosome number by looping over the
    ## multiPeaks list and repeating the chromosome number
    ## for each QTL in a chromosome.  if there are no QTLs
    ## on a chromosome, do nothing
    multi_chr_number <- vector()
    for(l in 1:length(multiPeaks)){
        if(!is.null(multiPeaks[[l]]))
            multi_chr_number <- c(multi_chr_number,
                                  rep(l, nrow(multiPeaks[[l]])))
    }

    multi_table$chr <- multi_chr_number

    ## these need to be 'if' clauses to account for
    ## situations where there are no multiPool peaks, #rra, but need to also allow when there is only 1 peak
    ## re-order the output 
    if(nrow(multi_table) >= 1)
        multi_table <- multi_table[, c(5, 2, 1, 3, 4)]
    
    ## re-name a couple of columns for easier reading
    if(nrow(multi_table) >= 1)
        colnames(multi_table) <- c("chr", "LOD", "max_Index",
                                   "left_Index", "right_Index")

    ## write the table to a .txt file 
    write.table(x = multi_table,
                file = paste(peaks_dir,
                    paste0(experimentFile[i, c(7, 8)],
                           collapse="_"),
                           "_multipool_peaks.txt",
                           sep = ""),
                row.names = F,
                ## otherwise puts quotes in col names 
                quote = F)

################  
## <<LOD_plots>>
################
## this code produces a two-panel plot
## of the allele frequency difference between
## the high and low populations, as well as the
## significant QTLs.  The QTLs are marked w/ 
## asterisks in the bottom plot.  I no longer  
## plot LYP1, since it's not really relevant
## and it's somewhat visually distracting 

    pdf(paste(resultsFolder, plotName, "_difference_withMultipool.pdf", sep=""),
        width=11,
        height=11)

    par(mfrow=c(2,1))

    plot(gcoords,
         rep(0, length(gcoords)),
         cex=.2, col="grey",
         main = mainPlotName,
         type="n",
         xaxt='n',
         xlab="chromosome",
         ylab="High - Low Population",
         ylim=c(-1,1))
    
    points(gcoords,
           theseCounts[,"highBYAF"] - theseCounts[,"lowBYAF"],
           cex=.2,
           col="#00000022")


    ## difference thresholds from the Albert 2014 nullSorts:
    abline(h = c(AFThres, -AFThres),
           lty = 2,
           lwd = 2,
           col = "red")

    for (j in unique(SNPs[,1])){
        points(gcoords[theseCounts[,"chr"] == j],
               (rollHigh - rollLow)[theseCounts[,"chr"] == j],
               type="l",
               lwd=2,
               col="black")
    }

    ## chromosome separators
    for (j in c(chrCutoffs, 13566086)){
    abline(v = j,
           lty=2,
           col="light blue")
    }

    ## 0 allele frequency difference line
    abline(h = 0,
           lty=2,
           col="light blue")

    axis(side = 1,
         at = chrLabels,
         labels = as.roman(1:16),
         tick = F)

    ylimMax = max(c(multiThres, sapply(multipoolOutput, function(x){max(x[[2]][,2])}))) + 1

    plot(gcoords, rep(0, length(gcoords)),
         main = mainPlotName,
         type="n",
         xaxt='n',
         xlab="chromosome",
         ylab="Multipool LOD",
         ylim=c(0, ylimMax))


    for (j in 1:16){
        points(getGcoords(paste0("chr", as.roman(j)),
                          multipoolOutput[[j]][[2]][,1],
                          sepBetweenChr),
           multipoolOutput[[j]][[2]][,2],
           type="l",
           lwd=2,
           col="black")
    
    ## add stars at peaks
    if(!is.null(multiPeaks[[j]])){
        for (thisPeak in 1:nrow(multiPeaks[[j]])){
            thisPeakPlotPos <- getGcoords(j,
                                          multiPeaks[[j]][thisPeak, "maxIndex"],
                                          sepBetweenChr)

            ## fix this so the asterisk height is proportional to peak height
            text(x = thisPeakPlotPos,
                 y = multiPeaks[[j]][thisPeak, "maxValue"] + 1.25,
                 labels="*",
                 col="#aa1111",
                 cex=2.5)
        }
    }

    }
    ## 0 LOD line 
    abline(h = 0,
           lty = 2,
           col="light blue")

    ## chromosome divider lines
    for (j in c(chrCutoffs, 13566086)){
    abline(v = j,
           lty = 2,
           col="light blue")
    }

    axis(side=1,
         at = chrLabels,
         labels = as.roman(1:16),
         tick=F)

    abline(h = multiThres,
           lty = 2,
           col = "red",
           lwd = 2)

dev.off()

}
