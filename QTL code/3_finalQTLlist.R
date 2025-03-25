#Code by Randi Avery for 
  #Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity. bioRxiv, 2024.11.21.624644. https://doi.org/10.1101/2024.11.21.624644

#combines peak tables from each sample into one df
#saves all peaks into a table, separated by replicate
#combines QTLs found in both replicates, but keeps QTLs present in only replicate for later analysis
#removes peaks that were used more than once
#averages replicates (for QTLs present in only one replicate, just keeps those values)
#this results in a final list of QTLs to use for further analyses and plots

base_dir <- "~/myproteasome/"
proj  <- "2024.07.05_allHL_GxE/"
proj_dir <- paste0(base_dir, proj, "peaks/")
setwd(proj_dir)


#make list of files with peak info
peaks = list.files(pattern = 'AFDs_bothrep.txt$')

#read in files and combine into one df
allpeaks = do.call(rbind, lapply(peaks, function(x) read.table(x, header = TRUE)))

#find any rows that have an NA and send that to a dataframe to look at later. Will skip those in the code below
na_table = allpeaks[rowSums(is.na(allpeaks)) > 0, ] 
allpeaks = na.omit(allpeaks)

#save peaks file
write.table(allpeaks, #should be 1150 peaks
            file = "all_peaks_useforheatmap24.07.12.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#find QTLs that are present in both replicates, but also keep in only in one replicate for later analyses - just note this
  # same chromosome
  # peaks w/in 100 kb of ea. other
  # same direction of allele frequency difference (positive or negative)

#reorder columns to match the next steps (which I wrote previously)
allpeaks = allpeaks[,c('reporter', 'environment', 'chr', 'LOD', 'AFD', 'left_Index', 'max_Index', 'right_Index', 'rep')]

#make empty data frame to be used for overlaps
#want to keep QTLs that are only in one rep
thecolnames = c('reporter', 'environment', 'chrm', 'rep_1_LOD', 'rep_1_delta_AF', 'rep_1_left_Index', 'rep_1_max_Index', 'rep_1_right_Index', 'replicate1',
                'reporter', 'environment', 'chrm', 'rep_2_LOD', 'rep_2_delta_AF', 'rep_2_left_Index', 'rep_2_max_Index', 'rep_2_right_Index', 'replicate2', 
                'number_of_reps')
merge_frame = data.frame(matrix(ncol = 19, nrow = 0))
colnames(merge_frame) = thecolnames

reporters = unique(allpeaks$reporter) #get list of all the reporters

#for loop to combine reps. Keep QTLs present in only one rep
for(repor in reporters) {
  #repor= 'Thr' #just for testing one reporter at a time
  sub_reporters = allpeaks[allpeaks$reporter == repor,] #subset to get unique reporters
  for(e in unique(sub_reporters$environment)){ #now subset the reporter by the environments
    #print(e) #used to check
    #e = 'SC' #for testing loop
    sub_environ = sub_reporters[sub_reporters$environment == e,] #this is the data for one reporter in one environment
    rep1 = sub_environ[sub_environ$rep == 1,] 
    rep2 = sub_environ[sub_environ$rep == 2,]
    #subset chromosomes
    #i = 15 #for testing loop, chromosome number
    for(i in 1:16) { 
      subset_rep_1 <- rep1[rep1$chr == i, ]
      subset_rep_2 <- rep2[rep2$chr == i, ]
      #j = 2 #for testing loop
      usedQTLvec = 0 #keep track of which QTLs are used from rep 2
      #check to see if both reps are empty
      if (dim(subset_rep_1)[1]==0 & dim(subset_rep_2)[1]==0) {
        print(c(i,e, repor, ' no QTLs on this chrm'))
      }else {
        #check to see if rep 1 is empty
        if (dim(subset_rep_1)[1]==0) {#add QTLs from rep2 to merge_frame
          df = cbind(subset_rep_2,subset_rep_2)
          df$number_of_reps = 1 #only found once
          colnames(df) = thecolnames
          merge_frame = rbind(merge_frame, df)
        }else {
          for(j in 1:nrow(subset_rep_1)) {
            
            ## find peaks w / in 100 kb of ea. other:
            ## make sure to use 'abs' here!
            if(min(abs(subset_rep_1[j,"max_Index"] - subset_rep_2$max_Index)) <= 1e5) { 
              ## if we find something w / in 100 kb, grab that row from replicate 2
              subset_2_QTL <- subset_rep_2[which.min(abs(subset_rep_1[j, "max_Index"] - subset_rep_2$max_Index)), ] 
              usedQTLrep2 = which.min(abs(subset_rep_1[j, "max_Index"] - subset_rep_2$max_Index))
              ## make sure they have the same direction of effect:
              if(subset_rep_1[j,'AFD'] * subset_2_QTL$AFD > 0) {
                ## use 'cbind' to make a dataframe row that contains the info from both replicates
                df = cbind(subset_rep_1[j, ], subset_2_QTL)
                df$number_of_reps = 2
                colnames(df) = thecolnames
                merge_frame = rbind(merge_frame, df) #
                usedQTLvec = c(usedQTLvec, usedQTLrep2)
              } else {#if there is a peak within 100 kb, but it's going in the opposite direction
                df = cbind(subset_rep_1[j, ],subset_rep_1[j, ]) #just fill out table twice
                df$number_of_reps = 1 #only found once
                colnames(df) = thecolnames
                merge_frame = rbind(merge_frame, df)
              } #else 2 #warnings should just be from finding min when there is an empty vector
            } else { #if there isn't a peak within 100kb - keep QTL but say only found once
              df = cbind(subset_rep_1[j, ],subset_rep_1[j, ]) #just fill out table twice
              df$number_of_reps = 1 #only found once
              colnames(df) = thecolnames
              merge_frame = rbind(merge_frame, df)
            }#else 1
            
          } 
        }#test per chrom here
        #when we're done with a chromosome, add leftover rows from rep 2
        #make sub_df of subset_rep_2 without the rows of the indices listed in usedQTLlist
        leftoverrep2 = subset_rep_2[-usedQTLvec,]
        leftoverrep2 = cbind(leftoverrep2, leftoverrep2)
        #sometimes there are no leftovers, so skip then
        if (nrow(leftoverrep2) >0 ) {leftoverrep2$number_of_reps = 1 #only found once
        colnames(leftoverrep2) = thecolnames
        merge_frame = rbind(merge_frame, leftoverrep2)}
      }}   #testing per environment here
    #will get a warning when rep1 has a chrm with QTLs but rep 2 doesn't have that chrm at all
  }
}


#Some QTLs might be used twice because there are two peaks within 100kb in its corresponding replicate
#so need to look at repeated peaks from row to row
#make string of each of the two columns
#then can look at how many unique ones there are - compare to how many rows
col1 = paste0(merge_frame$reporter, "_", merge_frame$environment, "_", merge_frame$chrm, "_", merge_frame$rep_1_LOD, merge_frame$rep_1_delta_AF, merge_frame$rep_1_left_Index, merge_frame$rep_1_max_Index, merge_frame$rep_1_right_Index, merge_frame$replicate1)
col2 = paste0(merge_frame$reporter, "_", merge_frame$environment, "_", merge_frame$chrm, "_", merge_frame$rep_2_LOD, merge_frame$rep_2_delta_AF, merge_frame$rep_2_left_Index, merge_frame$rep_2_max_Index, merge_frame$rep_2_right_Index, merge_frame$replicate2)
length(unique(col1)) #same number as col1, so no repeats
length(unique(col2)) #in my data, 3 fewer items, so 3 repeats. need to find
dups = duplicated(col2) #gives T/F vector of which are repeats
dupdf = data.frame(col2, dups) #combines the QTLs with the T/F data
dupdf = dupdf[dupdf$dups== 'TRUE', ] #saves only the rows where there are duplicates

#these are the 3 peaks to remove
#Phe_N-end_SC_4_19.690.1716292576121012510002847003045002
#Thr_N-end_4NQO_14_120.1556805789537634901005146005377002
#Thr_N-end_Low_Glucose_9_59.210.3403757645331932442002657002898002

#make a column in order to find these rows
merge_frame$removeID = paste0(merge_frame$reporter, "_", merge_frame$environment, "_", merge_frame$chrm, "_", merge_frame$rep_2_LOD, merge_frame$rep_2_delta_AF, merge_frame$rep_2_left_Index, merge_frame$rep_2_max_Index, merge_frame$rep_2_right_Index, merge_frame$replicate2)

#look at what they were partnered up with, and figure out which row to keep
mergedups = merge_frame[merge_frame$removeID == "Phe_N-end_SC_4_19.690.1716292576121012510002847003045002",]
mergedups = rbind(mergedups, merge_frame[merge_frame$removeID == "Thr_N-end_4NQO_14_120.1556805789537634901005146005377002",])
mergedups = rbind(mergedups, merge_frame[merge_frame$removeID == "Thr_N-end_Low_Glucose_9_59.210.3403757645331932442002657002898002",])

#keep df just in case I want to refer back to it
write.table(mergedups,
            file = paste0("duplicated_peaks_24.07.12.txt"),
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE,
            sep = '\t')

#can use row name/numbers to remove because they will match merge_frame
#keeping peaks that are more similar (location/LOD) to the peak from the other replicate
#Remove row 903: Phe, SC, 4, 17.19, 0.141600, 312400, 341700, 365400, 1, Phe, SC, 4, 19.69, 0.1716293, 251000, 284700, 304500, 2, 2
#Remove row 115: Thr, 4NQO, 14, 13.85, 0.09842249, 567000, 579100, 585900, 1, Thr, 4NQO, 14, 12.00, 0.15568058, 490100, 514600, 537700, 2, 2
#Remove row 669: Thr, Low_Glucose, 9, 42.11, 0.2537518, 217200, 236000, 262100, 1, Thr, Low_C, 9, 59.21, 0.3403758, 244200, 265700, 289800, 2, 2
#list the row names to remove
rem = c('903', '115', '669')
#make a new df without those rows
merge_frame_filtered = merge_frame[!rownames(merge_frame)%in%rem, ]

#keep final table of QTLs!
write.table(merge_frame_filtered,
            file = paste0(base_dir, proj, "GxE_QTLs_24.07.12.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#average the QTLs found in both replicates
#if the QTL was only found in one replicate, the data from that QTL is just copied in the df, so taking the average keeps the numbers the same

avg_table <- data.frame(reporter = merge_frame_filtered$reporter,
                        environment = merge_frame_filtered$environment,
                        chr = merge_frame_filtered$chrm,
                        LOD = (merge_frame_filtered$rep_1_LOD + merge_frame_filtered$rep_2_LOD) / 2,
                        delta_AF = (merge_frame_filtered$rep_1_delta_AF + merge_frame_filtered$rep_2_delta_AF) / 2,
                        left_Index = (merge_frame_filtered$rep_1_left_Index + merge_frame_filtered$rep_2_left_Index) / 2,
                        max_Index = (merge_frame_filtered$rep_1_max_Index + merge_frame_filtered$rep_2_max_Index) / 2 ,
                        right_Index = (merge_frame_filtered$rep_1_right_Index + merge_frame_filtered$rep_2_right_Index) / 2,
                        replicateA = merge_frame_filtered$replicate1,
                        replicateB = merge_frame_filtered$replicate2,
                        number_of_reps = merge_frame_filtered$number_of_reps)

write.table(avg_table, #416 QTLs found in both replicates. 694 total QTLs
            file = paste0(base_dir, proj, "QTLs_averaged_keptsinglereps24.07.12.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#move to next script to analyze/describe the QTLs and make plots

