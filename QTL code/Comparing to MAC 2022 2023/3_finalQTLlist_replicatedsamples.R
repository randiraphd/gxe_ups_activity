#Code by Randi Avery for 
  #Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity. bioRxiv, 2024.11.21.624644. https://doi.org/10.1101/2024.11.21.624644

#combines peak tables from each sample into one df
#saves all peaks into a table, separated by replicate
#combines QTLs found in both replicates, but keeps QTLs present in only replicate for later analysis
#removes peaks that were used more than once
#averages replicates (for QTLs present in only one replicate, just keeps those values)
#this results in a final list of QTLs to use for further analyses and plots


#used to compare my QTLs to those previously found in Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
# and Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570


base_dir <- "~/myproteasome/mac_compare/"
proj    <- "24.07.18_maccompare/"
proj_dir <- paste0(base_dir, proj, "peaks/")
setwd(proj_dir)


#make list of files with peak info
peaks = list.files(pattern = 'AFDs_4reps_7.28.24.txt$')

#read in files and combine into one df
allpeaks = do.call(rbind, lapply(peaks, function(x) read.table(x, header = TRUE)))

#find any rows that have an NA and send that to a dataframe to look at later. Will skip those in the code below
na_table = allpeaks[rowSums(is.na(allpeaks)) > 0, ] 
allpeaks = na.omit(allpeaks)

#Make names uniform
allpeaks[allpeaks == "rpn4_degron"] = "Rpn4"
allpeaks[allpeaks == "Asn_N-end"] = "Asn"
allpeaks[allpeaks == "Phe_N-end"] = "Phe"
allpeaks[allpeaks == "Thr_N-end"] = "Thr"

#save peaks file
write.table(allpeaks, #should be 193 peaks for comparing me and Mahlon
            file = "all_peaks_useforheatmapMAC_compare24.07.28.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#find QTLs that are present in both replicates, but also keep in only in one replicate for later analyses - just note this
  # same chromosome
  # peaks w/in 100 kb of ea. other
  # same direction of allele frequency difference (positive or negative)

#separate out RRA and MAC peaks
rra = allpeaks[allpeaks$rep_name == 'R1', ]
rra = rbind(rra, allpeaks[allpeaks$rep_name == 'R2', ])

mac = allpeaks[allpeaks$rep_name == 'M1', ]
mac = rbind(mac, allpeaks[allpeaks$rep_name == 'M2', ])


#reorder columns to match the next steps (which I wrote previously)
rra = rra[,c('reporter', 'environment', 'chr', 'LOD', 'AFD', 'left_Index', 'max_Index', 'right_Index', 'rep')]
mac = mac[,c('reporter', 'environment', 'chr', 'LOD', 'AFD', 'left_Index', 'max_Index', 'right_Index', 'rep')]


#make empty data frame to be used for overlaps
#want to keep QTLs that are only in one rep
thecolnames = c('reporter', 'environment', 'chrm', 'rep_1_LOD', 'rep_1_delta_AF', 'rep_1_left_Index', 'rep_1_max_Index', 'rep_1_right_Index', 'replicate1',
                'reporter', 'environment', 'chrm', 'rep_2_LOD', 'rep_2_delta_AF', 'rep_2_left_Index', 'rep_2_max_Index', 'rep_2_right_Index', 'replicate2', 
                'number_of_reps')
merge_frame = data.frame(matrix(ncol = 19, nrow = 0))
colnames(merge_frame) = thecolnames

reporters = unique(rra$reporter) 

#for loop to combine reps. Keep QTLs present in only one rep
for(repor in reporters) {
  #repor= 'Thr' #just for testing one reporter at a time
  sub_reporters = rra[rra$reporter == repor,] #subset to get unique reporters
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


#testing to see if it's the same as when just analyzing my own data: yes it is
# testdf = qtls[qtls$reporter == 'Asn', ]
# testdf = rbind(testdf, qtls[qtls$reporter == 'Phe', ])
# testdf = rbind(testdf, qtls[qtls$reporter == 'Rpn4', ])
# testdf = rbind(testdf, qtls[qtls$reporter == 'Thr', ])
# testdf = testdf[testdf$environment == 'SC', ]

#save df and run again for MAC
rra_merge_frame = merge_frame

thecolnames = c('reporter', 'environment', 'chrm', 'rep_1_LOD', 'rep_1_delta_AF', 'rep_1_left_Index', 'rep_1_max_Index', 'rep_1_right_Index', 'replicate1',
                'reporter', 'environment', 'chrm', 'rep_2_LOD', 'rep_2_delta_AF', 'rep_2_left_Index', 'rep_2_max_Index', 'rep_2_right_Index', 'replicate2', 
                'number_of_reps')
merge_frame = data.frame(matrix(ncol = 19, nrow = 0))
colnames(merge_frame) = thecolnames


#for loop to combine reps. Keep QTLs present in only one rep
for(repor in reporters) {
  #repor= 'Thr' #just for testing one reporter at a time
  sub_reporters = mac[mac$reporter == repor,] #subset to get unique reporters
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

#mac data is merge_frame

#Run this for MAC data
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

#one peak to remove from RRA
#Phe_N-end_SC_4_19.690.1716292576121012510002847003045002


#one peak to remove from MAC
#Phe_SC_7_20.980.2033982585177823973004184004473002

#make a column in order to find these rows
rra_merge_frame$removeID = paste0(rra_merge_frame$reporter, "_", rra_merge_frame$environment, "_", rra_merge_frame$chrm, "_", rra_merge_frame$rep_2_LOD, rra_merge_frame$rep_2_delta_AF, rra_merge_frame$rep_2_left_Index, rra_merge_frame$rep_2_max_Index, rra_merge_frame$rep_2_right_Index, rra_merge_frame$replicate2)
merge_frame$removeID = paste0(merge_frame$reporter, "_", merge_frame$environment, "_", merge_frame$chrm, "_", merge_frame$rep_2_LOD, merge_frame$rep_2_delta_AF, merge_frame$rep_2_left_Index, merge_frame$rep_2_max_Index, merge_frame$rep_2_right_Index, merge_frame$replicate2)

#look at what they were partnered up with, and figure out which row to keep
mergedups_rra = rra_merge_frame[rra_merge_frame$removeID == "Phe_SC_4_19.690.1716292576121012510002847003045002",]
mergedups = merge_frame[merge_frame$removeID == "Phe_SC_7_20.980.2033982585177823973004184004473002",]


#can use row name/numbers to remove because they will match merge_frame
#keeping peaks that are more similar (location/LOD) to the peak from the other replicate
#Remove row 903: Phe, SC, 4, 17.19, 0.141600, 312400, 341700, 365400, 1, Phe, SC, 4, 19.69, 0.1716293, 251000, 284700, 304500, 2, 2
#list the row names to remove
rem_rra = '35'
#make a new df without those rows
rra_merge_frame_filtered = rra_merge_frame[!rownames(rra_merge_frame)%in%rem_rra, ]


#now which to remove from MAC
mergedups$left_diff = abs(mergedups$rep_1_left_Index - mergedups$rep_2_left_Index)
mergedups$right_diff = abs(mergedups$rep_1_right_Index - mergedups$rep_2_right_Index)
mergedups$max_diff = abs(mergedups$rep_1_max_Index - mergedups$rep_2_max_Index)
mergedups$lod_diff = abs(mergedups$rep_1_LOD - mergedups$rep_2_LOD)
#Remove row 69: Phe          SC    7      9.82       0.119454           479800          494000            582200          1
#list the row names to remove
rem_mac = '69'
#make a new df without those rows
mac_merge_frame_filtered = merge_frame[!rownames(merge_frame)%in%rem_mac, ]

#keep final table of QTLs!
write.table(rra_merge_frame_filtered,
            file = paste0(base_dir, proj, "RRA_QTLs_24.07.29.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#keep final table of QTLs!
write.table(mac_merge_frame_filtered,
            file = paste0(base_dir, proj, "MAC_QTLs_24.07.29.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#average the QTLs found in both replicates
#if the QTL was only found in one replicate, the data from that QTL is just copied in the df, so taking the average keeps the numbers the same

avg_table_rra <- data.frame(reporter = rra_merge_frame_filtered$reporter,
                        environment = rra_merge_frame_filtered$environment,
                        chr = rra_merge_frame_filtered$chrm,
                        LOD = (rra_merge_frame_filtered$rep_1_LOD + rra_merge_frame_filtered$rep_2_LOD) / 2,
                        delta_AF = (rra_merge_frame_filtered$rep_1_delta_AF + rra_merge_frame_filtered$rep_2_delta_AF) / 2,
                        left_Index = (rra_merge_frame_filtered$rep_1_left_Index + rra_merge_frame_filtered$rep_2_left_Index) / 2,
                        max_Index = (rra_merge_frame_filtered$rep_1_max_Index + rra_merge_frame_filtered$rep_2_max_Index) / 2 ,
                        right_Index = (rra_merge_frame_filtered$rep_1_right_Index + rra_merge_frame_filtered$rep_2_right_Index) / 2,
                        replicateA = rra_merge_frame_filtered$replicate1,
                        replicateB = rra_merge_frame_filtered$replicate2,
                        number_of_reps = rra_merge_frame_filtered$number_of_reps)

avg_table_mac <- data.frame(reporter = mac_merge_frame_filtered$reporter,
                        environment = mac_merge_frame_filtered$environment,
                        chr = mac_merge_frame_filtered$chrm,
                        LOD = (mac_merge_frame_filtered$rep_1_LOD + mac_merge_frame_filtered$rep_2_LOD) / 2,
                        delta_AF = (mac_merge_frame_filtered$rep_1_delta_AF + mac_merge_frame_filtered$rep_2_delta_AF) / 2,
                        left_Index = (mac_merge_frame_filtered$rep_1_left_Index + mac_merge_frame_filtered$rep_2_left_Index) / 2,
                        max_Index = (mac_merge_frame_filtered$rep_1_max_Index + mac_merge_frame_filtered$rep_2_max_Index) / 2 ,
                        right_Index = (mac_merge_frame_filtered$rep_1_right_Index + mac_merge_frame_filtered$rep_2_right_Index) / 2,
                        replicateA = mac_merge_frame_filtered$replicate1,
                        replicateB = mac_merge_frame_filtered$replicate2,
                        number_of_reps = mac_merge_frame_filtered$number_of_reps)



write.table(avg_table_rra, #39 QTLs found in both replicates. 59 total QTLs
            file = paste0(base_dir, proj, "RRA_QTLs_averaged_keptsinglereps24.07.29.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


write.table(avg_table_mac, #31 QTLs found in both replicates. 59 total QTLs
            file = paste0(base_dir, proj, "MAC_QTLs_averaged_keptsinglereps24.07.29.txt"),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#move to next script to analyze/describe the QTLs and make plots
#use /Users/randia/myproteasome/mac_compare/maccompare_analyzingQTLs_Plots24.07.29.R
