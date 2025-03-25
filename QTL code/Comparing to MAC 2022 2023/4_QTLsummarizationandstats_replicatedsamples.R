#Code by Randi Avery for 
  #Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity. bioRxiv, 2024.11.21.624644. https://doi.org/10.1101/2024.11.21.624644

#Analyzes and describes numbers of QTLs per samples
#Performs overlap analysis to test for replicated QTLs between 
  #Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity. bioRxiv, 2024.11.21.624644. https://doi.org/10.1101/2024.11.21.624644
  #and
  #Collins, M., Avery, R., & Albert, F. W. (2023). Substrate-specific effects of natural genetic variation on proteasome activity. PLOS Genetics, 19(5), e1010734. https://doi.org/10.1371/journal.pgen.1010734
  #Collins, M., Mekonnen, G., & Albert, F. W. (2022). Variation in ubiquitin system genes creates substrate-specific effects on proteasomal protein degradation. eLife, 11, e79570. https://doi.org/10.7554/eLife.79570

#only considers QTLs that were present in Avery et al 2024 (and present or absent in MAC's 2 papers)
#Makes Supplementary Figure 2F

base_dir <- "~/myproteasome/mac_compare/"
proj    <- "24.07.18_maccompare/"
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

#read in QTL data from previous script: ~/myproteasome/paper_final_QTLs_24.07.12.R
qtls_rra = read.table(file = "RRA_QTLs_averaged_keptsinglereps24.07.29.txt", header = T) #Includes QTLs that are present in only 1 replicate
qtls_mac = read.table(file = "MAC_QTLs_averaged_keptsinglereps24.07.29.txt", header = T) #Includes QTLs that are present in only 1 replicate


#only want QTLs in my data that are present in both replicates to plot total QTLs
qtls2_rra = qtls_rra[qtls_rra$number_of_reps == 2,] #39 QTLs found in both replicates



####overlapping analysis####

#add column with reporter and whose data it is
qtls2_rra$rep = paste0(qtls2_rra$reporter, '_RRA')
qtls_mac$rep = paste0(qtls_mac$reporter, '_MAC')

#make df to save the overlapping QTLs

thecolnamesoverlap = c('reporter_1', 'environment_1', 'chrm', 'rep_1_LOD', 'rep_1_delta_AF', 'rep_1_left_Index', 'rep_1_max_Index', 'rep_1_right_Index', 'repA1', 'repB1', 'num_of_reps_1', 'rep1',#will avg values from reps A and B - will tell me if that's from rep 1, 2, or was present in both)
                       'reporter_2', 'environment_2', 'chrm', 'rep_2_LOD', 'rep_2_delta_AF', 'rep_2_left_Index', 'rep_2_max_Index', 'rep_2_right_Index', 'repA2', 'repB2', 'num_of_reps_2', 'rep2')
merge_frame_new = data.frame(matrix(ncol = 24, nrow = 0))
colnames(merge_frame_new) = thecolnamesoverlap

#also make a new df with:
#counts of QTLs in both conditions, total QTLs in both, total overlapping QTLs
cols = c('reporter1', 'environment1', 'total_QTLs1', 'reporter2', 'environment2','total_QTLs2', 'overlapping_QTLs')
QTL_info = data.frame(matrix(ncol = 7, nrow = 0))
colnames(QTL_info) = cols

#Function for overlap - keep QTLs even if only in one replicate, and present/absent
#will only look at consider_sign = TRUE (adapting an older version of code)
find_overlaps <- function(combolist, avgtable, consider_sign = TRUE){#output a list with two dfs is QTLinfo and merge frame, respectively
  for(combo in 1:nrow(combolist)){ #can control what combos to use by making the combolist later
    r1 = combolist[combo,1]
    r2 = combolist[combo,2]
    rep1 = avgtable[avgtable$rep == r1,]
    rep2 = avgtable[avgtable$rep == r2,]
    #record the number of QTLs for both
    #initialize overlapping qtls as 0
    newrow = data.frame(rep1$reporter[1], rep1$environment[1], nrow(rep1), rep2$reporter[1], rep2$environment[1], nrow(rep2), 0)
    colnames(newrow) = cols
    #subset chromosomes
    thiscombo = merge_frame_new[0,]
    #need to go in both directions. and if absent, record that
    for(i in 1:16) { 
      #load the QTL table for
      #replicate 1 as 'rep_1_peaks_table' 
      subset_rep_1 <- rep1[rep1$chr == i, ]
      subset_rep_2 <- rep2[rep2$chr == i, ]
      usedQTLvec = 0 #keep track of which QTLs are used from rep 2
      if (dim(subset_rep_1)[1]==0 & dim(subset_rep_2)[1]==0) {
        print(c(i,r1, r2, ' no QTLs on this chrm'))
      }else {
        
        #check to see if rep 1 is empty
        if (dim(subset_rep_1)[1]==0) {#add QTLs from rep2 to merge_frame
          #make empty df
          df1 = subset_rep_1[NA,]
          df = cbind(subset_rep_2,df1)
          colnames(df) = thecolnamesoverlap
          #change NA to 0 for numb of reps
          df$num_of_reps_2 = 0
          df$repA2 = r1 #want to know what we're comparing for presence/absence 
          df$repB2 = r2
          merge_frame_new = rbind(merge_frame_new, df)
        }else {
          for(j in 1:nrow(subset_rep_1)) {
            ## find peaks w / in 100 kb of ea. other:
            ## make sure to use 'abs' here!
            if(min(abs(subset_rep_1[j,"max_Index"] - subset_rep_2$max_Index)) <= 1e5) 
            { ##min cuz could be more than one peak within 100kb. Abs cuz could be to the left or the right
              ## if we find something w / in 100 kb, grab that row from replicate 2
              subset_2_QTL <- subset_rep_2[which.min(abs(subset_rep_1[j, "max_Index"] - subset_rep_2$max_Index)), ] 
              usedQTLrep2 = which.min(abs(subset_rep_1[j, "max_Index"] - subset_rep_2$max_Index))
              ## make sure they have the same direction of effect:
              #decide if diff direction or not
              if(consider_sign == TRUE){
                print('dont use TRUE') #just edited code previously written
              }
              if(consider_sign == FALSE){
                ## use 'cbind' to make a dataframe row that contains the info from both replicates
                df = cbind(subset_rep_1[j, ], subset_2_QTL)
                colnames(df) = thecolnamesoverlap
                merge_frame_new = rbind(merge_frame_new, df)
                usedQTLvec = c(usedQTLvec, usedQTLrep2)
                #add dataframe for just this combo to count overlaps
                thiscombo = rbind(thiscombo, df)}
              
              #warnings should just be from finding min when there is an empty vector
            } else 
            { #there were no peaks within 100kb, but wanna keep for QTL from rep1 for presence/absence
              df2 = subset_rep_1[NA,]
              df = cbind(subset_rep_1[j, ],df2[1, ]) 
              colnames(df) = thecolnamesoverlap
              #change NA to 0 for numb of reps
              df$num_of_reps_2 = 0
              df$repA2 = r1#want to know what we're comparing for presence/absence 
              df$repB2 = r2
              merge_frame_new = rbind(merge_frame_new, df)}
          } #for(j in 1:nrow(subset_rep_1)) 
          #when we're done with a chromosome, add leftover rows from rep 2
          #make sub_df of subset_rep_2 without the rows of the indices listed in usedQTLlist
          if (length(usedQTLvec) == 1) {
            if (usedQTLvec == 0) { #if usedQTLvec, means nothing from rep 2 was used
              if (nrow(subset_rep_2)>0 ){#had to add this cuz if no QTLs were used, it's gonna be 0, but wanna keep QTLs if there are any in rep2, they just weren't used yet.
                leftoverrep2 = subset_rep_2
                df2 = subset_rep_1[NA,]
                leftoverrep2 = cbind(leftoverrep2, df2[1, ]) #should work correctly even if leftovers have more than 1 row
                colnames(leftoverrep2) = thecolnamesoverlap
                leftoverrep2$num_of_reps_2 = 0
                leftoverrep2$repA2 = r1#want to know what we're comparing for presence/absence 
                leftoverrep2$repB2 = r2
                merge_frame_new = rbind(merge_frame_new, leftoverrep2)}
            }else { #this is what happens if usedQTLvec = 1, but it's not 0
              leftoverrep2 = subset_rep_2[-usedQTLvec,]
              if(dim(leftoverrep2)[1]>0){
                df2 = subset_rep_1[NA,]
                leftoverrep2 = cbind(leftoverrep2, df2[1, ])
                colnames(leftoverrep2) = thecolnamesoverlap
                leftoverrep2$num_of_reps_2 = 0
                leftoverrep2$repA2 = r1#want to know what we're comparing for presence/absence 
                leftoverrep2$repB2 = r2
                merge_frame_new = rbind(merge_frame_new, leftoverrep2)}}
          } else { #this is what happens if usedQTLvec is bigger than 1
            leftoverrep2 = subset_rep_2[-usedQTLvec,]
            if(dim(leftoverrep2)[1]>0){
              df2 = subset_rep_1[NA,]
              leftoverrep2 = cbind(leftoverrep2, df2[1, ])
              colnames(leftoverrep2) = thecolnamesoverlap
              leftoverrep2$num_of_reps_2 = 0
              leftoverrep2$repA2 = r1#want to know what we're comparing for presence/absence 
              leftoverrep2$repB2 = r2
              merge_frame_new = rbind(merge_frame_new, leftoverrep2)}}
        } 
      } #end of block where things happen if both reps have something in them
    }#all chrms
    
    #per combo
    #save info about this combination 
    overlap <- nrow(thiscombo)
    newrow$overlapping_QTLs <- overlap
    QTL_info <- rbind(QTL_info, newrow)
    
  } #all combos
  dflist = list(QTL_info, merge_frame_new)
  return(dflist)
}

#make a table of the comparisons you want - two columns, each row is one comparison

rra_mac_combos = read.table(file = "rra_mac_combos.txt", header = T, sep = '\t') 

#recombine my and Mahlon's qtls
qtls_mac_rra = rbind(qtls2_rra, qtls_mac)

overlapsMR = find_overlaps(rra_mac_combos, qtls_mac_rra, consider_sign = FALSE)
#separate out the dfs from the list
overlapsMR_QTLinfo = overlapsMR[[1]]
overlapsMR_merge = overlapsMR[[2]] #68 overlapping pairs


#now add 'total reps' and 'direction of effect' 
overlapsMR_merge$total_reps = overlapsMR_merge$num_of_reps_1 + overlapsMR_merge$num_of_reps_2
overlapsMR_merge$true_means_same_dir = overlapsMR_merge$rep_1_delta_AF *overlapsMR_merge$rep_2_delta_AF > 0

#change N/As to absent, cuz I wanna look at present/absent
overlapsMR_merge[is.na(overlapsMR_merge)] = "Absent"


#save QTL info table
write.table(overlapsMR_QTLinfo,
            file = "mac_rra_overlappingQTLinfo_24.07.29_savingQTLsin1rep.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#save combined qtls
write.table(overlapsMR_merge,
            file = "mac_rra_overlappingQTLs_forGxE_24.07.29.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#can reload here
#overlaps_not_filtered = read.table(file = "overlappingQTLs_forGxE_24.07.15.txt", header = TRUE, sep = '\t')


#filter out rows that don't have enough QTLs to matter
#remove total QTLs that = 1
overlapsMR_merge_filtered = overlapsMR_merge[overlapsMR_merge$total_reps != 1,] #50 overlapping pairs. Has P/A with 2 reps in one sample and 0 in the other - either in RRA or MAC!
#remove rows that show overlaps, but there is only 1 QTL per sample
#we don't care about rows that say 'true' or 'false' but total qtls = 2, because that means the QTL was only found in 1/2 replicates for each sample
overlapsMR_merge_filtered$remrow = paste0(overlapsMR_merge_filtered$total_reps, overlapsMR_merge_filtered$true_means_same_dir)
overlapsMR_merge_filtered = overlapsMR_merge_filtered[overlapsMR_merge_filtered$remrow != '2TRUE',  ] #still 50
overlapsMR_merge_filtered = overlapsMR_merge_filtered[overlapsMR_merge_filtered$remrow != '2FALSE',  ] #still 50


#2Absent means present in 2 reps, absent in 2 reps (P/A GxE),
#3FALSE means one direction in two reps, and the opposite direction in one rep (sign change GxE),
#3TRUE means one direction in two reps, and the same direction in one rep (no GxE),
#4FALSE means one direction in two reps, and the opposite direction in two reps (sign change GxE),
#4TRUE means one direction in two reps, and the same direction in two reps (no GxE)


print(c('Total overlapping pairs:', nrow(overlapsMR_merge_filtered)), quote = F)

#save overlapping qtls
write.table(overlapsMR_merge_filtered,
            file = "mac_rra_overlappingQTLs_forGxE_24.07.29_PAeitherdirection.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')



#DON'T KEEP QTLS THAT WERE IN 2 REPS OF MAHLONS AND NONE OF MINE - because we are using MY QTLs as the focal list
#because in the combolist I always listed my sample first, the only case this will be true is P/A and MAC is in the rep1 column
overlapsMR_merge_filtered2 = overlapsMR_merge_filtered[overlapsMR_merge_filtered$rep1 != 'Asn_MAC',] #47
overlapsMR_merge_filtered2 = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$rep1 != 'Phe_MAC',] #44
overlapsMR_merge_filtered2 = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$rep1 != 'Rpn4_MAC',] #40
overlapsMR_merge_filtered2 = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$rep1 != 'Thr_MAC',] #39 #should be the number of focal QTLs we start out with

#save overlapping qtls
write.table(overlapsMR_merge_filtered2,
            file = "mac_rra_overlappingQTLs_forGxE_24.07.29_RRAfocal.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

###reload point###
overlapsMR_merge_filtered2 = read.table(file = "mac_rra_overlappingQTLs_forGxE_24.07.29_RRAfocal.txt", header = TRUE, sep = '\t')



#GxE stat plots

####Want to know per REPORTER:####
# #"2Absent" "3FALSE" "4FALSE" = GxE 
# #"3TRUE" "4TRUE" = overlapping, but not GxE


#need to separate out 'Absent' since the info is set up differently
#1/3 DF#
mr_notabsent = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$true_means_same_dir != 'Absent',] #30
mr_absent = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$true_means_same_dir == 'Absent',] #9


#for sign change or no GxE ('notabsent') (will result in all 0's for P/A)

#initiate GxE function- use for all gxestats functions - can change 'Reporter' to 'Environment'
statscols = c('Reporter', 'Total Overlapping Pairs', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Presence/Absense GxE (2Absent)','Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
gxestats = data.frame(matrix(ncol = 7, nrow = 0))
colnames(gxestats) = statscols

#write a function
gxestats_fun = function(overlapdf, subset_by){
  for (repo in unique(overlapdf[,subset_by])){ 
    subrows = overlapdf[overlapdf[,subset_by] == repo,] #all overlapping pairs of that reporter or environment
    myrow = gxestats[0,]
    myrow = data.frame(repo, nrow(subrows), nrow(subrows[subrows$remrow == '4FALSE',]), nrow(subrows[subrows$remrow == '3FALSE',]), nrow(subrows[subrows$remrow == '2Absent',]),  nrow(subrows[subrows$remrow == '4TRUE',]), nrow(subrows[subrows$remrow == '3TRUE',]))
    colnames(myrow) = statscols
    gxestats = rbind(gxestats, myrow)
  } 
  return(gxestats)
}

#run function for not absent gxe
mr_notabsentstat = gxestats_fun(mr_notabsent, "rep1")

#now do P/A
#ONLY EVER PRESENT IN MY DATA AND ABSENT IN MAHLON'S DATA (not doing the other way around)

#run function
mr_PA = gxestats_fun(mr_absent, 'rep1')


#need to combine the 3 dfs, so reorganize with the data I want to keep, and rename the columns
keepsRep = c('Reporter', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
mr_notabsentstat1 = mr_notabsentstat[keepsRep]
mr_PA$'Present in RRA; Absent in MAC' = mr_PA$`Presence/Absense GxE (2Absent)`
keeps2mr = c('Reporter', 'Present in RRA; Absent in MAC')
mr_notabsentstat2 = mr_PA[keeps2mr]

mr_stats_done = merge(mr_notabsentstat1,mr_notabsentstat2, by = 'Reporter')
mr_stats_done$`Total Overlapping Pairs` = mr_stats_done$`Sign change GxE (4FALSE)` + mr_stats_done$`Sign change GxE (3FALSE)` + mr_stats_done$`Overlapping but NO GxE (4TRUE)` + mr_stats_done$`Overlapping but NO GxE (3TRUE)` + mr_stats_done$`Present in RRA; Absent in MAC`
mr_stats_done$`Missing %` = (mr_stats_done$`Present in RRA; Absent in MAC` / mr_stats_done$`Total Overlapping Pairs`) *100
mr_stats_done$`% Replicated` = 100-mr_stats_done$`Missing %`
#Total overlapping pairs = total QTLs we started with

write.table(mr_stats_done,
            file = "RRA_MACoverlapping_stats24.07.29.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')



#Want to see how the QTLs that didn't replicate compare those that did in terms of LOD and dAF
mr_notabsent = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$true_means_same_dir != 'Absent',] #30
mr_absent = overlapsMR_merge_filtered2[overlapsMR_merge_filtered2$true_means_same_dir == 'Absent',] #9


mr_notabsent$key1 = paste0(mr_notabsent$rep1, mr_notabsent$chrm, ceiling(mr_notabsent$rep_1_LOD), round(mr_notabsent$rep_1_delta_AF, digits = 3), mr_notabsent$rep_1_left_Index, mr_notabsent$rep_1_max_Index, mr_notabsent$rep_1_right_Index)
#have to force as.numeric for when it says 'absent'
mr_notabsent$key2 = paste0(mr_notabsent$rep2, mr_notabsent$chrm, ceiling(as.numeric(mr_notabsent$rep_2_LOD)), round(as.numeric(mr_notabsent$rep_2_delta_AF), digits = 3), mr_notabsent$rep_2_left_Index, mr_notabsent$rep_2_max_Index, mr_notabsent$rep_2_right_Index)

#keys for QTLs that were replicated
keys = c(mr_notabsent$key1, mr_notabsent$key2)

#make key in master list of QTLs
qtls_mac_rra$key = paste0(qtls_mac_rra$rep, qtls_mac_rra$chr, ceiling(qtls_mac_rra$LOD), round(qtls_mac_rra$delta_AF, digits = 3), qtls_mac_rra$left_Index, qtls_mac_rra$max_Index, qtls_mac_rra$right_Index)

replicatedQTLs = qtls_mac_rra[0,]
replicatedQTLs = qtls_mac_rra[qtls_mac_rra$key %in% keys,]
replicatedQTLs$replicated = 'TRUE'

#now do again for those QTLs that did not replicate
mr_absent$key1 = paste0(mr_absent$rep1, mr_absent$chrm, ceiling(mr_absent$rep_1_LOD), round(mr_absent$rep_1_delta_AF, digits = 3), mr_absent$rep_1_left_Index, mr_absent$rep_1_max_Index, mr_absent$rep_1_right_Index)

#keys for QTLs that were replicated
keys_absent = mr_absent$key1

notreplicatedQTLs = qtls_mac_rra[0,]
notreplicatedQTLs = qtls_mac_rra[qtls_mac_rra$key %in% keys_absent,]
notreplicatedQTLs$replicated = 'FALSE'

annotatedQTLs = rbind(replicatedQTLs, notreplicatedQTLs)
#but prob only want to look at stats for just my QTLs, cuz I only have my QTLs for P/A
rra_annotatedQTLs = annotatedQTLs[annotatedQTLs$rep == 'Asn_RRA',]
rra_annotatedQTLs = rbind(rra_annotatedQTLs, annotatedQTLs[annotatedQTLs$rep == 'Phe_RRA',])
rra_annotatedQTLs = rbind(rra_annotatedQTLs, annotatedQTLs[annotatedQTLs$rep == 'Thr_RRA',])
rra_annotatedQTLs = rbind(rra_annotatedQTLs, annotatedQTLs[annotatedQTLs$rep == 'Rpn4_RRA',])
#total should be 39 QTLs - the number of QTLs found in these 4 samples in both replicates

write.table(rra_annotatedQTLs,
            file = "RRA_QTLs_replicated24.07.29.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#reload
rra_annotatedQTLs = read.table(file = "RRA_QTLs_replicated24.07.29.txt", header = TRUE, sep = '\t')

rra_annotatedQTLs[rra_annotatedQTLs =='FALSE'] = 'Did not replicate'
rra_annotatedQTLs[rra_annotatedQTLs =='Did not replicate'] = 'Did not\nreplicate'
rra_annotatedQTLs[rra_annotatedQTLs =='TRUE'] = 'Replicated'

#Supplementary Figure 2F
#boxplot
ggplot(rra_annotatedQTLs, aes(x=replicated,  y=LOD)) +
  geom_boxplot(outlier.shape = NA, position = 'identity', alpha = 0.4, width=.5) +
  labs(x = NULL,  y="LOD")+
  geom_dotplot( binaxis='y',  dotsize= 10 , binwidth = 1, alpha = 0.3,  stackdir='center', position = position_jitter(w = 0.05, h = 0))+ #binwidth 0.01 gets rid of course bins. dotsize and binwidth are connected
  theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 21),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"))
#saved as 3x3.25

t.test(LOD ~ replicated, data = rra_annotatedQTLs)

#boxplot
ggplot(rra_annotatedQTLs, aes(x=replicated,  y=abs(delta_AF))) +
  geom_boxplot(outlier.shape = NA, position = 'identity', alpha = 0.4, width=.5) +
  labs(x = NULL,  y="|∆AF|")+
  geom_dotplot(binaxis='y',  dotsize= 30, binwidth = 0.001, alpha = 0.3,  stackdir='center', position = position_jitter(w = 0.09, h = 0))+ #binwidth 0.01 gets rid of course bins. dotsize and binwidth are connected
  theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 21),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"))
#saved as 3x3.25

t.test(abs(delta_AF) ~ replicated, data = rra_annotatedQTLs)
