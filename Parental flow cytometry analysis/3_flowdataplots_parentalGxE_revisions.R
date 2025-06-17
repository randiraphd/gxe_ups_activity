#Written by Randi Avery in 2024 for Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity.

## load all required packages for stats
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

setwd("~/myproteasome/reporter_characterization/environments/combined")

#reload data table produced from 2_upload_2.4.2025_flowanalysis.R
gated_medians = read.delim(file = 'filtered_gated_medians2.7.25.txt', header = T, sep = "\t")



#Do environments have extra autofluorescence? Check fluorescence of negative controls (strains with no TFTs)
#basically just seeing autofluroescence of the media
negctrls = gated_medians[gated_medians$environment != "SC", ] 
negctrls = negctrls[negctrls$environment != "YNB", ]
negctrls = negctrls[negctrls$environment != "LowN", ]
negctrls = negctrls[negctrls$environment != "LowG", ]
negctrls = negctrls[negctrls$environment != "AZC", ]
negctrls = negctrls[negctrls$environment != "4NQO", ]
negctrls = negctrls[negctrls$environment != "BTZ", ]
negctrls = negctrls[negctrls$environment != "LiAc", ]

#now plot all neg ctrls
#want x-axis to be reporter and environment
ggplot(negctrls, aes(x=grp2,  y=log_GFP, fill = strain)) +
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "Reporter and Environment", y = "log_GFP") +
  geom_dotplot(binaxis='y', dotsize= 1 , alpha = 0.4, binwidth = 0.01, stackdir='center')+ 
  #geom_jitter(color= 'pink', size=2, alpha=0.9) +
  #geom_point(position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10),) +
  ggtitle("Autofluorescence - Negative Controls from Flow in Environments") +
  scale_fill_manual(values = c(col_by, col_rm))

#ynb looked higher. Plotted YNB compared to SC for all the points I have
#If YNB points are always higher than SC points, it would be due to autofluorescence and not changes in TFT
ynb = gated_medians[gated_medians$environment == "YNB", ]
sc = gated_medians[gated_medians$environment == "SC", ]
ynb = rbind(ynb, sc)
ggplot(ynb, aes(x=grp2,  y=log_GFP, fill = reporter)) +
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "Environment", y = "log_GFP") +
  geom_dotplot(binaxis='y', dotsize= 1 , alpha = 0.4, binwidth = 0.01, stackdir='center')+ #seems like dotsize and binwidth are connected
  #geom_jitter(color= 'pink', size=2, alpha=0.9) +
  #geom_point(position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle("Is YNB always higher than SC?")
#doing a t-test to see if SC and YNB are diff with GFP
t.test(log_GFP ~ environment, data = ynb)
#no, it is not. two sample t-test p-value = 0.6911


####----

#Did each environment change degradation in the direction we expected? Compared TFT ratio from each environment to SC.

#now making plots to see if diff environments had the effect we expected
btz = gated_medians[gated_medians$environment == "BTZ", ] #for the other environments, just replaced variables in this section by hand
sc = gated_medians[gated_medians$environment == "SC", ]
btz = rbind(btz, sc)
#diff x-axis
btz$grp3 = paste0(btz$reporter, '_', btz$environment)

ggplot(btz, aes(x=grp3,  y=TFT_ratio, fill = strain)) +
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity', alpha = 0.7) + 
  labs( x = "Reporter and Environment", y = "-log2(RFP/GFP)") +
  geom_dotplot(binaxis='y', dotsize= 5 , alpha = 0.4, binwidth = 0.01, stackdir='center')+ #dotsize and binwidth are related
  #geom_jitter(color= 'pink', size=2, alpha=0.9) +
  #geom_point(position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  ggtitle("Did the environments do what we thought?") +
  scale_fill_manual(values = c(col_by, col_rm))
#for the other environments, just replaced variables in the above section by hand

####----

#####-----##### LINEAR MODEL

gated_medians = read.delim(file = 'filtered_gated_medians2.7.25.txt', header = T, sep = "\t")

#model = lmer(df$TFT_loess ~ df$strain_final * df$environment + (1 | df$chr_rep))

#For loop for doing all comparisons
###
##before loop:
#make function to extract values from model: https://stackoverflow.com/questions/29923986/how-to-reformat-an-r-data-frame-with-multiple-rows-into-one-row 
# Turn the above steps into a function
colToRow <- function(x) {
  melted <- melt(x[, -1], id.vars = NULL)
  row <- t(melted[, 2])
  colnames(row) <- paste0(rownames(x),"_",  melted[, 1])
  row
}


gated_medians[gated_medians == "LowG"] = "Low G"
gated_medians[gated_medians == "LowN"] = "Low N"

#make lists in order that we want the plots to be in
#environments = c("LowG", "LowN", "YNB", "LiAc", "4NQO", "AZC", "BTZ")
reporters = c("Asn", "Phe", "Thr", "UFD", "4xUb", "Rpn4")
environments = rev(c("BTZ",   "AZC",   "4NQO",  "LiAc",  "YNB",   "Low N", "Low G"))


#4.10.25 add means, new rank change, if rank change is the same, and SD to the table


#make df to save model results
dfcolnames = c("intercept_SE", "strain_SE", "env_SE", "interaction_SE", "intercept_df", "strain_df", "env_df", "interaction_df", "intercept_Tval", "strain_Tval", "env_Tval", "interaction_Tval", "intercept_Pval", "strain_Pval", "env_Pval", "interaction_Pval")
results_table_loess = data.frame(matrix(ncol = 16, nrow = 0))
colnames(results_table_loess) = dfcolnames

#df to save final results with info on rank change
ranknames = c("intercept_Pval", "strain_Pval", "env_Pval", "interaction_Pval", 
              "reporter", "environment1", "environment2", "Int_Sig", "medscBY", "medscRM", "medEnvBY", "medEnvRM", 
              "RMhigherinSC", "RMhigherinEnv", "rankchange", "strain_PvalSC", "sig", "strain_PvalEnv", "sig", "rankchange_sig", 
              "env_PvalBY", "env_PvalRM", 'BYsig80', 'RMsig80',
              
              "meanscBY", "stdevSCBY", "meanscRM", "stdevSCRM", "meanEnvBY", "stdevEnvBY", "meanEnvRM", "stdevEnvRM",#revision 4.29.25
              "RMhigherinSCmean", "RMhigherinEnvmean", "rankchangemean", "RC_med_mean_same?")#revision 4.29.25
rankchange = data.frame(matrix(ncol = 36, nrow = 0))
colnames(rankchange) = ranknames


##NEED TO MAKE NEW REPLICATE SO THAT IT'S NOT TREATED LIKE A NUMBER
#BY1 in SC is the same strain as BY1 in LiAc, etc
gated_medians$chr_rep = paste0(gated_medians$strain, gated_medians$replicate)

#pvalue cut off
p40 =0.05/40
p80 = 0.05/80


#4/10/25 commented out plot - uncomment to make plots again

#Fig. 2DEF
pdf(file = "GxE_TFTtimeloess_2.7.25_maxlabels.pdf", # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5)

#Which SC to use - for 4xUb and UFD I ran a second set of SC samples on the second plate. I will use those samples for the environments on those plates.

#saved individual plots for Fig 2 as 4x5 pdfs
#r = 'Phe'
#e2 = '4NQO'

#r = 'Asn'
#e2 = 'AZC'

#r = 'Thr'
#e2 = 'Low G'


for(r in reporters){ # subset the reporter by the environments
  
  print(r) #used to check
  sub_rep = gated_medians[gated_medians$reporter == r,] #this is the data for one reporter
  #for each comparison within that reporter
  for(env in environments){
    tryCatch({ #since 4xUb and UFD for Low N were thrown out, needed to add this to skip those tests
      e1 = "SC"
      e2 = env
      rep2 = sub_rep[sub_rep$environment == e2,]
      rep1 = sub_rep[sub_rep$environment == "SC",]
      if(rep2$plate[1] == 'plate02' & r == 'UFD' | rep2$plate[1] == 'plate02' & r == '4xUb'){rep1 = rep1[rep1$plate == 'plate02',]} else #when SC was ran on the second plate, and the other environment was on the second plate, use the data from the second plate
        {rep1 = rep1[rep1$plate == 'plate01',]}
      print(nrow(rep1))
      
      #make df to see sig diff due to environment on each strain. doing it here, so the right SC is used
      df = rbind(rep1, rep2)
      by = df[df$strain == 'BY',]
      rm = df[df$strain == 'RM',]
      
      
      #get medians to find rank changes
      medscBY = median(rep1[rep1$strain== 'BY', 'TFT_loess'])
      medscRM = median(rep1[rep1$strain== 'RM', 'TFT_loess'])
      medEnvBY = median(rep2[rep2$strain== 'BY', 'TFT_loess'])
      medEnvRM = median(rep2[rep2$strain== 'RM', 'TFT_loess'])
      rmhigherinSC = medscRM > medscBY
      rmhigherinEnv = medEnvRM > medEnvBY
      
      #t-tests between strains to see if the rank change is significant
      scpval = t.test(TFT_loess ~ strain, data = rep1) #SC
      envpval = t.test(TFT_loess ~ strain, data = rep2) #other env
      
      #t-tests between environments to see if the environment had a significant effect
      bypval = t.test(TFT_loess ~ environment, data = by) 
      rmpval = t.test(TFT_loess ~ environment, data = rm) 
           
      #df = rbind(rep1, rep2)
      #print(df$plate)
      #print(df$environment)
      #run LINEAR MODEL
      model = lmer(df$TFT_loess ~ df$strain_final * df$environment + (1 | df$chr_rep))
      #extract values as needed - use function
      tmodel = as.data.frame(coef(summary(model)))
      dfresult = as.data.frame(colToRow(tmodel)) #use function - save as data frame
      colnames(dfresult) = c('intercept_SE', 'strain_SE', 'env_SE', 'interaction_SE', 'intercept_df', 'strain_df', 'env_df', 'interaction_df', 'intercept_Tval', 'strain_Tval', 'env_Tval', 'interaction_Tval', 'intercept_Pval', 'strain_Pval', 'env_Pval', 'interaction_Pval')
      #keep only needed columns
      keeps = c('intercept_Pval', 'strain_Pval', 'env_Pval', 'interaction_Pval')
      dfresult = dfresult[keeps]
      #add columns with info on what was compared
      dfresult$reporter = r
      dfresult$environment1 = e1
      dfresult$environment2 = e2
      dfresult$Int_Sig = dfresult$interaction_Pval < p40
      dfresult$medscBY = medscBY
      dfresult$medscRM = medscRM
      dfresult$medEnvBY = medEnvBY
      dfresult$medEnvRM = medEnvRM
      dfresult$rmhigherinSC = rmhigherinSC
      dfresult$rmhigherinEnv = rmhigherinEnv
      dfresult$rankchange = rmhigherinEnv != rmhigherinSC
      dfresult$strain_PvalSC = scpval$p.value
      dfresult$pvalSCsig = dfresult$strain_PvalSC < p40
      dfresult$strain_PvalEnv = envpval$p.value
      dfresult$pvalEnvsig = dfresult$strain_PvalEnv < p40
      dfresult$rankchangesig = dfresult$rankchange == TRUE & dfresult$pvalSCsig == TRUE & dfresult$pvalEnvsig == TRUE
      dfresult$env_PvalBY = bypval$p.value
      dfresult$env_PvalRM = rmpval$p.value
      dfresult$BYsig = dfresult$env_PvalBY < p80
      dfresult$RMsig = dfresult$env_PvalRM < p80
      
      #following block is from revisions
      dfresult$meansscBY = mean(rep1[rep1$strain== 'BY', 'TFT_loess'])
      dfresult$stdevSCBY = sd(rep1[rep1$strain== 'BY', 'TFT_loess'])
      dfresult$meanscRM = mean(rep1[rep1$strain== 'RM', 'TFT_loess'])
      dfresult$stdevSCRM = sd(rep1[rep1$strain== 'RM', 'TFT_loess'])
      dfresult$meanEnvBY = mean(rep2[rep2$strain== 'BY', 'TFT_loess'])
      dfresult$stdevEnvBY = sd(rep2[rep2$strain== 'BY', 'TFT_loess'])
      dfresult$meanEnvRM = mean(rep2[rep2$strain== 'RM', 'TFT_loess'])
      dfresult$stdevEnvRM = sd(rep2[rep2$strain== 'RM', 'TFT_loess'])
      RMhigherinSCmean = mean(rep1[rep1$strain== 'RM', 'TFT_loess']) > mean(rep1[rep1$strain== 'BY', 'TFT_loess'])
      dfresult$RMhigherinSCmean = RMhigherinSCmean
      RMhigherinEnvmean = mean(rep2[rep2$strain== 'RM', 'TFT_loess']) > mean(rep2[rep2$strain== 'BY', 'TFT_loess'])
      dfresult$RMhigherinEnvmean = RMhigherinEnvmean
      medianrankchange = rmhigherinEnv != rmhigherinSC
      meanrankchange = RMhigherinEnvmean != RMhigherinSCmean
      dfresult$rankchangemean = meanrankchange
      dfresult$RC_med_mean_same = medianrankchange == meanrankchange #only two cases where this didn't match, but they were not significant rank changes

      
      #combine this row to the rows already made
      rankchange = rbind(rankchange, dfresult)

      # #now make GxE plot
      # p = 'TFT_loess'
      # plotdf = df
      # yvalue = plotdf[[p]]
      # mytitle = paste0(unique(plotdf$reporter), " Reporter  | GxE: ", dfresult$interaction_Pval < 0.00125,
      #                  '\n', "Interaction p-value: ", formatC(dfresult$interaction_Pval, format = 'e', digits = 2), ' | ', "env_p=", formatC(dfresult$env_Pval, format = 'e', digits = 2), ' | ', "strain_p=", formatC(dfresult$strain_Pval, format = 'e', digits = 2),
      #                  '\nSignificant rank change: ', dfresult$rankchangesig, " | T-test p = SC: ", formatC(scpval$p.value, format = 'e', digits = 2), ", ", e2, ": ", formatC(envpval$p.value, format = 'e', digits = 2))
      # myplot = ggplot(plotdf, aes(x=environment,  y=yvalue, fill = strain)) +
      #   geom_boxplot(outlier.size = 0.3, position = 'identity', alpha = 0.4, width=.4) +
      #   labs(x = expression(atop(bold("Environment"))), y=expression(atop(bold("UPS Activity"),atop(italic("-log2(RFP/GFP)")))))+
      #   geom_dotplot(binaxis='y', dotsize= 20 , alpha = 0.8, binwidth = 0.01, stackdir='center', position = position_jitter(w = 0.05, h = 0))+ #binwidth 0.01 gets rid of course bins. dotsize and binwidth are connected
      #   scale_x_discrete(limits = c("SC", e2)) + #e2 should maybe be environment2
      #   theme(axis.text.x = element_text(size = 18),
      #         axis.text.y = element_text(size = 18),
      #         axis.title = element_text(size = 15),
      #         plot.title = element_text(hjust = 0.5),
      #         panel.background = element_blank(),
      #         panel.grid.major.x = element_blank(),
      #         panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
      #         axis.line = element_line(colour = "grey")) +
      #         ggtitle(mytitle) +
      #   scale_fill_manual(values = c(col_by, col_rm)) +
      #   #guides(fill=guide_legend(title="Strain"))+
      #   ylim(-2,4)+ #gets all plots to fit on one scale
      #   stat_summary(
      #     fun = median,
      #     geom = 'line',
      #     aes(group = strain, color = strain)) +
      #   scale_color_manual(values = c(col_by, col_rm))+
      #   labs(fill='Strain') +
      #   labs(color='Strain') 
      # print(myplot)

    #-end of for loop-
    }, error=function(e){cat("ERROR :",conditionMessage(e), env, "\n")})
  }}

dev.off()

##Counting up how many higher/lower due to environment
rankchange$BYdecreased = rankchange$medscBY > rankchange$medEnvBY
rankchange$RMdecreased = rankchange$medscRM > rankchange$medEnvRM
rankchange$BYdecSig = paste0(rankchange$BYdecreased, "_", rankchange$BYsig)
rankchange$RMdecSig = paste0(rankchange$RMdecreased, "_", rankchange$RMsig)
#TRUE_TRUE = environment caused a significant decrease in UPS activity
#FALSE_TRUE = environment caused a significant increase in UPS activity
table(rankchange$BYdecSig)
table(rankchange$RMdecSig)

#counting changes per reporter
effectofenvBY = rankchange %>%
  group_by(reporter, BYdecSig) %>%
  summarize(Freq=n())

effectofenvBY$strain = 'BY'
colnames(effectofenvBY) <- c("reporter", "Sig_decrease", "Freq",     "strain")

effectofenvRM = rankchange %>%
  group_by(reporter, RMdecSig) %>%
  summarize(Freq=n())

effectofenvRM$strain = 'RM'
colnames(effectofenvRM) <- c("reporter", "Sig_decrease", "Freq",     "strain")

effectofenv = rbind(effectofenvBY, effectofenvRM)
#Result: Sig Decreases by Reporter - 4xUb 6, Asn 8, Phe 5, Rpn4 4, Thr 8, UFD 11


#output df into a table file
# result table 
write.table(rankchange,
            file = "results_lmer_40tests_04.29.25_TFTtime_loess_revisions.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#reload data in new sessions
###!!! "rankchange" is the new name of "results_table_loess" - use "results_table_loess" for the rest of the code, that was written before the updated "rankchange"
results_table_loess = read.delim(file = 'results_lmer_40tests_2.7.25_TFTtime_loess.txt', header = T, sep = "\t")


#want to see if sig diff of presence of GxE between categories
results_table_loess$ubsystem = ifelse(results_table_loess$rep == '4xUb' |results_table_loess$rep == 'Rpn4', 'Independent', 'Dependent')

boxplot(interaction_Pval ~ ubsystem, results_table_loess)
wilcox.test(interaction_Pval ~ ubsystem, results_table_loess, exact = FALSE)

results_table_loess$starve = ifelse(results_table_loess$environment2 == 'Low N' |results_table_loess$environment2 == 'Low G'|results_table_loess$environment2 == 'YNB' , 'Starvation', 'Not Starvation')
boxplot(interaction_Pval ~ starve, results_table_loess)
wilcox.test(interaction_Pval ~ starve, results_table_loess, exact = FALSE)


#make a df with the number of tests passing/not passing significance for plotting
intdf = results_table_loess %>% group_by(reporter, Int_Sig) %>% tally()


#Supplementary Fig. 1B
#order based on number of significant tests (made plot first, then made order based on first/unordered plot)
level_order2 <- c('Asn', 'UFD', 'Rpn4', 'Thr', '4xUb', 'Phe')
ggplot(intdf, aes(fill=Int_Sig, y=n, x=factor(reporter, level = level_order2))) +
  geom_bar(position="fill", stat="identity")+ #position = stack for non-percentage. fill for percentage
  scale_fill_manual(values = c(mymaroon, mygold))+
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.title=element_text(size=15)) + 
  labs( x = "Reporter", y ='Proportion of significant tests for GxE') +
  guides(fill=guide_legend(title="Bonferroni \nsignificance"))
#saved plot in 4.25x5.5 pdf


#diff version of plot above
#only want to plot straight numbers, not proportion
intdf_true = intdf[intdf$Int_Sig == TRUE,]
ggplot(intdf_true, aes(y=n, x=reorder(reporter, n))) +
  geom_bar(position="stack", stat="identity")+ #position = stack for non-percentage. fill for percentage
  #scale_fill_manual(values = c(mymaroon, mygold))+
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey")) + 
  labs( x = "Reporter", y ='Number of significant tests for GxE') 
#UFD and 4xUb have n=6


#now the same, but grouped by environment
#make a df with the number of tests passing/not passing significance
intdf_env = results_table_loess %>% group_by(environment2, Int_Sig) %>% tally()

#Supplementary Fig. 1C
level_order3 <- c('BTZ', '4NQO','Low N', 'Low G', 'AZC', 'LiAc', 'YNB')
ggplot(intdf_env, aes(fill=Int_Sig, y=n, x=factor(environment2, level = level_order3))) +
  geom_bar(position="fill", stat="identity")+ #position = stack for non-percentage. fill for percentage
  scale_fill_manual(values = c(mymaroon, mygold))+
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.title=element_text(size=15)) + 
  labs( x = "Environment", y ='Proportion of significant tests for GxE') +
  guides(fill=guide_legend(title="Bonferroni \nsignificance"))
#saved plot in 4.25x6 pdf


#diff version of plot above
#only want to plot straight numbers, not proportion
intdf_env_true = intdf_env[intdf_env$Int_Sig == TRUE,]
ggplot(intdf_env_true, aes(y=n, x=reorder(environment2, n))) +
  geom_bar(position="stack", stat="identity")+ #position = stack for non-percentage. fill for percentage
  #scale_fill_manual(values = c(mymaroon, mygold))+
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey")) + 
  labs( x = "Environment", y ='Number of significant tests for GxE') 



#Now creating a plot that summarizes all parental GxE (flow data)
#values of replicates are very tight, so using the average

#make df to save averages for each reporter and environment
dfcolnames_means = c("reporter", "environment", "strain", "plate", "avg_TFT_ratio", "n")
means_table = data.frame(matrix(ncol = 6, nrow = 0))
colnames(means_table) = dfcolnames_means

#want it separated based on strain and plate as well
gated_medians$grp4 = paste0(gated_medians$grp2, "_", gated_medians$strain, "_", gated_medians$plate)

#calculate the average loess time corrected -log2(RFP/GFP) for each set of replicates
for(rep in unique(gated_medians$grp4)){ # subset the reporter by the environments
  #print(rep) #used to check
  sub_rep = gated_medians[gated_medians$grp4 == rep,] #this is the data for one reporter in one environment in one plate
  avgTFT = mean(sub_rep$TFT_loess) #calculate mean degradation activity (time corrected)
  newrow = c(sub_rep$reporter[1], sub_rep$environment[1], sub_rep$strain[1], sub_rep$plate[1], avgTFT, nrow(sub_rep))
  means_table = rbind(means_table, newrow)
  colnames(means_table) = dfcolnames_means
}

#revisions adding error
dfcolnames_means = c("reporter", "environment", "strain", "plate", "avg_TFT_ratio", "n", "SD")
means_table_error = data.frame(matrix(ncol = 7, nrow = 0))
colnames(means_table_error) = dfcolnames_means
for(rep in unique(gated_medians$grp4)){ # subset the reporter by the environments
  #print(rep) #used to check
  sub_rep = gated_medians[gated_medians$grp4 == rep,] #this is the data for one reporter in one environment in one plate
  avgTFT = mean(sub_rep$TFT_loess) #calculate mean degradation activity (time corrected)
  sd = sd(sub_rep$TFT_loess)
  newrow = c(sub_rep$reporter[1], sub_rep$environment[1], sub_rep$strain[1], sub_rep$plate[1], avgTFT, nrow(sub_rep), sd)
  means_table_error = rbind(means_table_error, newrow)
  colnames(means_table_error) = dfcolnames_means
}



#also do for medians
dfcolnames_medians = c("reporter", "environment", "strain", "plate", "median_TFT_ratio", "n")
medians_table = data.frame(matrix(ncol = 6, nrow = 0))
colnames(medians_table) = dfcolnames_medians

for(rep in unique(gated_medians$grp4)){ # subset the reporter by the environments
  #print(rep) #used to check
  sub_rep = gated_medians[gated_medians$grp4 == rep,] #this is the data for one reporter in one environment in one plate
  medianTFT = median(sub_rep$TFT_loess) #calculate median degradation activity (time corrected)
  newrow = c(sub_rep$reporter[1], sub_rep$environment[1], sub_rep$strain[1], sub_rep$plate[1], medianTFT, nrow(sub_rep))
  medians_table = rbind(medians_table, newrow)
  colnames(medians_table) = dfcolnames_medians
}




#don't need negative controls in plot
means_table = means_table[means_table$environment != 'negctrl4NQO',]
means_table = means_table[means_table$environment != 'negctrlLiAc',]
means_table = means_table[means_table$environment != 'negctrlLowG',]
means_table = means_table[means_table$environment != 'negctrlLowN',]
means_table = means_table[means_table$environment != 'negctrlSC',]
means_table = means_table[means_table$environment != 'negctrlYNB',]


#don't need negative controls in plot
means_table_error = means_table_error[means_table_error$environment != 'negctrl4NQO',]
means_table_error = means_table_error[means_table_error$environment != 'negctrlLiAc',]
means_table_error = means_table_error[means_table_error$environment != 'negctrlLowG',]
means_table_error = means_table_error[means_table_error$environment != 'negctrlLowN',]
means_table_error = means_table_error[means_table_error$environment != 'negctrlSC',]
means_table_error = means_table_error[means_table_error$environment != 'negctrlYNB',]

#do again for medians table
medians_table = medians_table[medians_table$environment != 'negctrl4NQO',]
medians_table = medians_table[medians_table$environment != 'negctrlLiAc',]
medians_table = medians_table[medians_table$environment != 'negctrlLowG',]
medians_table = medians_table[medians_table$environment != 'negctrlLowN',]
medians_table = medians_table[medians_table$environment != 'negctrlSC',]
medians_table = medians_table[medians_table$environment != 'negctrlYNB',]

##create df to match corresponding SC's with each environmental sample##

#want to group based on reporter and strain
means_table$grp = paste0(means_table$reporter, '_', means_table$strain)
#in order to match which SC to subtract, use grp5 (only for 4xUb and UFD)
means_table$grp5 = paste0(means_table$grp, '_', means_table$plate)

means_table_error$grp = paste0(means_table_error$reporter, '_', means_table_error$strain)
#in order to match which SC to subtract, use grp5 (only for 4xUb and UFD)
means_table_error$grp5 = paste0(means_table_error$grp, '_', means_table_error$plate)

#i don't think i needed to add anything to means_table_error from here!

#save means and medians tables
write.table(means_table,
            file = "parentalGxE_meansof8reps.2.7.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

means_table = read.table(file = "parentalGxE_meansof8reps.2.7.25.txt", header = T, sep = '\t') 

write.table(means_table_error,
            file = "parentalGxE_meansof8reps_error4.14.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

means_table_error = read.table(file = "parentalGxE_meansof8reps_error4.14.25.txt", header = T, sep = '\t') 




write.table(medians_table,
            file = "parentalGxE_mediansof8reps.2.7.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

medians_table = read.table(file = "parentalGxE_mediansof8reps.2.7.25.txt", header = T, sep = '\t') 
#want to group based on reporter and strain
medians_table$grp = paste0(medians_table$reporter, '_', medians_table$strain)
#in order to match which SC to subtract, use grp5 (only for 4xUb and UFD)
medians_table$grp5 = paste0(medians_table$grp, '_', medians_table$plate)

#separate SC - will look at difference between SC and other environments for TFT_loess
scmeans = means_table[means_table$environment == 'SC',]
othermeans = means_table[means_table$environment != 'SC',]

scmeans_error = means_table_error[means_table_error$environment == 'SC',]
othermeans_error = means_table_error[means_table_error$environment != 'SC',]


scmedians = medians_table[medians_table$environment == 'SC',]
othermedians = medians_table[medians_table$environment != 'SC',]

#since we are grouping which SC to use in a different way for 4xUb and UFD, we need to separate those
ufd_scmeans = scmeans[scmeans$reporter == 'UFD',]
fourx_scmeans = scmeans[scmeans$reporter == '4xUb',]
ufd4x_scmeans = rbind(ufd_scmeans,fourx_scmeans) #contains SC values for 4xUb and UFD
ufd_othermeans = othermeans[othermeans$reporter == 'UFD',]
fourx_othermeans = othermeans[othermeans$reporter == '4xUb',]
ufd4x_othermeans = rbind(ufd_othermeans,fourx_othermeans) #contains environment values for 4xUb and UFD

reporters4_scmeans = scmeans[scmeans$reporter != 'UFD',]
reporters4_scmeans = reporters4_scmeans[reporters4_scmeans$reporter != '4xUb',] #contains SC values for the other 4 reporters
reporters4_othermeans = othermeans[othermeans$reporter != 'UFD',]
reporters4_othermeans = reporters4_othermeans[reporters4_othermeans$reporter != '4xUb',] #contains environment values for the other 4 reporters


#medians chunk
ufd_scmedians = scmedians[scmedians$reporter == 'UFD',]
fourx_scmedians = scmedians[scmedians$reporter == '4xUb',]
ufd4x_scmedians = rbind(ufd_scmedians,fourx_scmedians) #contains SC values for 4xUb and UFD
ufd_othermedians = othermedians[othermedians$reporter == 'UFD',]
fourx_othermedians = othermedians[othermedians$reporter == '4xUb',]
ufd4x_othermedians = rbind(ufd_othermedians,fourx_othermedians) #contains environment values for 4xUb and UFD

reporters4_scmedians = scmedians[scmedians$reporter != 'UFD',]
reporters4_scmedians = reporters4_scmedians[reporters4_scmedians$reporter != '4xUb',] #contains SC values for the other 4 reporters
reporters4_othermedians = othermedians[othermedians$reporter != 'UFD',]
reporters4_othermedians = reporters4_othermedians[reporters4_othermedians$reporter != '4xUb',] #contains environment values for the other 4 reporters


#now combine df's so that i can extract the correct value
#need to combine SC with other environments. And do separately for 4xUb/UFD and the 4 other reporters (since 4xUb/UFD have an extra SC from plate02 I can use)
#library(dplyr)
#for 4xUb and UFD match up using grp5
ufd4x = left_join(ufd4x_othermeans, ufd4x_scmeans, by=c("grp5", "reporter", "strain", "grp"))
#for other reporters, match up using grp
reporters4 = left_join(reporters4_othermeans, reporters4_scmeans, by=c("grp", "reporter", "strain"))
#now need to edit ufd4x so I can combine it with reporters4
names(ufd4x)[names(ufd4x) == 'grp5'] <- 'grp5.x'
ufd4x$grp5.y = ufd4x$grp5.x #just need to repeat the column so i can combine with reporters4

normalizedmeans = rbind(ufd4x, reporters4)

#now subtract TFT_loess of SC from TFT_loess of the given environment
normalizedmeans$delta_TFT_ratio = as.numeric(normalizedmeans$avg_TFT_ratio.x) - as.numeric(normalizedmeans$avg_TFT_ratio.y) #y is the SC value, based on the way the dfs were built


#error chunk - didn't change variable names here! just rerun! #added for revisions
#since we are grouping which SC to use in a different way for 4xUb and UFD, we need to separate those
ufd_scmeans = scmeans_error[scmeans_error$reporter == 'UFD',]
fourx_scmeans = scmeans_error[scmeans_error$reporter == '4xUb',]
ufd4x_scmeans = rbind(ufd_scmeans,fourx_scmeans) #contains SC values for 4xUb and UFD
ufd_othermeans = othermeans_error[othermeans_error$reporter == 'UFD',]
fourx_othermeans = othermeans_error[othermeans_error$reporter == '4xUb',]
ufd4x_othermeans = rbind(ufd_othermeans,fourx_othermeans) #contains environment values for 4xUb and UFD

reporters4_scmeans = scmeans_error[scmeans_error$reporter != 'UFD',]
reporters4_scmeans = reporters4_scmeans[reporters4_scmeans$reporter != '4xUb',] #contains SC values for the other 4 reporters
reporters4_othermeans = othermeans_error[othermeans_error$reporter != 'UFD',]
reporters4_othermeans = reporters4_othermeans[reporters4_othermeans$reporter != '4xUb',] #contains environment values for the other 4 reporters
#now combine df's so that i can extract the correct value
#need to combine SC with other environments. And do separately for 4xUb/UFD and the 4 other reporters (since 4xUb/UFD have an extra SC from plate02 I can use)
#library(dplyr)
#for 4xUb and UFD match up using grp5
ufd4x = left_join(ufd4x_othermeans, ufd4x_scmeans, by=c("grp5", "reporter", "strain", "grp"))
#for other reporters, match up using grp
reporters4 = left_join(reporters4_othermeans, reporters4_scmeans, by=c("grp", "reporter", "strain"))
#now need to edit ufd4x so I can combine it with reporters4
names(ufd4x)[names(ufd4x) == 'grp5'] <- 'grp5.x'
ufd4x$grp5.y = ufd4x$grp5.x #just need to repeat the column so i can combine with reporters4

normalizedmeans_error = rbind(ufd4x, reporters4)

#now subtract TFT_loess of SC from TFT_loess of the given environment
normalizedmeans$delta_TFT_ratio = as.numeric(normalizedmeans$avg_TFT_ratio.x) - as.numeric(normalizedmeans$avg_TFT_ratio.y) #y is the SC value, based on the way the dfs were built

normalizedmeans_error$delta_TFT_ratio = as.numeric(normalizedmeans_error$avg_TFT_ratio.x) - as.numeric(normalizedmeans_error$avg_TFT_ratio.y) #y is the SC value, based on the way the dfs were built


#medians chunk #added for revisions
#now combine df's so that i can extract the correct value
#need to combine SC with other environments. And do separately for 4xUb/UFD and the 4 other reporters (since 4xUb/UFD have an extra SC from plate02 I can use)
#library(dplyr)
#for 4xUb and UFD match up using grp5
mediansufd4x = left_join(ufd4x_othermedians, ufd4x_scmedians, by=c("grp5", "reporter", "strain", "grp"))
#for other reporters, match up using grp
mediansreporters4 = left_join(reporters4_othermedians, reporters4_scmedians, by=c("grp", "reporter", "strain"))
#now need to edit ufd4x so I can combine it with reporters4
names(mediansufd4x)[names(mediansufd4x) == 'grp5'] <- 'grp5.x'
mediansufd4x$grp5.y = mediansufd4x$grp5.x #just need to repeat the column so i can combine with reporters4

normalizedmedians = rbind(mediansufd4x, mediansreporters4)

#now subtract TFT_loess of SC from TFT_loess of the given environment
normalizedmedians$delta_TFT_ratio = as.numeric(normalizedmedians$median_TFT_ratio.x) - as.numeric(normalizedmedians$median_TFT_ratio.y) #y is the SC value, based on the way the dfs were built



## -- now plot: summary of GxE tests -- ##

#want to highlight points that are significant for GxE
#use results_table, and create a column to relate to the normalizedmeans
results_table_loess$grp3 = paste0(results_table_loess$reporter, '_', results_table_loess$environment2)
normalizedmeans$grp3 = paste0(normalizedmeans$reporter, '_', normalizedmeans$environment.x)
#can combine dfs, and results_table will just repeat rows accordingly
gxesummary = left_join(results_table_loess, normalizedmeans, by = c('grp3', 'reporter'))
#keep only columns I want
keeps = c('reporter', 'interaction_Pval', 'Int_Sig', 'environment.x', 'strain', 'delta_TFT_ratio')
gxesummary = gxesummary[keeps]

#revisions Medians chunk
#want to highlight points that are significant for GxE
#use results_table, and create a column to relate to the normalizedmeans
results_table_loess$grp3 = paste0(results_table_loess$reporter, '_', results_table_loess$environment2)
normalizedmedians$grp3 = paste0(normalizedmedians$reporter, '_', normalizedmedians$environment.x)
#can combine dfs, and results_table will just repeat rows accordingly
gxesummary_medians = left_join(results_table_loess, normalizedmedians, by = c('grp3', 'reporter'))
#keep only columns I want
keeps = c('reporter', 'interaction_Pval', 'Int_Sig', 'environment.x', 'strain', 'delta_TFT_ratio')
gxesummary_medians = gxesummary_medians[keeps]


#save df
#output df into a table file
write.table(gxesummary,
            file = "gxesummary_2.7.25_TFTtime_loess.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#reload data in new sessions
gxesummary = read.delim(file = 'gxesummary_2.7.25_TFTtime_loess.txt', header = T, sep = "\t")

write.table(gxesummary_medians,
            file = "gxesummary_medians_3.26.25_TFTtime_loess.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#reload data in new sessions
gxesummary_medians = read.delim(file = 'gxesummary_medians_3.26.25_TFTtime_loess.txt', header = T, sep = "\t")



#need to get environment p-values for dots from results_table_loess
#BYsig and RMsig are the columns i want
newcols = results_table_loess[,c("reporter", "environment2", "BYsig", "RMsig")]

newcolslong = pivot_longer(newcols, cols = c('BYsig', 'RMsig'),
                                 names_to = 'Strain',
                                 values_to = 'Environment Effect')

newcolslong[newcolslong =='BYsig'] = 'BY'
newcolslong[newcolslong =='RMsig'] = 'RM'
#match column names to gxesummary
colnames(newcolslong) <- c('reporter', 'environment.x', 'strain', 'Environment Effect')


#add these new columns to gxesummary
gxesummary = merge(gxesummary, newcolslong, by = c('reporter', 'environment.x', 'strain'))

gxesummary_medians = merge(gxesummary_medians, newcolslong, by = c('reporter', 'environment.x', 'strain'))


#Final parental GxE summary plot
#use to assign alpha
#first alpha is based on interaction
alpha <- ifelse(gxesummary$Int_Sig, 1, 0.65) #i think these numbers are irrelevant - changes based on scale_alpha later
alpha_env <-ifelse(gxesummary$`Environment Effect`, 1, 0.65)

alpha <- ifelse(gxesummary_medians$Int_Sig, 1, 0.65) #i think these numbers are irrelevant - changes based on scale_alpha later
alpha_env <-ifelse(gxesummary_medians$`Environment Effect`, 1, 0.65)

#change order in the way the reporters are displayed
reporter_factor = c('Asn', 'Phe', 'Thr', 'UFD', '4xUb', 'Rpn4')

#environment colors, will need to make a new one that includes SC for those plots
env_colors = c("4NQO" = '#BC71C9', "AZC" = '#6878D2', 'BTZ' = '#00B6EB',
               'LiAc' = '#00C094', 'Low G' = '#C49A02', 'Low N' = '#FFB512', 'YNB' = '#F85735')


#Fig. 2A
ggplot(gxesummary, aes(y=delta_TFT_ratio, x=strain)) +
  geom_hline(yintercept = 0, color = "grey", size = 1)+ #make 0 thick
  stat_summary( #connects medians
    fun = median,
    geom = 'line', size = 2, lineend = 'round',
    aes(group = environment.x, color = environment.x, alpha = alpha))+
  scale_color_manual(values = env_colors)+
  geom_dotplot(aes(fill = environment.x, alpha = alpha_env), color = NA,  binaxis='y', stackdir='center', dotsize=1.1)+ #color = environment.x for circle outline
  scale_fill_manual(values = env_colors) +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Environment Effect on UPS Activity"),atop(italic("Difference in average -log2(RFP/GFP) from baseline")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(reporter, level = reporter_factor)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  scale_alpha(range = c(.3,1)) +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Environment', color = 'Environment') #changes legend title. Have to do both fill and color for dots and lines

#Fig. 2A revisions. Did not use
ggplot(gxesummary_medians, aes(y=delta_TFT_ratio, x=strain)) +
  geom_hline(yintercept = 0, color = "grey", size = 1)+ #make 0 thick
  stat_summary( #connects medians
    fun = median,
    geom = 'line', size = 2, lineend = 'round',
    aes(group = environment.x, color = environment.x, alpha = alpha))+
  scale_color_manual(values = env_colors)+
  geom_dotplot(aes(fill = environment.x, alpha = alpha_env), color = NA,  binaxis='y', stackdir='center', dotsize=1.1)+ #color = environment.x for circle outline
  scale_fill_manual(values = env_colors) +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Environment Effect on UPS Activity"),atop(italic("Difference in average -log2(RFP/GFP) from baseline")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(reporter, level = reporter_factor)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  scale_alpha(range = c(.3,1)) +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Environment', color = 'Environment') #changes legend title. Have to do both fill and color for dots and lines




level_order_e = rev(c("BTZ",   "AZC",   "4NQO",  "LiAc",  "YNB",   "Low N", "Low G"))

#have to add reporter order to the DF
gxesummary$reporter <-factor(gxesummary$reporter, levels = reporter_factor)

gxesummary_medians$reporter <-factor(gxesummary_medians$reporter, levels = reporter_factor)


#Fig. 2B
#switch reporters/environments
ggplot(gxesummary, aes(y=delta_TFT_ratio, x=strain)) +
  geom_hline(yintercept = 0, color = "grey", size = 1)+ #make 0 thick
  stat_summary( #connects medians
    fun = median,
    geom = 'line', size = 2, lineend = 'round',
    aes(group = reporter, color = reporter, alpha = alpha))+
  scale_color_brewer(palette = "Dark2") +
  geom_dotplot(aes(fill = reporter, alpha = alpha_env), color = NA,  binaxis='y', stackdir='center', dotsize=1.1)+ #color = environment.x for circle outline
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Environment Effect on UPS Activity"),atop(italic("Difference in average -log2(RFP/GFP) from baseline")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(environment.x, level = level_order_e)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  scale_alpha(range = c(.3,1)) +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Reporter', color = 'Reporter') #changes legend title. Have to do both fill and color for dots and lines


###More numbers to report on###

#number of negative and positive changes compared to SC
table(sign(gxesummary$delta_TFT_ratio))

#want to know if RM decreased more often than BY
by = gxesummary[gxesummary$strain == 'BY', ]
rm = gxesummary[gxesummary$strain == 'RM', ]

table(sign(by$delta_TFT_ratio))
table(sign(rm$delta_TFT_ratio))


#plotting just straight up avg TFT_loess, so you can see how BY and RM compare

#plotting straight up avg TFT_loess, but switching reporter and environment
sc_env_colors = c("SC" = 'black', "4NQO" = '#BC71C9', "AZC" = '#6878D2', 'BTZ' = '#00B6EB',
                  'LiAc' = '#00C094', 'Low G' = '#C49A02', 'Low N' = '#FFB512', 'YNB' = '#F85735')

#change order of environments
env_factor = c("SC", "4NQO", "AZC", 'BTZ', 'LiAc' , 'Low G' , 'Low N' , 'YNB')
means_table$environment <-factor(means_table$environment, levels = env_factor)
#have to add a column to make size of SC line bigger
means_table$thelines = FALSE
means_table[means_table$environment == 'SC', 'thelines'] <- TRUE
#average the SC's for 4xUb and UFD
means_table$grp6 = paste0(means_table$environment, '_', means_table$plate)


medians_table$environment <-factor(medians_table$environment, levels = env_factor)
medians_table$grp6 = paste0(medians_table$environment, '_', medians_table$plate)

alpha2 <- ifelse(means_table$thelines, 1, 0.65)
reporter_factor = c('Asn', 'Phe', 'Thr', 'UFD', '4xUb', 'Rpn4')

#Supplementary Fig. 1A
ggplot(means_table, aes(y = as.numeric(avg_TFT_ratio), x=strain)) +
  stat_summary( #connects medians
    fun = median,
    geom = 'line', lineend = 'round', size = 1.5, 
    aes(group = grp6, color = environment, alpha = alpha2))+ #have to do group six cuz there's two sets of SC for 4xUb and UFD
  scale_color_manual(values = sc_env_colors)+
  geom_dotplot(aes(fill = environment), color = NA,  binaxis='y', stackdir='center', dotsize=1.1, alpha = alpha2)+ #color = environment.x for circle outline, , alpha = alpha_env
  scale_fill_manual(values = sc_env_colors) +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Mean UPS Activity"),atop(italic("Mean -log2(RFP/GFP)")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(reporter, level = reporter_factor)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  scale_alpha(range = c(.6,1), guide = 'none') +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Environment', color = 'Environment') #changes legend title. Have to do both fill and color for dots and lines
#saved as 4x6 inches pdf

#revisions - plot with SD
means_table_error$environment <-factor(means_table_error$environment, levels = env_factor)
means_table_error$grp6 = paste0(means_table_error$environment, '_', means_table_error$plate)

#need value for error bars to be added/subtracted from actual point
means_table_error$sdmin = as.numeric(means_table_error$avg_TFT_ratio) - as.numeric(means_table_error$SD)
means_table_error$sdmax = as.numeric(means_table_error$avg_TFT_ratio) + as.numeric(means_table_error$SD)

#revision - New Fig. 1A. Also changing alpha to be the same, not based on values
#USE THIS ONE!
ggplot(means_table_error, aes(y = as.numeric(avg_TFT_ratio), x=strain)) +
  stat_summary( #connects medians. Changing to mean
    fun = mean,
    geom = 'line', lineend = 'round', size = 1.5, 
    aes(group = grp6, color = environment, alpha = 0.4))+ #have to do group six cuz there's two sets of SC for 4xUb and UFD
  scale_color_manual(values = sc_env_colors)+
  geom_dotplot(aes(fill = environment), color = NA,  binaxis='y', stackdir='center', dotsize=1.1, alpha = 0.6)+ #color = environment.x for circle outline, , alpha = alpha_env
  geom_errorbar(aes(ymin = sdmin, ymax = sdmax, color=factor(environment)), width = 0.25,size=0.25)+ #revision
  scale_fill_manual(values = sc_env_colors) +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Mean UPS Activity"),atop(italic("Mean -log2(RFP/GFP)")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(reporter, level = reporter_factor)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  scale_alpha(range = c(.6,1), guide = 'none') +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Environment', color = 'Environment') #changes legend title. Have to do both fill and color for dots and lines
#saved as 4x7 inches pdf into new folder

#revisions - reordered by environment
#New plot New Fig 2B

level_order_new = c("Low G", "Low N", "YNB", "LiAc", "4NQO", "AZC", "BTZ", "SC")
means_table_error$reporter <-factor(means_table_error$reporter, levels = reporter_factor)


means_table_error$lines = paste0(means_table_error$reporter, "_", means_table_error$plate) #so the 2 diff SC's for 4xUb and UFD get separate lines

ggplot(means_table_error, aes(y = as.numeric(avg_TFT_ratio), x=strain)) +
  stat_summary( #connects medians. Changing to mean. Doesn't matter - based on one point.
    fun = mean,
    geom = 'line', lineend = 'round', size = 1.5, 
    aes(group = lines, color = reporter))+ #, alpha = 0.9
  #scale_color_manual(values = sc_env_colors)+
  scale_color_brewer(palette = "Dark2") +
  geom_dotplot(aes(fill = reporter), color = NA,  binaxis='y', stackdir='center', dotsize=1.1, alpha = 0.7)+ #color = environment.x for circle outline, , alpha = alpha_env
  geom_errorbar(aes(ymin = sdmin, ymax = sdmax, color=factor(reporter)), width = 0.25,size=0.25)+ #revision
  #scale_fill_manual(values = sc_env_colors) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Mean UPS Activity"),atop(italic("Mean -log2(RFP/GFP)")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(environment, level = level_order_new)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  #scale_alpha(range = c(.6,1), guide = 'none') +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Reporter', color = 'Reporter') #changes legend title. Have to do both fill and color for dots and lines
#saved as 4x7.5


# ggplot(gxesummary, aes(y=delta_TFT_ratio, x=strain)) +
#   geom_hline(yintercept = 0, color = "grey", size = 1)+ #make 0 thick
#   stat_summary( #connects medians
#     fun = median,
#     geom = 'line', size = 2, lineend = 'round',
#     aes(group = reporter, color = reporter, alpha = alpha))+
#   scale_color_brewer(palette = "Dark2") +
#   geom_dotplot(aes(fill = reporter, alpha = alpha_env), color = NA,  binaxis='y', stackdir='center', dotsize=1.1)+ #color = environment.x for circle outline
#   scale_fill_brewer(palette = "Dark2") +
#   theme(axis.text=element_text(size=13), panel.background = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
#         axis.line = element_line(colour = "grey"),
#         strip.text.x = element_text(size = 13)) + #for facet labels
#   labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Environment Effect on UPS Activity"),atop(italic("Difference in average -log2(RFP/GFP) from baseline")))))+
#   theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
#   facet_grid(. ~ factor(environment.x, level = level_order_e)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
#   scale_alpha(range = c(.3,1)) +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
#   labs(fill = 'Reporter', color = 'Reporter') #changes legend title. Have to do both fill and color for dots and lines




#revisions - plot with MEDIANS
ggplot(medians_table, aes(y = as.numeric(median_TFT_ratio), x=strain)) +
  stat_summary( #connects medians
    fun = median,
    geom = 'line', lineend = 'round', size = 1.5, 
    aes(group = grp6, color = environment, alpha = alpha2))+ #have to do group six cuz there's two sets of SC for 4xUb and UFD
  scale_color_manual(values = sc_env_colors)+
  geom_dotplot(aes(fill = environment), color = NA,  binaxis='y', stackdir='center', dotsize=1.1, alpha = alpha2)+ #color = environment.x for circle outline, , alpha = alpha_env
  scale_fill_manual(values = sc_env_colors) +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        strip.text.x = element_text(size = 13)) + #for facet labels
  labs(x = expression(atop(bold("Strain"))), y=expression(atop(bold("Median UPS Activity"),atop(italic("Median -log2(RFP/GFP)")))))+
  theme(axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))+
  facet_grid(. ~ factor(reporter, level = reporter_factor)) + #separates into 6 subplots. can add switch = 'x' to put reporter labels on the bottom
  scale_alpha(range = c(.6,1), guide = 'none') +#makes the insignificant tests pale, and ", guide = 'none' " removes alpha legend
  labs(fill = 'Environment', color = 'Environment') #changes legend title. Have to do both fill and color for dots and lines
#saved as 4x6 inches pdf



##Want to run a t-test for every line in this plot to see if strains are sig diff, and which is higher##
#each t-test will be made up of 8 data points for each of the 2 groups
#take out negative controls
ttests_gated_medians = gated_medians[gated_medians$environment != 'negctrl4NQO',]
ttests_gated_medians = ttests_gated_medians[ttests_gated_medians$environment != 'negctrlLiAc',]
ttests_gated_medians = ttests_gated_medians[ttests_gated_medians$environment != 'negctrlLowG',]
ttests_gated_medians = ttests_gated_medians[ttests_gated_medians$environment != 'negctrlLowN',]
ttests_gated_medians = ttests_gated_medians[ttests_gated_medians$environment != 'negctrlSC',]
ttests_gated_medians = ttests_gated_medians[ttests_gated_medians$environment != 'negctrlYNB',]
#4xUb and UFD have two sets of SC, so will just keep both together - comparing SC to SC anyway
ttestnames = c("reporter", "environment", "mean of x (BY)", "mean of y (RM)", "p.value")
ttests = data.frame(matrix(ncol = 5, nrow = 0))
colnames(ttests) = ttestnames

for(env in unique(ttests_gated_medians$environment)){
  thisenv = ttests_gated_medians[ttests_gated_medians$environment == env,]
  for(r in unique(thisenv$reporter)){
    thisreporter = thisenv[thisenv$reporter == r,]
    bydf = thisreporter[thisreporter$strain == 'BY', ]
    rmdf = thisreporter[thisreporter$strain == 'RM', ]
    myttest = t.test(bydf$TFT_loess, rmdf$TFT_loess)
    myrow = c(r, env, myttest$estimate[1], myttest$estimate[2], myttest$p.value)
    ttests = rbind(ttests,myrow)
    colnames(ttests) = ttestnames
  }
}

ttestpval = 0.05/46 #cutoff value for significance
#make sure columns are numeric
ttests$`mean of x (BY)` = as.numeric(ttests$`mean of x (BY)`)
ttests$`mean of y (RM)` = as.numeric(ttests$`mean of y (RM)`)
ttests$p.value = as.numeric(ttests$p.value)

#add column if test was significant or not
ttests$Pval_sig = ttests$p.value < ttestpval
#add column if BY or RM was higher
ttests$rmhigher = ttests$`mean of y (RM)` > ttests$`mean of x (BY)`
table(ttests$Pval_sig)
#add column for significant and which strain is higher
ttests$results = paste0(ttests$Pval_sig, '_RMhigher_', ttests$rmhigher)
table(ttests$results)

#save table
write.table(ttests,
            file = "parentalstrains_ttests_DiffBtwnStrains_2.7.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#Didn't include these final 2 plots in final version of the paper
## RM looks like it is usually higher than BY - plot:
means_table$grp2 = paste0(means_table$reporter, '_', means_table$environment)
bymeans = means_table[means_table$strain == 'BY', ]
rmmeans = means_table[means_table$strain == 'RM', ]
strainmeans = merge(bymeans,rmmeans, by = c('grp2','reporter','environment', 'plate'))

#calculate which strain is higher by reporter. Must include as.numeric, otherwise it doesn't handle negatives correctly
strainmeans$rmhigher = as.numeric(strainmeans$avg_TFT_ratio.x) < as.numeric(strainmeans$avg_TFT_ratio.y) #x = BY, y = RM. TRUE will mean RM is higher
table(strainmeans$rmhigher) #times RM was higher in total
rmhigher = strainmeans %>% group_by(rmhigher, reporter) %>% tally()
rmhigher = as.data.frame(rmhigher)
#if there is a row that needs to be added because there is a 0:
rmhigher[nrow(rmhigher) + 1,] = c("FALSE","Thr", 0)
#now only keep the TRUE values and plot those
#rmhigher = rmhigher[rmhigher$rmhigher == 'TRUE', ]
#change false and true to by and rm
rmhigher[rmhigher == "FALSE"] = "BY"
rmhigher[rmhigher == "TRUE"] = "RM"

#can tell it to order based on number, not alphabetical. Do based on times RM is higher
level_order4 <- c('UFD','4xUb', 'Phe', 'Asn', 'Rpn4','Thr') 
ggplot(rmhigher, aes(fill = rmhigher, y=as.numeric(n), x=factor(reporter, level = level_order4))) + 
  geom_bar(position=position_dodge(), stat="identity")+ #position = stack for non-percentage. fill for percentage. postion dodge
  #ylim(0,8)+
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey")) + 
  scale_fill_manual(values = c(col_by, col_rm))+
  labs( x = "Reporter", y ='Times Each Strain Had Higher Degradation Activity') +
  ggtitle('Samples each strain had higher degradation activity')+
  labs(fill = 'Strain')


#now tally by environment
rmhigherenv = strainmeans %>% group_by(rmhigher, environment) %>% tally()
rmhigherenv = as.data.frame(rmhigherenv)
rmhigherenv[nrow(rmhigherenv) + 1,] = c("FALSE","Low G", 0)
#change false and true to by and rm
rmhigherenv[rmhigherenv == "FALSE"] = "BY"
rmhigherenv[rmhigherenv == "TRUE"] = "RM"
#rmhigherenv = rmhigherenv[rmhigherenv$rmhigher == 'TRUE', ]

level_order5 <- c('AZC', 'LiAc', 'LowN', 'YNB', '4NQO', 'SC', 'BTZ', 'LowG')
#can tell it to order based on number, not alphabetical
ggplot(rmhigherenv, aes(fill = rmhigher, y=as.numeric(n), x=factor(environment, level = level_order5))) +
  geom_bar(position=position_dodge(), stat="identity")+ #position = stack for non-percentage. fill for percentage
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey")) + 
  scale_fill_manual(values = c(col_by, col_rm))+
  labs( x = "Environment", y ='Times Each Strain Had Higher Degradation Activity')+
  ggtitle('Samples each strain had higher degradation activity')+
  labs(fill = 'Strain')
