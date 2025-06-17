#Code by Randi Avery for 
  #Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity. bioRxiv, 2024.11.21.624644. https://doi.org/10.1101/2024.11.21.624644

#Analyzes and describes numbers of QTLs per samples
#Performs overlap analysis to test for GxE in the QTLs
  #tests GxE per reporter, per sample, and per environment
#Analyzes 'uniqueness' of QTLs (if it was found in at least one other environment)
#Looks to see if GxE QTLs are 'close' to hotspots or known loci
#Makes lots of plots

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

#read in QTL data from previous script: ~/myproteasome/paper_final_QTLs_24.07.12.R
qtls = read.table(file = "QTLs_averaged_keptsinglereps24.07.12.txt", header = T) #694 QTLs. Includes QTLs that are present in only 1 replicate

#at this point, change name of rpn4_redo to rpn4 so I can find overlaps, and how I want it in the plot
qtls[qtls == "rpn4_degron_redo"] = "Rpn4"
qtls[qtls == "rpn4_degron"] = "Rpn4"
#change names of reporters and environments to how I want them in the plot:
qtls[qtls == "4x_Ub"] = "4xUb"
qtls[qtls == "Asn_N-end"] = "Asn"
qtls[qtls == "Phe_N-end"] = "Phe"
qtls[qtls == "Thr_N-end"] = "Thr"
qtls[qtls == "Low_Glucose"] = "Low G"
qtls[qtls == "Low_Nitrogen"] = "Low N"
qtls[qtls == "Bortezomib"] = "BTZ"


####plots summarizing total QTLs####
#only want QTLs that are present in both replicates to plot total QTLs
qtls2 = qtls[qtls$number_of_reps == 2,] #416 QTLs found in both replicates

#add categories for reporters and environments
qtls2cats = qtls2
qtls2cats$starve = ifelse(qtls2cats$environment == 'Low N' |qtls2cats$environment == 'Low G'|qtls2cats$environment == 'YNB' , 'Starvation', 'Not Starvation')
qtls2cats$ubsystem = ifelse(qtls2cats$reporter == '4xUb' |qtls2cats$reporter == 'Rpn4', 'Independent', 'Dependent')
qtls2cats$cats = paste0(qtls2cats$starve, "_", qtls2cats$ubsystem)
#qtls2cats$chem = ifelse(qtls2cats$environment == '4NQO' |qtls2cats$environment == 'AZC'|qtls2cats$environment == 'BTZ' , 'Chemical', 'Starvation')
#qtls2cats[qtls2cats$environment == 'SC', 'chem'] <- 'Other'
#qtls2cats[qtls2cats$environment == 'LiAc', 'chem'] <- 'Other'

#do a wilcox for starvation vs Not Starvation
#chemstarvetotqtls = chemstarve %>% count(reporter, environment)
#chemstarvetotqtls$chem = ifelse(chemstarvetotqtls$environment == '4NQO' |chemstarvetotqtls$environment == 'AZC'|chemstarvetotqtls$environment == 'BTZ' , 'Chemical', 'Starvation')
#boxplot(n ~ chem, chemstarvetotqtls)
#wilcox.test(n ~ chem, data = chemstarvetotqtls, exact = FALSE)
#sum(chemstarvetotqtls[which(chemstarvetotqtls$chem == 'Starvation'), "n"]) 

write.table(qtls2cats,
            file = "416QTLs_10.15.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

qtls2cats = read.table(file = "416QTLs_10.15.24.txt", header = T, sep = '\t')


#revisions - adding additional info to QTL table - supplementary Table 4
#want to add whether the QTL was found previously in either MAC paper

#import table previously made on comparing data to our two previous papers
replicated = read.table(file = "/Users/randia/myproteasome/mac_compare/24.07.18_maccompare/RRA_QTLs_replicated24.07.29.txt", header = T, sep = '\t')

#make keys to find matching QTLs
qtls2cats$key = paste0(qtls2cats$reporter, qtls2cats$chr, ceiling(qtls2cats$LOD), round(qtls2cats$delta_AF, digits = 3), qtls2cats$left_Index, qtls2cats$max_Index, qtls2cats$right_Index)
replicated$key = paste0(replicated$reporter, replicated$chr, ceiling(replicated$LOD), round(replicated$delta_AF, digits = 3), replicated$left_Index, replicated$max_Index, replicated$right_Index)

#keep only info that I want to merge
newcol = replicated[c("key", "replicated")]

#merge DFs with info about being replicated
newQTLtable = merge(x = qtls2cats, y = newcol, by = "key", all = TRUE)

length(na.omit(newQTLtable$replicated)) #check. should = 39
newQTLtable$key = NULL #remove "key" column


#another revision is to include the genes within a CI and the total number of BY-RM variants in that CI
variants = read.table(file = "RM-BYvariants.txt", header = T, sep = '\t')

variants = separate(data = variants, col = Uploaded_variation, into = c('chrm', 'info'), sep = ':')
variants$chrm[variants$chrm == "chrI"] = 1
variants$chrm[variants$chrm == "chrII"] = 2
variants$chrm[variants$chrm == "chrIII"] = 3
variants$chrm[variants$chrm == "chrIV"] = 4
variants$chrm[variants$chrm == "chrV"] = 5
variants$chrm[variants$chrm == "chrVI"] = 6
variants$chrm[variants$chrm == "chrVII"] = 7
variants$chrm[variants$chrm == "chrVIII"] = 8
variants$chrm[variants$chrm == "chrIX"] = 9
variants$chrm[variants$chrm == "chrX"] = 10
variants$chrm[variants$chrm == "chrXI"] = 11
variants$chrm[variants$chrm == "chrXII"] = 12
variants$chrm[variants$chrm == "chrXIII"] = 13
variants$chrm[variants$chrm == "chrXIV"] = 14
variants$chrm[variants$chrm == "chrXV"] = 15
variants$chrm[variants$chrm == "chrXVI"] = 16
variants = separate(data = variants, col = info, into = c('position', 'variant'), sep = '_')

#now add up how many variants within each CI. Not sure how to do that yet
  #for each QTL
  #subset variants by chrm
  #variant position < CI right and > CI left
  #if position of variant is within CI, then add 1 to a count
  #record total count


#new table
newQTLtable2 = newQTLtable
newQTLtable2$variantcount = 0
newQTLtable2 = newQTLtable2[0,]

#count how many variants (including non-coding) are within the CI
for(qtl1 in 1:nrow(newQTLtable)){
  thischrm = variants[variants$chrm == newQTLtable[qtl1, "chr"],] #variants on the same chromosome as this QTL
  withinCI = thischrm[as.numeric(thischrm$position) < newQTLtable[qtl1, "right_Index"], ]
  withinCI = withinCI[as.numeric(withinCI$position) > newQTLtable[qtl1, "left_Index"], ]
  count = nrow(withinCI) 
  newrow = newQTLtable[qtl1,]
  newrow$variantcount = count
  newQTLtable2 = rbind(newQTLtable2, newrow)
}

#now genes
#Downloaded SGD_Chr_Features_Tab.csv on 7/9/2023 from https://yeastmine.yeastgenome.org/yeastmine/results.do?trail=%257Cquery#

#import gene info
info = read.csv("SGD_Chr_Features_Tab.csv", header = T)
#change chromosome name to number
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrI"] = 1
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrII"] = 2
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrIII"] = 3
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrIV"] = 4
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrV"] = 5
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrVI"] = 6
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrVII"] = 7
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrVIII"] = 8
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrIX"] = 9
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrX"] = 10
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrXI"] = 11
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrXII"] = 12
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrXIII"] = 13
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrXIV"] = 14
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrXV"] = 15
info$Sequence.Feature...Chromosome...Identifier[info$Sequence.Feature...Chromosome...Identifier == "chrXVI"] = 16

#remove tRNA genes
filteredinfo = info[info$Sequence.Feature...Feature.Type != "tRNA gene", ]
#remove rows with no systematic names
filteredinfo = filteredinfo[filteredinfo$Sequence.Feature...Systematic.Name != "",]
#remove rows with no Sequence.Feature...Qualifier
filteredinfo = filteredinfo[filteredinfo$Sequence.Feature...Qualifier != "",]
#same result as just filtering for all "ORFs" in Sequence.Feature...Feature.Type...

#make a column that has standard names, and if there's no standard name, it gets filled in with systemic name
info_noname = filteredinfo[filteredinfo$Sequence.Feature...Standard.Name == "",]
info_name = filteredinfo[filteredinfo$Sequence.Feature...Standard.Name != "",]

info_noname$label = info_noname$Sequence.Feature...Systematic.Name
info_name$label = info_name$Sequence.Feature...Standard.Name

filteredinfo = rbind(info_name, info_noname)

#new table with genes
QTLtable_genes = newQTLtable2
QTLtable_genes$genecount = 0
QTLtable_genes$genes = ""
QTLtable_genes = QTLtable_genes[0,]

allorfs = "" #want to get a full list of ORFs to see how many unique ones there are total

# #sanity check to see if length of CI is related to variants/number of ORFs
# QTLtable_genes2 = QTLtable_genes
# QTLtable_genes2$CIlength = 0
# QTLtable_genes2 = QTLtable_genes2[0,]

#count and list how many genes are within the CI
for(qtl1 in 1:nrow(newQTLtable2)){
  thischrm = filteredinfo[filteredinfo$Sequence.Feature...Chromosome...Identifier == newQTLtable2[qtl1, "chr"],] #genes on the same chromosome as this QTL
  withinCI = thischrm[thischrm$Sequence.Feature...Chromosome.Location...Start < newQTLtable2[qtl1, "right_Index"], ]
  withinCI = withinCI[withinCI$Sequence.Feature...Chromosome.Location...End > newQTLtable2[qtl1, "left_Index"], ]
  newrow = newQTLtable2[qtl1,]
  newrow$genecount = nrow(withinCI)
  genes = withinCI$label #gets all ORFs within the CI
  allorfs = c(allorfs, genes) #want to get a full list of ORFs to see how many unique ones there are total
  newrow$genes = paste(genes, collapse = ", ") #pastes all the genes into one string
  newrow$CIlength = newQTLtable2[qtl1, "right_Index"] -  newQTLtable2[qtl1, "left_Index"] #sanity check to see if length of CI is related to variants/number of ORFs
  QTLtable_genes2 = rbind(QTLtable_genes2, newrow)
  
}

plot(QTLtable_genes2$genecount, QTLtable_genes2$CIlength)
#number of all orfs within all the QTLs
length(unique(allorfs)) #4181 minus one for the empty value I initiated with = 4180

hist(QTLtable_genes$genecount)

write.table(QTLtable_genes,
            file = "QTLlist_revisions_4.1.25.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#reload

QTLtable_genes = read.table(file = "/Users/randia/myproteasome/2024.07.05_allHL_GxE/QTLlist_revisions_4.1.25.txt", header = T, sep = '\t')

#code below was written before revisions

#count up how many QTLs per sample
totqtls = qtls2 %>% count(reporter, environment)

write.table(totqtls,
            file = "totalQTLs8.07.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

totqtls = read.table(file = "totalQTLs8.07.24.txt", header = T, sep = '\t')


#test to see if number of QTLs for each reporter is sig diff
#reorder levels - did Asn cuz it's in the middle for # of QTLs
totqtls$reporter = factor(totqtls$reporter, levels = c('Asn', 'Phe', 'Thr', 'UFD', '4xUb', 'Rpn4'))
nb = glm.nb(n ~ reporter, data = totqtls)
summary(nb)

#test to see if number of QTLs for each environment is sig diff
#reorder levels to compare to SC
totqtls$environment = factor(totqtls$environment, levels = c('SC', 'AZC', 'LiAc', '4NQO', 'BTZ', 'Low N', 'YNB', 'Low G'))
nbenv = glm.nb(n ~ environment, data = totqtls)
summary(nbenv)


#starvation environments have most QTLs, so compare those to the rest
#make new column
totqtls$type = ifelse(totqtls$environment == 'Low N' |totqtls$environment == 'Low G'|totqtls$environment == 'YNB' , 'Starvation', 'Not Starvation')
totqtls$ubsystem = ifelse(totqtls$reporter == '4xUb' |totqtls$reporter == 'Rpn4', 'Independent', 'Dependent')
sum(totqtls[which(totqtls$type == 'Starvation'), "n"]) 
sum(totqtls[which(totqtls$ubsystem == 'Independent'), "n"]) 
#t.test(n ~ type, totqtls) #might not be normally distributed, so do wilcoxon instead
#DO WILCOXON
boxplot(n ~ type, totqtls) #starvation vs not
wilcox.test(n ~ type, data = totqtls, exact = FALSE)
boxplot(n ~ ubsystem, totqtls) #ub ind vs dep
wilcox.test(n ~ ubsystem, data = totqtls, exact = FALSE)

median(totqtls[which(totqtls$type == 'Not Starvation'), "n"]) 
median(totqtls[which(totqtls$ubsystem == 'Dependent'), "n"]) 

#use env_factor from parental data
env_factor = c("SC", "4NQO", "AZC", 'BTZ',
               'LiAc' , 'Low G' , 'Low N' , 'YNB')
totqtls$environment <-factor(totqtls$environment, levels = env_factor)

#same colors as parental plots
sc_env_colors = c("SC" = 'black', "4NQO" = '#BC71C9', "AZC" = '#6878D2', 'BTZ' = '#00B6EB',
                  'LiAc' = '#00C094', 'Low G' = '#C49A02', 'Low N' = '#FFB512', 'YNB' = '#F85735')

#Figure 3B
##plot by reporter
ggplot(totqtls, aes(fill=environment, y=n, x=reorder(reporter, n))) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = sc_env_colors)+
  theme(axis.text=element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.title=element_text(size=15)) + #,face="bold"
  labs( x = "Reporter", y ='Number of QTLs') +
  ggtitle('QTLs per reporter and environment')+ 
  scale_y_continuous(breaks=seq(0,120,20))
#saved as 4x5 pdf


#Figure 3C
#order by environment
level_order <- c('AZC', 'LiAc', 'SC', '4NQO', 'BTZ', 'Low N', 'YNB', 'Low G') 

#change order of reporters
reporter_factor = c('Asn', 'Phe', 'Thr', 'UFD', '4xUb', 'Rpn4')
totqtls$reporter <-factor(totqtls$reporter, levels = reporter_factor)
ggplot(totqtls, aes(fill=reporter, y=n, x=factor(environment, level = level_order))) + 
  scale_fill_brewer(palette = "Dark2") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text=element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        panel.grid.minor.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.title=element_text(size=15)) +
  labs( x = "Environment", y ='Number of QTLs') +
  ggtitle('QTLs per reporter and environment')+ 
  scale_y_continuous(breaks=seq(0,80,20))
#saved as 4x6.5 pdf


#Figure 3D
#plot heatmap of total QTLs
#level_order_env <- c('BTZ', 'AZC', '4NQO', 'LiAc', 'YNB', 'Low N', 'Low G', 'SC') 
level_order_rep <- c('Asn', 'Phe', 'Thr', 'UFD', '4xUb', 'Rpn4')

#number of QTLs heatmap
ggplot(totqtls, aes(x=factor(reporter, level = level_order_rep), y=factor(environment, level = rev(env_factor)), fill = n)) +
  geom_tile(color = 'black') +
  labs(x = 'Reporter', y = 'Environment')+
  #facet_grid(~ubsystem, switch = "x", scales = "free_x", space = "free_x") +
  theme(axis.text=element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "white"),
        axis.title=element_text(size=15),
        axis.ticks = element_blank()) +
  geom_text(aes(label = n, size = 16)) +
  scale_fill_gradient(high = 'maroon', low = 'beige') +
  scale_y_discrete(position = "left", guide = guide_axis(position = 'right'))+
  scale_x_discrete(position = "top", guide = guide_axis(position = 'bottom'))
#saved 4.75x5.25 inches pdf

  #facet_grid(. ~ factor(reporter, level = reporter_factor))
#add labels for categories



  
#plot where fill is whether BY or RM increased deg
totalqtlsBYRM = qtls2 %>% dplyr::count(reporter, sign(delta_AF))
cols = c('reporter', 'AFsign', 'n')
colnames(totalqtlsBYRM) = cols

sum(totalqtlsBYRM[which(totalqtlsBYRM$AFsign == 1), "n"]) 

totalqtlsBYRM[totalqtlsBYRM$AFsign == -1, "AFsign"] = 'BY' 
totalqtlsBYRM[totalqtlsBYRM$AFsign == 1, "AFsign"] = 'RM' 

#Supplementary Figure 4B
ggplot(totalqtlsBYRM, aes(fill=as.character(AFsign), y=n, x=reorder(reporter, n))) + 
  geom_bar(position=position_dodge(), stat="identity")+ #position="stack" for stacked bars
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=15),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.title=element_text(size=15)) + #,face="bold"
  scale_fill_manual(values = c(col_by, col_rm))+
  labs( x = "Reporter", y ='Number of QTLs') +
  #ggtitle('Total QTLs by which allele increased degradation')+ 
  scale_y_continuous(breaks=seq(0,120,20)) +
  guides(fill=guide_legend(title="Allele that\nincreased\nUPS activity"))
#saved as 3x4.5


#Stats for BY v RM for reporter

#Two groups are RM increased or decreased deg. Data points are each environment, and how many QTLs increased or decreased deg for that env
totalqtlsBYRMttest = qtls2 %>% count(reporter, environment, sign(delta_AF))
cols = c('reporter', 'environment', 'AFsign', 'n')
colnames(totalqtlsBYRMttest) = cols
totalqtlsBYRMttest$alleleincreaseddeg = NA
totalqtlsBYRMttest[totalqtlsBYRMttest$AFsign == -1, "alleleincreaseddeg"] = 'BY' 
totalqtlsBYRMttest[totalqtlsBYRMttest$AFsign == 1, "alleleincreaseddeg"] = 'RM' 

#need to add rows with 0's for samples that are missing in the previous df due to being zero
totalqtlsBYRMttest = rbind(totalqtlsBYRMttest, c('UFD', 'Low N', 1, 0, 'RM'))
totalqtlsBYRMttest = rbind(totalqtlsBYRMttest, c('4xUb', 'LiAc', 1, 0, 'RM'))
#for some reason adding these rows made n a chr, even tho I added numbers
totalqtlsBYRMttest$n = as.numeric(totalqtlsBYRMttest$n)

#binomial test. This method collapses the multiple data points we have, so prob won't use
sum(totalqtlsBYRM[totalqtlsBYRM$AFsign == 1, 'n'])
binom.test(x = 229,  ## n. 'successes' - in this context, n. QTLs where RM allele increases degradation (number of pos dAF QTLs = 'success')
           n = 416,  ## n. total QTLs
           p = 0.5,  ## expected ratio: BY / RM allele of QTL equally likely to increase degradation
           alternative = "two.sided",
           conf.level = 0.95)$p.value
## 0.0442813

#six binomials for each reporter
#assuming 50-50 ratio, and 95% confidence interval
for(x in unique(totalqtlsBYRMttest$reporter)){
   df = totalqtlsBYRMttest[totalqtlsBYRMttest$reporter == x,]
   tot = sum(df$n)
   success = sum(df[df$alleleincreaseddeg == 'RM', 'n'])
   result = binom.test(x = success, n = tot)
   print(x)
   print(result$p.value)}

#binomials for categories of reporters and environments
totalqtlsBYRMttest$starve = ifelse(totalqtlsBYRMttest$environment == 'Low N' |totalqtlsBYRMttest$environment == 'Low G'|totalqtlsBYRMttest$environment == 'YNB' , 'Starvation', 'Not Starvation')
totalqtlsBYRMttest$ubsystem = ifelse(totalqtlsBYRMttest$reporter == '4xUb' |totalqtlsBYRMttest$reporter == 'Rpn4', 'Independent', 'Dependent')
for(x in unique(totalqtlsBYRMttest$starve)){
  df = totalqtlsBYRMttest[totalqtlsBYRMttest$starve == x,]
  tot = sum(df$n)
  success = sum(df[df$alleleincreaseddeg == 'RM', 'n'])
  result = binom.test(x = success, n = tot)
  print(x)
  print(c("tot:", tot, "success", success))
  print(c("percent:", success/tot))
  print(c("p-value:", result$p.value))}

for(x in unique(totalqtlsBYRMttest$ubsystem)){
  df = totalqtlsBYRMttest[totalqtlsBYRMttest$ubsystem == x,]
  tot = sum(df$n)
  success = sum(df[df$alleleincreaseddeg == 'RM', 'n'])
  result = binom.test(x = success, n = tot)
  print(x)
  print(c("tot:", tot, "success", success))
  print(c("percent:", success/tot))
  print(c("p-value:", result$p.value))}




#one paired ttest for all reporters
#need to make sure the 'paired' data points are matched up/in order

totalqtlsBYRMttestpaired = totalqtlsBYRMttest[c('reporter', 'environment', 'n', 'alleleincreaseddeg')]
totalqtlsBYRMttestpaired = pivot_wider(totalqtlsBYRMttestpaired, values_from = n, names_from = alleleincreaseddeg)

t.test(totalqtlsBYRMttestpaired$RM, totalqtlsBYRMttestpaired$BY, paired = TRUE)

#six t-tests, one per reporter.
totalqtlsBYRMttestpaired %>% group_by(reporter) %>%
  summarise(p.value = t.test(RM, BY, paired = TRUE)$p.value) %>%
  ungroup()


# #plot coloring bars by environment - split bars in to -1 and 1. Prob won't use
# ggplot(totalqtlsBYRMttest, aes(y=n, x=reorder(AFsign, n))) +
#   geom_bar(aes(fill = environment), position='stack', stat="identity")+ #position="stack" for stacked bars
#   theme(axis.text=element_text(size=14),
#         panel.background = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
#         axis.line = element_line(colour = "grey"),
#         axis.title=element_text(size=15)) + #,face="bold"
#   facet_grid(. ~ reporter)



##plot by environment
#reorder didn't work with one Low N sample missing (4xUb)
# level_order <- c('AZC', 'LiAc', 'SC', '4NQO', 'BTZ', 'Low N', 'YNB', 'Low G') 
# ggplot(totqtls, aes(fill=reporter, y=n, x=factor(environment, level = level_order))) + 
#   scale_fill_brewer(palette = "Dark2") +
#   geom_bar(position="stack", stat="identity")+
#   theme(axis.text=element_text(size=14),
#         panel.background = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
#         panel.grid.minor.y = element_line(color="light grey", linewidth = 0.2),
#         axis.line = element_line(colour = "grey"),
#         axis.title=element_text(size=15)) +
#   labs( x = "Environment", y ='Total QTLs') +
#   ggtitle('QTLs per reporter and environment')+ 
#   scale_y_continuous(breaks=seq(0,80,20)) +
#   scale_fill_brewer(palette="RdBu")


#plot whether BY or RM increased deg
totalqtlsBYRMenv = qtls2 %>% dplyr::count(environment, sign(delta_AF))
cols2 = c('environment', 'AFsign', 'n')
colnames(totalqtlsBYRMenv) = cols2

totalqtlsBYRMenv[totalqtlsBYRMenv$AFsign == -1, "AFsign"] = 'BY' 
totalqtlsBYRMenv[totalqtlsBYRMenv$AFsign == 1, "AFsign"] = 'RM' 

#Supplementary Figure 4C
ggplot(totalqtlsBYRMenv, aes(fill=as.character(AFsign), y=n, x=reorder(environment, n))) + 
  geom_bar(position=position_dodge(), stat="identity")+ #position="stack" for stacked bars
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=15),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.title=element_text(size=15)) + #,face="bold"
  scale_fill_manual(values = c(col_by, col_rm))+
  labs( x = "Environment", y ='Number of QTLs') +
  #ggtitle('Total QTLs by which allele increased degradation')+ 
  scale_y_continuous(breaks=seq(0,120,10)) +
  guides(fill=guide_legend(title="Allele that\nincreased\nUPS activity"))
#saved as 3x6


#Stats for BY v RM for environment

#8 t-tests, one per environment. Two groups are RM increased or decreased deg. Data points are each reporter, and how many QTLs increased or decreased deg for that reporter

totalqtlsBYRMttestpaired %>% group_by(environment) %>%
  summarise(p.value = t.test(RM, BY, paired = TRUE)$p.value) %>%
  ungroup()



#do eight binomials - per environment
#assuming 50-50 ratio, and 95% confidence interval
for(x in unique(totalqtlsBYRMttest$environment)){
  df = totalqtlsBYRMttest[totalqtlsBYRMttest$environment == x,]
  tot = sum(df$n)
  success = sum(df[df$alleleincreaseddeg == 'RM', 'n'])
  result = binom.test(x = success, n = tot)
  print(x)
  print(result$p.value)}


# #plot coloring bars by environment - split bars in to -1 and 1. Prob won't use
# ggplot(totalqtlsBYRMttest, aes(y=n, x=reorder(AFsign, n))) + 
#   geom_bar(aes(fill = environment), position='stack', stat="identity")+ #position="stack" for stacked bars
#   theme(axis.text=element_text(size=14),
#         panel.background = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
#         axis.line = element_line(colour = "grey"),
#         axis.title=element_text(size=15)) + #,face="bold"
#   facet_grid(. ~ reporter) 

##Done describing numbers of QTLs per sample##



####GxE aka overlapping analysis####

#Find QTLs that overlap SC and another ENVIRONMENT, WITHIN REPORTERS
##Just comparing to SC##

#add column with reporter and environment
qtls$rep = paste0(qtls$reporter, '_in_', qtls$environment)

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
      ##load the QTL table for
      ## replicate 1 as 'rep_1_peaks_table' 
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
            { ##min because there could be more than one peak within 100kb. Abs because there could be a QTL to the left or the right
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
                merge_frame_new = rbind(merge_frame_new, df) #rra changed i to j
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
sc_combos = read.table(file = "5_combos_QTL_overlaps.txt", header = T, sep = '\t') 

overlaps = find_overlaps(sc_combos, qtls, consider_sign = FALSE)
#separate out the dfs from the list
overlaps_QTLinfo = overlaps[[1]]
overlaps_merge = overlaps[[2]] #849 overlapping pairs

#now add 'total reps' and 'direction of effect' 
overlaps_merge$total_reps = overlaps_merge$num_of_reps_1 + overlaps_merge$num_of_reps_2
overlaps_merge$true_means_same_dir = overlaps_merge$rep_1_delta_AF *overlaps_merge$rep_2_delta_AF > 0

#change N/As to absent, cuz I wanna look at present/absent
overlaps_merge[is.na(overlaps_merge)] = "Absent"


#save QTL info table
write.table(overlaps_QTLinfo,
            file = "overlappingQTLinfo_24.07.24_savingQTLsin1rep.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#save combined qtls
write.table(overlaps_merge,
            file = "overlappingQTLs_forGxE_24.07.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#can reload here
#overlaps_not_filtered = read.table(file = "overlappingQTLs_forGxE_24.07.15.txt", header = TRUE, sep = '\t')


#filter out rows that don't have enough QTLs
#remove total QTLs that = 1
overlaps_merge_filtered = overlaps_merge[overlaps_merge$total_reps != 1,] #547 overlapping pairs
#remove rows that show overlaps, but there is only 1 QTL per sample
#we don't care about rows that say 'true' or 'false' but total qtls = 2, because that means the QTL was only found in 1/2 replicates for each sample
overlaps_merge_filtered$remrow = paste0(overlaps_merge_filtered$total_reps, overlaps_merge_filtered$true_means_same_dir)
overlaps_merge_filtered = overlaps_merge_filtered[overlaps_merge_filtered$remrow != '2TRUE',  ] #512
overlaps_merge_filtered = overlaps_merge_filtered[overlaps_merge_filtered$remrow != '2FALSE',  ] #507

#2Absent means present in 2 reps, absent in 2 reps (P/A GxE),
#3FALSE means one direction in two reps, and the opposite direction in one rep (sign change GxE),
#3TRUE means one direction in two reps, and the same direction in one rep (no GxE),
#4FALSE means one direction in two reps, and the opposite direction in two reps (sign change GxE),
#4TRUE means one direction in two reps, and the same direction in two reps (no GxE)

print(c('Total overlapping pairs:', nrow(overlaps_merge_filtered)), quote = F) #should be 507 for my data

#save overlapping qtls
write.table(overlaps_merge_filtered,
            file = "overlappingQTLs_forGxE_24.07.24_filtered.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


###reload point###
overlaps_merge_filtered = read.table(file = "overlappingQTLs_forGxE_24.07.24_filtered.txt", header = TRUE, sep = '\t')




#Look at sign changes and rank based on LOD and dAF. Including overlaps, even if only in 3 replicates total (was only 4 out of 17 sign changes)
strongGxE = overlaps_merge_filtered[overlaps_merge_filtered$remrow == '4FALSE',]
false3 = overlaps_merge_filtered[overlaps_merge_filtered$remrow == '3FALSE',]
strongGxE = rbind(strongGxE, false3) #these are the Sign Change GxE QTLs. n = 17
strongGxE$lodAdded = strongGxE$rep_1_LOD + as.numeric(strongGxE$rep_2_LOD)
strongGxE$dAFAdded = abs(strongGxE$rep_1_delta_AF) + abs(as.numeric(strongGxE$rep_2_delta_AF)) #needs to be abs since some dAF are neg

#can see if any stand out in particular
hist(strongGxE$lodAdded, breaks = 50)
hist(strongGxE$dAFAdded, breaks = 50)

write.table(strongGxE,
            file = "strongestGxE24.07.14.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

strongGxE = read.table(file = "strongestGxE24.07.14.txt", header = TRUE, sep = '\t', check.names=FALSE)

# #Want to see where LOD and dAF lay for the top sign change loci within all QTLs. The top 3 pairs are made of 5 QTLs
# hist(qtls2$LOD, breaks = 100)
# abline(v = c(47.3,134.5,6.1,106.8,22,7), col = 'red', lwd = 0.5)
#
# #deciding if I should rank out of qtls or qtls2. Can rank out of qtls2, since QTLs present in 1 rep overlaps well with those in 2 reps
# ggplot(qtls, aes(x = LOD)) +
#        geom_histogram(aes(color = as.character(number_of_reps),
#                           fill = as.character(number_of_reps)),  alpha = 0.4, position = 'stack', bins = 100)
#
# boxplot(rankLOD ~ as.character(number_of_reps), data = qtls)
##I think it makes more sense to report on overlapping pairs
#
# hist(abs(qtls2$delta_AF), breaks = 100)
# abline(v = c(0.286,0.493,0.101,0.44,0.194), col = 'red', lwd = 0.5)

#Want to see where the 5 QTL's confidence interval widths lay for all QTLs. The top 3 pairs are made of 5 QTLs
qtls$CIwidth = qtls$right_Index - qtls$left_Index
qtls$rankCIwidth = rank(qtls$CIwidth)

#deciding if I should rank out of qtls or qtls2. Can rank out of qtls2, since QTLs present in 1 rep overlaps well with those in 2 reps
#can just color by both
#only showing the 3 QTLs that map to the HAP1 locus
ggplot(qtls, aes(x = CIwidth)) +
       geom_histogram(aes(color = as.character(number_of_reps),
                          fill = as.character(number_of_reps)),  alpha = 0.4, position = 'stack', bins = 100) +#position is stack or fill or dodge
  geom_vline(xintercept = c(25900, 23550, 40600, 90700, 22150), linewidth = 0.25)

#want to rank the 5 LODs and dAFs in qtls and see where they are in the 694 QTLs
qtls$rankLOD = rank(qtls$LOD)
qtls$rankdAF = rank(abs(qtls$delta_AF))

#make a table to of the top QTLs
topoverlaps = qtls[qtls$LOD == 47.285,]
topoverlaps = rbind(topoverlaps, qtls[qtls$LOD == 134.485,])
topoverlaps = rbind(topoverlaps, qtls[qtls$LOD == 22.745,])
topoverlaps = rbind(topoverlaps, qtls[qtls$LOD == 6.105,])
topoverlaps = rbind(topoverlaps, qtls[qtls$LOD == 106.795,])
#what % is the rank
topoverlaps = transform(topoverlaps, rankLODp = 1-(rankLOD / 694)) #Higher LOD means stronger
topoverlaps = transform(topoverlaps, rankdAFp = 1-(rankdAF / 694)) #Higher dAF means stronger
topoverlaps = transform(topoverlaps, rankdCIwp = rankCIwidth / 694) #smaller CI width means stronger


#save this additional info
write.table(topoverlaps,
            file = "strongestGxE24.07.26_more_info_outof694qtls.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#Want to see where added LOD and dAF lay for the top sign change pairs. The top 3 pairs are made of 5 QTLs
overlapping = overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir != 'Absent', ] #253 (same as 'notabsent' below)
overlapping$lodAdded = overlapping$rep_1_LOD + as.numeric(overlapping$rep_2_LOD)
overlapping$dAFAdded = abs(overlapping$rep_1_delta_AF) + abs(as.numeric(overlapping$rep_2_delta_AF)) #needs to be abs since some dAF are neg
hist(overlapping$lodAdded, breaks = 75)
abline(v = c(181.77, 112.9, 70.03), col = 'red', lwd = 0.5)
hist(overlapping$dAFAdded, breaks = 75)
abline(v = c(0.779,0.541,0.480), col = 'red', lwd = 0.5)

#want to rank the 3 strongest overlaps, where they are in all 253 overlaps
overlapping$rankLOD = rank(overlapping$lodAdded)
overlapping$rankdAF = rank(abs(overlapping$dAFAdded))

#what % is the rank
topoverlappairs = overlapping[overlapping$lodAdded == 181.77,]
topoverlappairs = rbind(topoverlappairs, overlapping[overlapping$lodAdded == 112.9,])
topoverlappairs = rbind(topoverlappairs, overlapping[overlapping$lodAdded == 70.03,])

topoverlappairs = transform(topoverlappairs, rankLODp = 1 - (rankLOD / 253))
topoverlappairs = transform(topoverlappairs, rankdAFp = 1 - (rankdAF / 253))

#save this additional info
write.table(topoverlappairs,
            file = "strongestGxE24.07.26_more_info.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#showing where the 17 sign changes are within all 253 overlapping pairs of QTLs, in terms of LOD and dAF
myhist = data.frame(overlapping$lodAdded, overlapping$true_means_same_dir, overlapping$dAFAdded)
ggplot(myhist, aes(x = overlapping.lodAdded)) +
       geom_histogram(aes(color = overlapping.true_means_same_dir,
                          fill = overlapping.true_means_same_dir,  alpha = 0.4),  
                      position = 'stack', bins = 100)
ggplot(myhist, aes(x = overlapping.dAFAdded)) +
  geom_histogram(aes(color = overlapping.true_means_same_dir,
                     fill = overlapping.true_means_same_dir,  alpha = 0.4),  
                 position = 'stack', bins = 100)



##Want to see how many QTLs make up the 17 sign change pairs, because it's not 34! It's fewer!

strongGxE$key1 = paste0(strongGxE$reporter_1, strongGxE$environment_1, strongGxE$chrm, ceiling(strongGxE$rep_1_LOD), round(strongGxE$rep_1_delta_AF, digits = 3), strongGxE$rep_1_left_Index, strongGxE$rep_1_max_Index, strongGxE$rep_1_right_Index)
#have to force as.numeric for when it says 'absent'
strongGxE$key2 = paste0(strongGxE$reporter_2, strongGxE$environment_2, strongGxE$chrm, ceiling(as.numeric(strongGxE$rep_2_LOD)), round(as.numeric(strongGxE$rep_2_delta_AF), digits = 3), strongGxE$rep_2_left_Index, strongGxE$rep_2_max_Index, strongGxE$rep_2_right_Index)

strongkeys = c(strongGxE$key1, strongGxE$key2)
length(unique(strongkeys)) #29 QTLs
duplicated(strongkeys)
strongkeysdf = data.frame(strongkeys, duplicated(strongkeys))
#df of the 29 QTLs in the sign change pairs
strongkeysdf = strongkeysdf[strongkeysdf$duplicated.strongkeys. == FALSE,]
qtlssc = qtls
qtlssc$key = paste0(qtlssc$reporter, qtlssc$environment, qtlssc$chr, ceiling(qtlssc$LOD), round(qtlssc$delta_AF, digits = 3), qtlssc$left_Index, qtlssc$max_Index, qtlssc$right_Index)

signchangeQTLs = qtlssc[0,]
signchangeQTLs = qtlssc[qtlssc$key %in% strongkeysdf$strongkeys,]

write.table(signchangeQTLs,
            file = "signchangeQTLs_29.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#next I plotted 29 sign change QTLs where there was more than one pair to a chromosome on UCSC genome browser and saved screenshots
#see /Users/randia/myproteasome/2024.07.05_allHL_GxE/peaks for UCSC_sign changes.xlsx


#GxE stat plots

####Want to know per REPORTER:####
# #"2Absent" "3FALSE" "4FALSE" = GxE 
# #"3TRUE" "4TRUE" = overlapping, but not GxE


#need to separate out 'Absent' since the info is set up differently
#1/3 DF#
notabsent = overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir != 'Absent',] #253 (same as 'overlapping' above)
absent = overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir == 'Absent',] #254

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
notabsent_gxe = gxestats_fun(notabsent, "reporter_1")

#now do P/A
#need to separate out if the QTL was present in the environment (and absent in SC) or absent in the environment (and present in SC)
#'environment_1' is where the QTL was present. repB2 has the sample (rep_env) where the QTL was absent - need to separate out that environment
sc_present = absent[absent$environment_1 == 'SC', ] #this only works if SC is always the first environment in combolist!
#2/3 DF#
sc_absent = absent[absent$environment_1 != 'SC', ] #the environment listed had the QTL, and the QTL was absent in SC
#3/3 DF#
sc_present$environment2absent = str_split_i((sc_present$repB2), '_', 3) #shows the environment that the QTL was absent in

#run function
scpresent_gxe = gxestats_fun(sc_present, 'reporter_1')
scabsent_gxe = gxestats_fun(sc_absent, 'reporter_1')


#need to combine the 3 dfs, so reorganize with the data I want to keep, and rename the columns
keepsRep = c('Reporter', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
gxestats1 = notabsent_gxe[keepsRep]
scpresent_gxe$'Present in SC' = scpresent_gxe$`Presence/Absense GxE (2Absent)`
keeps2Rep = c('Reporter', 'Present in SC')
gxestats2 = scpresent_gxe[keeps2Rep]
scabsent_gxe$'Absent in SC' = scabsent_gxe$`Presence/Absense GxE (2Absent)`
keeps3Rep = c('Reporter', 'Absent in SC')
gxestats3 = scabsent_gxe[keeps3Rep]

gxestats_done = merge(gxestats1,gxestats2, by = 'Reporter')
gxestats_done = merge(gxestats_done,gxestats3, by = 'Reporter')
gxestats_done$`Total Overlapping Pairs` = gxestats_done$`Sign change GxE (4FALSE)` + gxestats_done$`Sign change GxE (3FALSE)` + gxestats_done$`Overlapping but NO GxE (4TRUE)` + gxestats_done$`Overlapping but NO GxE (3TRUE)` + gxestats_done$`Absent in SC` + gxestats_done$`Present in SC`


#GxE stat plots for numbers separated by reporter
#add up sign change columns
gxestats_done$'Sign Change GxE' = gxestats_done$`Sign change GxE (4FALSE)` + gxestats_done$`Sign change GxE (3FALSE)` 
gxestats_done$'Total GxE' = gxestats_done$`Sign change GxE (4FALSE)` + gxestats_done$`Sign change GxE (3FALSE)` + gxestats_done$`Absent in SC` + gxestats_done$`Present in SC`


write.table(gxestats_done,
            file = "GxEstatsReporter24.07.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


#reload
gxestats_done = read.table(file = "GxEstatsReporter24.07.24.txt", header = TRUE, sep = '\t', check.names=FALSE)





#format df for plotting (did not combine the 3 and 4 replicates of No GxE here, like I did for the reporters)
gxestatslonger = pivot_longer(gxestats_done, cols = c('Absent in SC', 'Present in SC', 'Sign Change GxE', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)'),
                                 names_to = 'GxE Status',
                                 values_to = 'Overlapping Pairs')

#change labels for keys
gxestatslonger[gxestatslonger == "Overlapping but NO GxE (4TRUE)"] = "No GxE (4 samples)"
gxestatslonger[gxestatslonger == "Overlapping but NO GxE (3TRUE)"] = "No GxE (3 samples)"

#make factor for fill column so that we can reorder the stacks
#gxestatslonger$`GxE Status` = factor(gxestatslonger$`GxE Status`, levels = c('Overlapping but NO GxE (3TRUE)', 'Overlapping but NO GxE (4TRUE)', 'Sign Change GxE', 'Absent in SC', 'Present in SC'))
gxestatslonger$`GxE Status` = factor(gxestatslonger$`GxE Status`, levels = c('No GxE (3 samples)', 'No GxE (4 samples)', 'Sign Change GxE', 'Absent in SC', 'Present in SC'))

#Figure 4C
#plot (other plot reordered by environment is on line ~1133)
ggplot(gxestatslonger, aes(fill=factor(`GxE Status`), y=`Overlapping Pairs`, x=reorder(Reporter, `Total Overlapping Pairs`))) + 
  geom_bar(position="stack", stat="identity")+ #position fill for proportion
  scale_fill_manual(values = c('light grey', 'grey', mygold, mymaroon, 'brown'))+
  labs(fill = 'GxE Status') +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.3),
        axis.line = element_line(colour = "grey"))+ 
  labs( x = "Reporter", y = 'Number of comparisons (n = 507)') +
  scale_y_continuous(breaks=seq(0,150,25), limits = c(0,150))

#simple colors for defense ppt
ggplot(gxestatslonger, aes(fill=factor(`GxE Status`), y=`Overlapping Pairs`, x=reorder(Reporter, `Total Overlapping Pairs`))) + 
  geom_bar(position="stack", stat="identity")+ #position fill for proportion
  scale_fill_manual(values = c('grey', 'grey', mygold, mymaroon, mymaroon))+
  labs(fill = 'GxE Status') +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.3),
        axis.line = element_line(colour = "grey"))+ 
  labs( x = "Reporter", y = 'Number of comparisons (n = 507)') +
  scale_y_continuous(breaks=seq(0,150,25), limits = c(0,150))



#plotting proportion
ggplot(gxestatslonger, aes(fill=factor(`GxE Status`), y=`Overlapping Pairs`, x=reorder(Reporter, `Total Overlapping Pairs`))) + 
  geom_bar(position="fill", stat="identity")+ #position fill for proportion
  scale_fill_manual(values = c('light grey', 'grey', mygold, mymaroon, 'brown'))+
  labs(fill = 'GxE Status') +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.3),
        axis.line = element_line(colour = "grey"))+ 
  labs( x = "Reporter", y= 'Proportion of Overlapping Pairs') 


#heatmap of proportion of GxE pairs
env_factor_noSC = c("4NQO", "AZC", 'BTZ',
               'LiAc' , 'Low G' , 'Low N' , 'YNB')

#Supplementary Figure 5
ggplot(stats_done, aes(x=factor(rep, level = level_order_rep), y=factor(environment, level = rev(env_factor_noSC)), fill = proportionGxE)) +
  geom_tile(color = 'black') +
  labs(x = 'Reporter', y = 'Environment')+
  #facet_grid(~ubsystem, switch = "x", scales = "free_x", space = "free_x") +
  theme(axis.text=element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "white"),
        axis.title=element_text(size=15),
        axis.ticks = element_blank()) +
  geom_text(aes(label = round(proportionGxE, digits = 2), size = 14)) +
  scale_fill_gradient(high = 'maroon', low = 'beige') +
  scale_y_discrete(position = "left", guide = guide_axis(position = 'right'))+
  scale_x_discrete(position = "top", guide = guide_axis(position = 'bottom'))
#saved 4.75x5.75 inches pdf



#make a pie chart for totals
piedf = gxestats_done %>% adorn_totals('row')
piedf = tail(piedf, n=1)
keepspie = c('Present in SC', 'Absent in SC', 'Sign Change GxE', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
piedf = piedf[keepspie]
piedf = as.data.frame(t(piedf))
piedf$GxE_Status = rownames(piedf)
pienames = c('Overlapping Pairs', 'GxE Status')
colnames(piedf) = pienames

#make factor for fill column so that we can reorder the stacks
piedf$`GxE Status` = factor(piedf$`GxE Status`, levels = c('Overlapping but NO GxE (3TRUE)', 'Overlapping but NO GxE (4TRUE)', 'Absent in SC', 'Present in SC','Sign Change GxE'))


ggplot(piedf, aes(fill=`GxE Status`, y=`Overlapping Pairs`, x="")) + 
  geom_bar(width = 1, stat="identity", color = NA)+ #position = stack for non-percentage
  scale_fill_manual(values = rev(c(mygold, 'brown', mymaroon, 'grey', 'light grey')))+
  coord_polar("y", start=0) +
  geom_text(aes(label = `Overlapping Pairs`, x = 1.3), position = position_stack(vjust = 0.5))+
  theme(axis.line = element_blank()) +
  theme_void()

## Stats for GxE per Reporter

#test to see if amount of GxE across reporters is different
#gxe stats per SAMPLE (n = 47)
#initiate GxE function- use for all gxestats functions - can change 'Reporter' to 'Environment'
statscols = c('Sample', 'Total Overlapping Pairs', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Presence/Absense GxE (2Absent)','Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
gxestats = data.frame(matrix(ncol = 7, nrow = 0))
colnames(gxestats) = statscols


#Next, want to know gxe per SAMPLE in order to actually run stat tests

#run function for 3 overlapping dfs
#notabsent, sc_present, sc_absent

stat_notabsent = gxestats_fun(notabsent, 'rep2') 
stat_sc_present = gxestats_fun(sc_present, 'repB2')
stat_sc_absent = gxestats_fun(sc_absent, 'rep1')

#need to combine the 3 dfs, so reorganize with the data I want to keep, and rename the columns
keepsRep = c('Sample', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
stats1 = stat_notabsent[keepsRep]
stat_sc_present$'Present in SC' = stat_sc_present$`Presence/Absense GxE (2Absent)`
keeps2Rep = c('Sample', 'Present in SC')
stats2 = stat_sc_present[keeps2Rep]
stat_sc_absent$'Absent in SC' = stat_sc_absent$`Presence/Absense GxE (2Absent)`
keeps3Rep = c('Sample', 'Absent in SC')
stats3 = stat_sc_absent[keeps3Rep]

#not all samples were present in all the df's, so need to make rows with 0's for each sample that isn't there
samples = unique(qtls$rep)
for(s in samples){
  row = stats1[stats1$Sample == s,]
  if(nrow(row)>0){
    print(c(s, 'is present'))
  } else {
    newrow = stats1[0,]
    newrow = data.frame(s, 0, 0, 0, 0)
    colnames(newrow) = colnames(stats1)
    stats1 = rbind(stats1, newrow)
  }
}
#stats2
for(s in samples){
  row = stats2[stats2$Sample == s,]
  if(nrow(row)>0){
    print(c(s, 'is present'))
  } else {
    newrow = stats2[0,]
    newrow = data.frame(s, 0)
    colnames(newrow) = colnames(stats2)
    stats2 = rbind(stats2, newrow)
  }
}
#stats3
for(s in samples){
  row = stats3[stats3$Sample == s,]
  if(nrow(row)>0){
    print(c(s, 'is present'))
  } else {
    newrow = stats3[0,]
    newrow = data.frame(s, 0)
    colnames(newrow) = colnames(stats3)
    stats3 = rbind(stats3, newrow)
  }
}


stats_done = merge(stats1,stats2, by = 'Sample')
stats_done = merge(stats_done,stats3, by = 'Sample')
stats_done$`Total Overlapping Pairs` = stats_done$`Sign change GxE (4FALSE)` + stats_done$`Sign change GxE (3FALSE)` + stats_done$`Overlapping but NO GxE (4TRUE)` + stats_done$`Overlapping but NO GxE (3TRUE)` + stats_done$`Absent in SC` + stats_done$`Present in SC`
#double checked totals for each column, and they are all correct

stats_done$rep = str_split_i(stats_done$Sample, '_', 1)
stats_done$environment = str_split_i(stats_done$Sample, '_', 3)

#now switch to long form, and total GxE and no GxE
stats_done$'Total GxE' = stats_done$`Sign change GxE (4FALSE)` + stats_done$`Sign change GxE (3FALSE)` + stats_done$`Absent in SC` + stats_done$`Present in SC`
stats_done$'No GxE' = stats_done$`Overlapping but NO GxE (3TRUE)` + stats_done$`Overlapping but NO GxE (4TRUE)`

#all the SC rows are empty (cuz comparing to SC)

write.table(stats_done,
            file = "GxEstatsSample24.07.31.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

###reload point###
stats_done = read.table(file = "GxEstatsSample24.07.31.txt", header = TRUE, sep = '\t', check.names=FALSE)


stats_done_long = pivot_longer(stats_done, cols = c('Total GxE', 'No GxE'),
                                 names_to = 'GxE Status',
                                 values_to = 'Overlapping Pairs')

#this makes all SC rows blank, because everything is in relation to SC, so remove those rows
stats_done_long = stats_done_long[stats_done_long$environment != 'SC', ] #82 data points


#Amount of GxE stats##

#the SC rows just have zeros, so remove those (everything is compared to SC)
stats_done_glm = stats_done[stats_done$environment != 'SC', ]
#want to calculate % so it's normalized against total number of QTLs
stats_done_glm$percentGxE = stats_done_glm$`Total GxE` / stats_done_glm$`Total Overlapping Pairs`
summary(glm.nb(percentGxE ~ rep, data = stats_done_glm))

#try it again with amount of GxE (disregarding no GxE)
summary(glm.nb(`Total GxE` ~ rep, data = stats_done_glm))


#ttest across all data - comparing amount of GxE vs No GxE
t.test(`Overlapping Pairs` ~ `GxE Status`, stats_done_long, paired = TRUE)

#6 ttests - one per reporter, to see if there is a sig diff between GxE and No GxE per reporter
stats_done_long %>% group_by(rep) %>%
  summarise(p.value = t.test(`Overlapping Pairs` ~ `GxE Status`, paired = TRUE)$p.value) %>%
  ungroup()

#next test we're trying is *logistic regression* - glm( family = 'binomial')
#for glm, each row needs to be 1 aka 'hit' aka 'GxE' or a 0 aka 'miss' aka 'No GxE'
  #already have a df where each row is a 'close pair of QTLs' (formerly 'overlapping pair')
  #so just need to assign 1 or 0 to the row

overlaps_merge_filtered$gxeis1 = NA
overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir == 'Absent', "gxeis1"] = 1 
overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir == 'FALSE', "gxeis1"] = 1 
overlaps_merge_filtered[overlaps_merge_filtered$true_means_same_dir == 'TRUE', "gxeis1"] = 0

mylogreg = glm(gxeis1 ~ reporter_1, data = overlaps_merge_filtered, family = binomial())
summary(mylogreg)

#need to run multiple, so we get all pairwise combos
#make a list of pairwise combos to record data in Excel and my Results Ppt
reporters = unique(overlaps_merge_filtered$reporter_1)
rep_combos = combn(reporters, 2, simplify = F)


#now change factor level order to get other combinations. Changed the level order by hand and re-did the logistic regression to get through all the iterations
overlaps_merge_filtered$reporter_1 = factor(overlaps_merge_filtered$reporter_1, levels = c("Thr",  "UFD", "4xUb", "Asn", "Phe", "Rpn4"))
mylogreg2 = glm(gxeis1 ~ reporter_1, data = overlaps_merge_filtered, family = binomial())
summary(mylogreg2)

##done looking at QTL GxE per reporter## (Just comparing to SC)


#adding categories
#did reporter categories and environment categories together here

stats_done$ubsystem = ifelse(stats_done$rep == '4xUb' |stats_done$rep == 'Rpn4', 'Independent', 'Dependent')
#stats_done$chem = ifelse(stats_done$environment == '4NQO' |stats_done$environment == 'AZC'|stats_done$environment == 'BTZ' , 'Chemical', 'Starvation')
#stats_done[stats_done$environment == 'LiAc', 'chem'] <- 'Other'
stats_done = stats_done[stats_done$environment != 'SC', ] #everything is compared to SC, so this was just 0's

#want to see if there's a sig diff in *proportion* of GxE btwn the categories
#get proportion of GxE to normalize for amount of QTLs for that category
stats_done$proportionGxE = stats_done$`Total GxE` / (stats_done$`No GxE` + stats_done$`Total GxE`) #oops, had a column with the total already
boxplot(proportionGxE ~ ubsystem, stats_done)
wilcox.test(proportionGxE ~ ubsystem, stats_done, exact = FALSE)

boxplot(proportionGxE ~ chem, stats_done) #other is just LiAc

stats_done$starve = ifelse(stats_done$environment == 'Low N' |stats_done$environment == 'Low G'|stats_done$environment == 'YNB' , 'Starvation', 'Not Starvation')
boxplot(proportionGxE ~ starve, stats_done)
wilcox.test(proportionGxE ~ starve, stats_done, exact = FALSE)

#Now just for PA
#it's a proportion, putting sign change into non-P/A, with the non GxE pairs
stats_done$PA = stats_done$`Present in SC` + stats_done$`Absent in SC`
stats_done$PAproportion = stats_done$PA / stats_done$`Total Overlapping Pairs`
stats_done$PAproportion2 = stats_done$PA / (stats_done$`No GxE` + stats_done$PA) #leaving out sign change pairs
wilcox.test(PAproportion ~ ubsystem, stats_done, exact = FALSE)
wilcox.test(PAproportion2 ~ starve, stats_done, exact = FALSE)

#now for the 17 sign changes...
stats_done$signchangeall = stats_done$`Sign change GxE (3FALSE)` + stats_done$`Sign change GxE (4FALSE)`
wilcox.test(signchangeall ~ ubsystem, stats_done, exact = FALSE)

#proportion
stats_done$signchangeprop = stats_done$signchangeall/ (stats_done$`No GxE` + stats_done$`signchangeall`)
wilcox.test(signchangeprop ~ ubsystem, stats_done, exact = FALSE)


####Now want to know per ENVIRONMENT:####

#run function to get stats per environment
#update column names
statscols = c('Environment', 'Total Overlapping Pairs', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Presence/Absense GxE (2Absent)','Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
colnames(gxestats) = statscols

env_notabsent_gxe = gxestats_fun(notabsent, 'environment_2')
env_scpresent_gxe = gxestats_fun(sc_present, 'environment2absent')
env_scabsent_gxe = gxestats_fun(sc_absent, 'environment_1')


#need to combine the 3 dfs, so reorganize with the data I want to keep, and rename the columns
keeps = c('Environment', 'Sign change GxE (4FALSE)', 'Sign change GxE (3FALSE)', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)')
gxestatsenv1 = env_notabsent_gxe[keeps]
env_scpresent_gxe$'Absent in Condition; Present in SC' = env_scpresent_gxe$`Presence/Absense GxE (2Absent)`
keeps2 = c('Environment', 'Absent in Condition; Present in SC')
gxestatsenv2 = env_scpresent_gxe[keeps2]
env_scabsent_gxe$'Present in Condition; Absent in SC' = env_scabsent_gxe$`Presence/Absense GxE (2Absent)`
keeps3 = c('Environment', 'Present in Condition; Absent in SC')
gxestatsenv3 = env_scabsent_gxe[keeps3]

env_gxestats_done = merge(gxestatsenv1,gxestatsenv2, by = 'Environment')
env_gxestats_done = merge(env_gxestats_done,gxestatsenv3, by = 'Environment')
env_gxestats_done$`Total Overlapping Pairs` = env_gxestats_done$`Sign change GxE (4FALSE)` + env_gxestats_done$`Sign change GxE (3FALSE)` + env_gxestats_done$`Overlapping but NO GxE (4TRUE)` + env_gxestats_done$`Overlapping but NO GxE (3TRUE)` + env_gxestats_done$`Present in Condition; Absent in SC` + env_gxestats_done$`Absent in Condition; Present in SC`

write.table(env_gxestats_done,
            file = "GxEstatsEnvironment24.07.24.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

##reload point##
env_gxestats_done = read.table(file = "GxEstatsEnvironment24.07.24.txt", header = TRUE, sep = '\t', check.names=FALSE)


#GxE stat plots for numbers separated by environment
#add up sign change columns
env_gxestats_done$'Sign Change GxE' = env_gxestats_done$`Sign change GxE (4FALSE)` + env_gxestats_done$`Sign change GxE (3FALSE)` 
env_gxestats_done$'Total GxE' = env_gxestats_done$`Sign change GxE (4FALSE)` + env_gxestats_done$`Sign change GxE (3FALSE)` + env_gxestats_done$`Present in Condition; Absent in SC` + env_gxestats_done$`Absent in Condition; Present in SC`

#format df for plotting
gxestatsENVlonger = pivot_longer(env_gxestats_done, cols = c('Present in Condition; Absent in SC', 'Absent in Condition; Present in SC', 'Sign Change GxE', 'Overlapping but NO GxE (4TRUE)', 'Overlapping but NO GxE (3TRUE)'),
                              names_to = 'GxE Status',
                              values_to = 'Overlapping Pairs')



gxestatsENVlonger[gxestatsENVlonger == "Overlapping but NO GxE (4TRUE)"] <- "No GxE (4 samples)"
gxestatsENVlonger[gxestatsENVlonger == "Overlapping but NO GxE (3TRUE)"] <- "No GxE (3 samples)"

#make factor for fill column so that we can reorder the stacks
gxestatsENVlonger$`GxE Status` = factor(gxestatsENVlonger$`GxE Status`, levels = c('No GxE (3 samples)', 'No GxE (4 samples)', 'Sign Change GxE', 'Present in Condition; Absent in SC', 'Absent in Condition; Present in SC'))

#Figure 4D
#plot stack
ggplot(gxestatsENVlonger, aes(fill=factor(`GxE Status`), y=`Overlapping Pairs`, x=reorder(Environment, `Total GxE`))) + 
  geom_bar(position="stack", stat="identity")+ #position is fill or stack
  scale_fill_manual(values = c('light grey', 'grey', mygold, mymaroon, 'brown'))+
  labs(fill = 'GxE Status') +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.3),
        axis.line = element_line(colour = "grey"))+ 
  labs( x = "Environment", y = 'Number of comparisons (n = 507)') +
  scale_y_continuous(breaks=seq(0,110,20), limits = c(0,110))
#saved as 5x8


#plot proportion
ggplot(gxestatsENVlonger, aes(fill=factor(`GxE Status`), y=`Overlapping Pairs`, x=reorder(Environment, `Total GxE`))) + 
  geom_bar(position="fill", stat="identity")+ #position is fill or stack
  scale_fill_manual(values = c('light grey', 'grey', mygold, mymaroon, 'brown'))+
  labs(fill = 'GxE Status') +
  theme(axis.text=element_text(size=13), panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="grey", linewidth = 0.3),
        axis.line = element_line(colour = "grey"))+ 
  labs( x = "Environment", y = 'Proportion of Overlapping Pairs')

###now run stats on per environment

#7 t.tests - one per environment (that was compared to SC), to see if there is a sig diff between GxE and No GxE per reporter
stats_done_long %>% group_by(environment) %>%
  summarise(p.value = t.test(`Overlapping Pairs` ~ `GxE Status`, paired = TRUE)$p.value) %>%
  ungroup()


#negative binomial
summary(glm.nb(percentGxE ~ environment, data = stats_done_glm))

#try it again with amount of GxE (disregarding no GxE)
summary(glm.nb(`Total GxE` ~ environment, data = stats_done_glm))

#next test we're trying is *logistic regression* - glm( family = 'binomial')
#need to report the environment that is NOT SC - everything is compared to SC
  #change column name and make matching columns for the 3 absent/notabsent dfs, and combine

colnames(sc_present)[which(names(sc_present) == 'environment2absent')] <- 'env_comparedto_SC'
notabsent$env_comparedto_SC = notabsent$environment_2
sc_absent$env_comparedto_SC = sc_absent$environment_1

#now combine dfs
overlaps_merge_filtered_envs = rbind(sc_present, notabsent)
overlaps_merge_filtered_envs = rbind(overlaps_merge_filtered_envs, sc_absent)

#now record 1's and 0's for hits/misses of GxE
overlaps_merge_filtered_envs$gxeis1 = NA
overlaps_merge_filtered_envs[overlaps_merge_filtered_envs$true_means_same_dir == 'Absent', "gxeis1"] = 1 
overlaps_merge_filtered_envs[overlaps_merge_filtered_envs$true_means_same_dir == 'FALSE', "gxeis1"] = 1 
overlaps_merge_filtered_envs[overlaps_merge_filtered_envs$true_means_same_dir == 'TRUE', "gxeis1"] = 0

#run logistic regression
mylogregenv = glm(gxeis1 ~ env_comparedto_SC, data = overlaps_merge_filtered_envs, family = binomial())
summary(mylogregenv)

#need to run multiple, so we get all pairwise combos
#make a list of pairwise combos to record data in Excel and my Results Ppt
environments = unique(overlaps_merge_filtered_envs$env_comparedto_SC)
env_combos = combn(environments, 2, simplify = F)
env_combos = as.data.frame(env_combos)
env_combos = as.data.frame(t(env_combos))
rownames(env_combos) <- NULL


#now change factor level order to get other combinations. Changed the level order by hand and re-did the logistic regression to get through all the iterations
overlaps_merge_filtered_envs$env_comparedto_SC = factor(overlaps_merge_filtered_envs$env_comparedto_SC, levels = c("YNB", "Low N", "AZC", "4NQO", "Low G", "LiAc", "BTZ"))
mylogregenv2 = glm(gxeis1 ~ env_comparedto_SC, data = overlaps_merge_filtered_envs, family = binomial())
summary(mylogregenv2)




#### Is a QTL found in at least one other environment?####


#we want to use the 416 QTLs that were found in 2 replicates to see if they are found at least once in another environment for that reporter
#if a query QTL that was in SC is in a row with a pair (not P/A row), then yes, it does exist in another environment.

#output the 416 QTLs, how many times they were found, and if it's >0


#Find ALL PAIRWISE overlaps to look to see if a QTL is unique to each of the rest of environments
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

#reload the qtls df here since I edited it above!
#lines 28-40!!!


#need 'rep' for function
#may have already added above, if not using a reload point
qtls$rep = paste0(qtls$reporter, '_in_', qtls$environment)

#using combolist where SC is always in the first column
combolist_all_SCfirst = read.table(file = "combolist_allpairwise.7.22.24.txt", header = TRUE, sep = '\t')

overlaps_all = find_overlaps(combolist_all_SCfirst, qtls, consider_sign = FALSE)
#separate out the dfs from the list
overlaps_all_QTLinfo = overlaps_all[[1]]
overlaps_all_merge = overlaps_all[[2]] #3656 overlapping pairs

#now add 'total reps' and 'direction of effect' 
overlaps_all_merge$total_reps = overlaps_all_merge$num_of_reps_1 + overlaps_all_merge$num_of_reps_2
overlaps_all_merge$true_means_same_dir = overlaps_all_merge$rep_1_delta_AF *overlaps_all_merge$rep_2_delta_AF > 0

#change N/As to absent, cuz I wanna look at present/absent
overlaps_all_merge[is.na(overlaps_all_merge)] = "Absent"


#save QTL info table
write.table(overlaps_all_QTLinfo,
            file = "overlapping_ALLPAIRWISE_QTLinfo_24.07.23_savingQTLsin1rep.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#save combined qtls
write.table(overlaps_all_merge,
            file = "overlappingQTLs_ALLPAIRWISE_forGxE_24.07.23.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#can reload here
#overlaps_all_merge_not_filtered = read.table(file = "overlappingQTLs_ALLPAIRWISE_forGxE_24.07.23.txt", header = TRUE, sep = '\t')


#filter out rows that don't have enough QTLs
#remove total QTLs that = 1
overlaps_all_merge_filtered = overlaps_all_merge[overlaps_all_merge$total_reps != 1,] #2344 overlapping pairs
#remove rows that show overlaps, but there is only 1 QTL per sample
#we don't care about rows that say 'true' or 'false' but total qtls = 2, because that means the QTL was only found in 1/2 replicates for each sample
overlaps_all_merge_filtered$remrow = paste0(overlaps_all_merge_filtered$total_reps, overlaps_all_merge_filtered$true_means_same_dir)
overlaps_all_merge_filtered = overlaps_all_merge_filtered[overlaps_all_merge_filtered$remrow != '2TRUE',  ] #2242
overlaps_all_merge_filtered = overlaps_all_merge_filtered[overlaps_all_merge_filtered$remrow != '2FALSE',  ] #2221 (includes the 507 from SC)

#2Absent means present in 2 reps, absent in 2 reps (P/A GxE),
#3FALSE means one direction in two reps, and the opposite direction in one rep (sign change GxE),
#3TRUE means one direction in two reps, and the same direction in one rep (no GxE),
#4FALSE means one direction in two reps, and the opposite direction in two reps (sign change GxE),
#4TRUE means one direction in two reps, and the same direction in two reps (no GxE)

print(c('Total overlapping pairs:', nrow(overlaps_all_merge_filtered)), quote = F) #should be 2221 for my data

#save overlapping qtls
write.table(overlaps_all_merge_filtered,
            file = "overlappingQTLs_ALLPAIRWISE_forGxE_24.07.23_filtered.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')


###reload point###
overlaps_all_merge_filtered = read.table(file = "overlappingQTLs_ALLPAIRWISE_forGxE_24.07.23_filtered.txt", header = TRUE, sep = '\t')
####Done finding overlaps in all pairwise comparisons


#query: qtls2
#search df: overlaps_all_merge_filtered

#make keys
#for query (round the same as search df)
qtls2$key = paste0(qtls2$reporter, qtls2$environment, qtls2$chr, ceiling(qtls2$LOD), round(qtls2$delta_AF, digits = 3), qtls2$left_Index, qtls2$max_Index, qtls2$right_Index)


#for place to search (round the same as query df)
df = overlaps_all_merge_filtered 
df$key1 = paste0(df$reporter_1, df$environment_1, df$chrm, ceiling(df$rep_1_LOD), round(df$rep_1_delta_AF, digits = 3), df$rep_1_left_Index, df$rep_1_max_Index, df$rep_1_right_Index)
#have to force as.numeric for when it says 'absent'. Ok for warning that produced NA's since some say 'absent'
df$key2 = paste0(df$reporter_1, df$environment_2, df$chrm, ceiling(as.numeric(df$rep_2_LOD)), round(as.numeric(df$rep_2_delta_AF), digits = 3), df$rep_2_left_Index, df$rep_2_max_Index, df$rep_2_right_Index)

#initiate
founddf = qtls2[1,]
founddf$found_other_env = NA
founddf$in_other_env = NA
founddf = founddf[0,]

for (qtl in 1:nrow(qtls2)){ #QUERY
  #get query values
  key = qtls2[qtl, 'key']
  repor = qtls2[qtl, 'reporter']
  #env = qtls2[qtl, 'environment']
  #subset search df
  findin = df[df$reporter_1 == repor,]
  found1= findin[findin$key1 == key,]
  found2= findin[findin$key2 == key,]
  foundboth = rbind(found1, found2)
  founds = nrow(foundboth[foundboth$true_means_same_dir != 'Absent',]) 
  myrow = qtls2[qtl,] #make a new row to add to the found data frame
  myrow$found_other_env = founds
  myrow$in_other_env = as.character(founds >= 1)
  founddf = rbind(founddf, myrow)
}


write.table(founddf,
            file = "QTLs_FoundInOtherEnvs_24.07.31.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#reload 
#founddf = read.table(file = "QTLs_FoundInOtherEnvs_24.07.31.txt", header = TRUE, sep = '\t', check.names=FALSE)
#when you reload, remove logical vector
#founddf$in_other_env = as.character(founddf$in_other_env)

#reorganize data for plotting
found_plot = founddf %>% group_by(environment, in_other_env) %>% summarize(Freq=n())
#found_plot = founddf %>% group_by(environment, reporter, in_other_env) %>% summarize(Freq=n())
#found_plot$xaxis = paste0(found_plot$reporter, '_', found_plot$in_other_env)
#found_plot$xaxis2 = paste0(found_plot$environment, '_', found_plot$in_other_env)
#change data labels for plot

found_plot[found_plot == 'FALSE'] = 'QTLs unique to this environment'
found_plot[found_plot == 'TRUE'] = 'QTLs found in at least one other environment'

ggplot(found_plot, aes(fill = in_other_env, y=Freq, x= reorder(environment, Freq))) + #fill=environment,  #x=xaxis
  geom_bar(position=position_dodge(), stat="identity")+ #position = stack for non-percentage position_dodge()
  scale_fill_manual(values = c(mymaroon, mygold))+
  theme(axis.text=element_text(size=13),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color="light grey", linewidth = 0.2),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs( x = "Environment", y ='Number of QTLs (n = 416)')  +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25)

#get percentages for Results text
foundpercent = pivot_wider(found_plot, values_from = Freq, names_from = in_other_env)
foundpercent$percent = foundpercent$`QTLs found in at least one other environment` / (foundpercent$`QTLs found in at least one other environment` + foundpercent$`QTLs unique to this environment`)



####Hot spot analysis of the 29 QTLs in the sign change pairs####
#signchangeQTLs

#For Mahlon's 5 fine-mapped loci, I used the start and end positions of the genes
#For known hotspots, I used the peak and the start and end positions of the most likely causal genes
#rpt6 - chrm 7
#ubc6 - chrm 5
#ubr1 - chrm 7
#doa10 - chrm 9
#nta1 - chrm 10
#hap1 - chrm 12
#mkt1 - chrm 14
#ira2 - chrm 15

knownloci = read.table(file = "known_loci_7.28.24.txt", header = T)

#function to see if max Index of a QTL is within ##bp from a particular locus
inalocus = function(qtl_list, hotspotinfo, dist){
  newinfo = qtl_list #need to get the column info
  newinfo$inahotspot = "NA" #want to save a column of true/false if the peak is within a hotspot
  newinfo$whichone = "NA" #want to save which hotspot it is
  newinfo$distance = 'NA' #want to know how far the QTL peak is from the locus
  newinfo = newinfo[0,] #makes new df to save the info
  for(i in unique(qtl_list$chr)){ #go through everything by each chrm
    #subset the pairs
    overlaps <- qtl_list[qtl_list$chr == i, ] #gets all the qtls for that chromosome
    hs <- hotspotinfo[hotspotinfo$chrm == i,] #gets the loci for that chrm from the known loci
    for(j in 1:nrow(overlaps)) {
      ## find peaks w / in ## kb of ea. other:
      ## make sure to use 'abs' here!
      if(min(abs(overlaps[j,"max_Index"] - hs$position)) <= dist) {
        ## if we find something w / in ## kb, grab that row from the hotspot info that is closest
        thehs <- hs[which.min(abs(overlaps[j, "max_Index"] - hs$position)), ]
        ## use 'cbind' to make a dataframe row that contains the info from both replicates
        df = overlaps[j, ]
        df$inahotspot = 'TRUE'
        df$whichone = thehs$gene
        df$distance = min(abs(overlaps[j,"max_Index"] - hs$position))
        newinfo = rbind(newinfo, df)} else {
       #now need to add info if there isn't an overlapping hotspot
        df = overlaps[j, ]
        df$inahotspot = 'FALSE'
        df$whichone = 'NA'
        df$distance = 'NA'
        newinfo = rbind(newinfo, df)}
    }
  } 
  return(newinfo)
}

inalocus125kb = inalocus(signchangeQTLs, knownloci, 125000)
inalocus50kb = inalocus(signchangeQTLs, knownloci, 50000)
inalocus75kb = inalocus(signchangeQTLs, knownloci, 75000)
inalocus100kb = inalocus(signchangeQTLs, knownloci, 100000)

#tells us how many there is in or out of a known locus

table(inalocus50kb$inahotspot)
table(inalocus75kb$inahotspot)
table(inalocus100kb$inahotspot)
table(inalocus125kb$inahotspot)

#used a big distance (125kb), but many QTLs map much closer.
hist(as.numeric(inalocus125kb$distance), breaks = 10)
sort(as.numeric(inalocus125kb$distance))


#want to see how many of the 416 QTLs are in these loci
inalocus100kb_allQTLs = inalocus(qtls2, knownloci, 100000)
table(inalocus100kb_allQTLs$inahotspot)
table(inalocus100kb_allQTLs$whichone)

write.table(inalocus100kb_allQTLs, "inalocus100kb_allQTLs_8.28.24.txt", 
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")


#tells us how many are in/out a locus for reporter, or each category
for(i in unique(inalocus100kb_allQTLs$reporter)){
  mysub = subset(inalocus100kb_allQTLs, reporter == i)
  print(i)
  print(table(mysub$inahotspot))
}


write.table(inalocus125kb, "in_a_known_locus125kb_7.28.24.txt", 
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")


