#Written by Randi Avery as a response to reviewers

#Comparing reporter/environment combinations and whether they showed GxE in the parents to fraction of QTLs in SC and the same environment (and reporter)

library(ggplot2)

parents = read.table(file = "/Users/randia/myproteasome/reporter_characterization/environments/combined/gxesummary_2.7.25_TFTtime_loess.txt", header = T, sep = '\t')


QTLs = read.table(file = "/Users/randia/myproteasome/2024.07.05_allHL_GxE/GxEstatsSample24.07.31.txt", header = T, sep = '\t')

#for the parents, split up reporter/environment combo whether or not it showed GxE
#don't need all the data from the parent table. Also need to collapse by reporter/env combo
keeps = c("reporter", "interaction_Pval","Int_Sig", "environment.x" )
"Int_Sig"

parents = parents[keeps]

parents = unique(parents)
#Int_Sig groups whether there was GxE or not
parents$Sample = paste0(parents$reporter, "_in_", parents$environment.x)

parentsnQTLs = merge(parents, QTLs, by = "Sample")

#want fraction of QTLs that show GxE per sample
parentsnQTLs$QTLGxEfraction = parentsnQTLs$Total.GxE / (parentsnQTLs$Total.GxE + parentsnQTLs$No.GxE)


ggplot(parentsnQTLs, aes(x=Int_Sig,  y=QTLGxEfraction)) + #, fill = reporter
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "GxE between BY and RM?", y = "Fraction of QTLs with GxE") +
  geom_dotplot(binaxis='y', dotsize= 1 , alpha = 0.4, binwidth = 0.05, #aes(fill = reporter),
               stackdir = "center", 
               method="dotdensity",
               stackgroups = T,
               binpositions="all")+ #when i fill based on a group, the dots overlap, so this kinda gets around it
  #geom_dotplot(aes(fill = environment, alpha = 0.7), color = NA,  binaxis='y', stackdir='center', dotsize=.1, binwidth = 0.3)+ #color = environment.x for circle outline
  #geom_jitter(color= 'pink', size=2, alpha=0.9) +
  #geom_point(position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x = element_text(size = 10), #angle = 90, vjust = 1, hjust=1
        axis.text.y = element_text(size = 10)) +
  ggtitle("Wilcoxon p-value = 0.851")
#saved 3x4
#wilcoxon on this data
wilcox.test(QTLGxEfraction ~ Int_Sig, parentsnQTLs, exact = FALSE)



ggplot(parentsnQTLs, aes(x=Int_Sig,  y=QTLGxEfraction)) + #, fill = reporter
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "GxE between BY and RM?", y = "Fraction of QTLs with GxE") +
  geom_dotplot(aes(fill = environment),binaxis='y', dotsize= 1 , alpha = 0.4, binwidth = 0.05, 
               stackdir = "center", 
               method="dotdensity",
               stackgroups = T,
               binpositions="all")+ #when i fill based on a group, the dots overlap, so this kinda gets around it
  #geom_dotplot(aes(fill = environment, alpha = 0.7), color = NA,  binaxis='y', stackdir='center', dotsize=.1, binwidth = 0.3)+ #color = environment.x for circle outline
  #geom_jitter(color= 'pink', size=2, alpha=0.9) +
  #geom_point(position = position_jitter(w = 0.1, h = 0)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10),) +
  ggtitle("Relating parents to QTLs")
#saved 3x4

#now do the same thing but do whole counts of QTLs
#dots colored by reporter
#with GxE
ggplot(parentsnQTLs, aes(x=Int_Sig,  y=Total.GxE)) + #, fill = reporter
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "GxE between BY and RM?", y = "Total Number of QTL pairs with GxE") +
  geom_dotplot(binaxis='y', alpha = 0.4, dotsize= 1.5,#aes(fill = reporter), # alpha = 0.4, binwidth = 0.05, 
               stackdir = "center", 
               method="dotdensity",
               stackgroups = T,
               binpositions="all")+ #when i fill based on a group, the dots overlap, so this kinda gets around it
  #geom_dotplot(aes(fill = environment, alpha = 0.7), color = NA,  binaxis='y', stackdir='center', dotsize=.1, binwidth = 0.3)+ #color = environment.x for circle outline
  theme(axis.text.x = element_text(size = 10), #, angle = 90, vjust = 1, hjust=1
        axis.text.y = element_text(size = 10)) +
  ggtitle("Wilcoxon p-value = 0.64")

#wilcoxon on this data
wilcox.test(Total.GxE ~ Int_Sig, parentsnQTLs, exact = FALSE)



#not GxE
ggplot(parentsnQTLs, aes(x=Int_Sig,  y=No.GxE)) + #, fill = reporter
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "GxE between BY and RM?", y = "Total Number of QTL pairs NOT showing GxE") +
  geom_dotplot(binaxis='y', dotsize= 1.5 , alpha = 0.4, #aes(fill = reporter),  binwidth = 0.05, color = NA, 
               stackdir = "center", 
               method="dotdensity",
               stackgroups = T,
               binpositions="all")+ #when i fill based on a group, the dots overlap, so this kinda gets around it
  #geom_dotplot(aes(fill = environment, alpha = 0.7), color = NA,  binaxis='y', stackdir='center', dotsize=.1, binwidth = 0.3)+ #color = environment.x for circle outline
  theme(axis.text.x = element_text(size = 10), #, angle = 90, vjust = 1, hjust=1
        axis.text.y = element_text(size = 10)) +
  ggtitle("Wilcoxon p-value = 0.91")
#saved 3.5x4
#wilcoxon on this data
wilcox.test(No.GxE ~ Int_Sig, parentsnQTLs, exact = FALSE)

#dots colored by environment
#with GxE
ggplot(parentsnQTLs, aes(x=Int_Sig,  y=Total.GxE)) + #, fill = reporter
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "Parents showing GxE?", y = "Total Number of QTL pairs with GxE") +
  geom_dotplot(aes(fill = environment),binaxis='y', color = NA, #dotsize= 1 , alpha = 0.4, binwidth = 0.05, 
               stackdir = "center", 
               method="dotdensity",
               stackgroups = T,
               binpositions="all")+ #when i fill based on a group, the dots overlap, so this kinda gets around it
  #geom_dotplot(aes(fill = environment, alpha = 0.7), color = NA,  binaxis='y', stackdir='center', dotsize=.1, binwidth = 0.3)+ #color = environment.x for circle outline
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10),) +
  ggtitle("Relating parents to QTLs")

#not GxE
ggplot(parentsnQTLs, aes(x=Int_Sig,  y=No.GxE)) + #, fill = reporter
  #x axis can be reporter or grp2 or rep_strain, to get stacked or separate
  geom_boxplot(outlier.size = 0.3, position = 'identity') + 
  labs( x = "Parents showing GxE?", y = "Total Number of QTL pairs NOT showing GxE") +
  geom_dotplot(aes(fill = environment),binaxis='y', color = NA, #dotsize= 1 , alpha = 0.4, binwidth = 0.05, 
               stackdir = "center", 
               method="dotdensity",
               stackgroups = T,
               binpositions="all")+ #when i fill based on a group, the dots overlap, so this kinda gets around it
  #geom_dotplot(aes(fill = environment, alpha = 0.7), color = NA,  binaxis='y', stackdir='center', dotsize=.1, binwidth = 0.3)+ #color = environment.x for circle outline
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10),) +
  ggtitle("Relating parents to QTLs")





#Next analysis - relating parents to QTLs, but no GxE

rankchange = read.table(file = "results_lmer_40tests_04.29.25_TFTtime_loess_revisions.txt", header = T, sep = '\t')

#need to suss out SC values
sc = rankchange[c("reporter", "environment1", "meansscBY", "meanscRM")]
sc = unique(sc)
#now tweeze out environment values I want
env = rankchange[c("reporter", "environment2", "meanEnvBY", "meanEnvRM")]
#match column names
colnames(sc) = c("reporter", "environment", "meanBY", "meanRM")
colnames(env) = c("reporter", "environment", "meanBY", "meanRM")

parentalmeans = rbind(sc, env)
#find the diff between BY and RM
parentalmeans$BYminusRM = parentalmeans$meanBY - parentalmeans$meanRM

#now get total number of QTLs per sample
totQTLs = read.table(file = "/Users/randia/myproteasome/2024.07.05_allHL_GxE/totalQTLs8.07.24.txt", header = T, sep = '\t')
colnames(totQTLs) = c("reporter", "environment", "totQTLs")
totQTLs$sample = paste0(totQTLs$reporter,"_in_", totQTLs$environment)

parentalmeans$sample = paste0(parentalmeans$reporter,"_in_", parentalmeans$environment)

qtls = merge(parentalmeans, totQTLs)
#make sure this is right with the diff number of samples and then extra SC's
#only makes sense as absolute value difference
ggplot(qtls, aes(y = abs(BYminusRM), x = totQTLs)) +
  geom_point()+ #aes(color = ubsystem.x, shape = starve)
  ylab('UPS activity: BY mean minus RM mean') +
  xlab('Total QTLs per sample')

#colored
ggplot(qtls, aes(y = abs(BYminusRM), x = totQTLs)) +
  geom_point(aes(color = reporter), alpha = 0.7)+ #aes(color = ubsystem.x, shape = starve)
  ylab('UPS activity: absolute value of BY mean minus RM mean') +
  xlab('Total QTLs per sample')

cor.test(qtls$totQTLs, abs(qtls$BYminusRM), method = 'pearson')
