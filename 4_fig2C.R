#Written by Frank Albert in 2024 for Avery, R. R., Collins, M. A., & Albert, F. W. (2024). Genotype-by-environment interactions shape ubiquitin-proteasome system activity.

#Fig. 2C

library(tidyverse)
library(readxl)
library(tidyHeatmap)
library(wesanderson)

setwd("~/Google Drive/My Drive/UMN_projects/Randi/TFT_GxE/parentalGxE/")
dat <- read_excel("GxE_Results_in_BY_and_RM_Strains_Flow.xlsx")

# for each combination of environment and reporter, calculate the differences between BY and RM differences from SC
dat <- dat %>% mutate(
  BY_ReactionToE = medEnvBY - medscBY,
  RM_ReactionToE = medEnvRM - medscRM)

dat <- dat %>% mutate(strainDifferenceInReactionToE = BY_ReactionToE - RM_ReactionToE)
dat$reporter <- fct_relevel(dat$reporter, "Asn", "Phe", "Thr", "UFD", "4xUb", "Rpn4") 

ggplot(dat, aes(x = environment, y = strainDifferenceInReactionToE)) + 
  geom_boxplot(outliers=FALSE) + geom_jitter(aes(color=reporter), size=2) + 
  scale_color_manual(values = wes_palette("AsteroidCity2")) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() + ylab("Difference in strain response to environment\nBY - RM") + xlab("Environment")

t.test(
    dat$strainDifferenceInReactionToE[dat$environment %in% c("4NQO", "AZC", "BTZ", "LiAc")],
    dat$strainDifferenceInReactionToE[dat$environment %in% c("LowG", "LowN", "YNB")]
  )
# p = 0.24, wilcox = 0.06
# wait, we need absolute:
t.test(
  abs(dat$strainDifferenceInReactionToE[dat$environment %in% c("4NQO", "AZC", "BTZ", "LiAc")]),
  abs(dat$strainDifferenceInReactionToE[dat$environment %in% c("LowG", "LowN", "YNB")])
)
# p = 0.17. wilcox 0.19


# by type of reporter
ggplot(dat, aes(x = reporter, y = strainDifferenceInReactionToE)) + 
  geom_boxplot(outliers=FALSE) + geom_jitter(aes(color=environment), size=2) + 
#  scale_color_manual(values = wes_palette("AsteroidCity2")) + 
  geom_hline(yintercept = 0, lty = 2) +
  theme_bw() + ylab("Difference in strain response to environment\nBY - RM") + xlab("Reporter")

t.test(
  abs(dat$strainDifferenceInReactionToE[dat$reporter %in% c("Asn", "Phe", "Thr", "UFD")]),
  abs(dat$strainDifferenceInReactionToE[dat$reporter %in% c("4xUb", "Rpn4")])
)
# p = 0.15, p = 0.41
# nope

dat |> heatmap(reporter, environment, strainDifferenceInReactionToE,
               palette_value = c("red", "white", "blue"),
               cluster_rows = FALSE, cluster_columns = FALSE,
               rect_gp = grid::gpar(col = "#161616", lwd = 0.5),
               heatmap_legend_param = list(title = "Strain difference\nin environment effect\n(BY - RM)")
               ) |> 
  layer_diamond(
      interaction_Pval < (0.05 / nrow(dat))
  )
  