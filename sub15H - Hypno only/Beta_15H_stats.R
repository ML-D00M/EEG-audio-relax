### Load some packages that we will need
library(rstatix)
library(ggplot2)
library(ggpubr)
library(R.matlab)
library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(dplyr)


### DATA PREPARATION ----


## set working directory
setwd("C:/Users/User/Desktop/EEG biolabs/audio - Relax/PSD EEG Estimation/C, CP, P, PO, O beta sub16-10/a_Sub15_hypno")

## Load the data saved from matlab

## Hypnosis group (15 subjects)
# 1 condition - Before
beta_15Hsub_1=readMat('beta_15Hsub_cond_1.mat')$beta.sub
beta_15Hsub_1=as.vector(beta_15Hsub_1)
# 2.1 condition - Hypno: Intro
beta_15Hsub_2_1=readMat('beta_15Hsub_cond_2_1.mat')$beta.sub
beta_15Hsub_2_1=as.vector(beta_15Hsub_2_1)
# 2.2 condition - Hypno: Body
beta_15Hsub_2_2=readMat('beta_15Hsub_cond_2_2.mat')$beta.sub
beta_15Hsub_2_2=as.vector(beta_15Hsub_2_2)
# 2.3 condition - Hypno: Mind
beta_15Hsub_2_3=readMat('beta_15Hsub_cond_2_3.mat')$beta.sub
beta_15Hsub_2_3=as.vector(beta_15Hsub_2_3)
# 2.4 condition - Hypno: Negative
beta_15Hsub_2_4=readMat('beta_15Hsub_cond_2_4.mat')$beta.sub
beta_15Hsub_2_4=as.vector(beta_15Hsub_2_4)
# 2.5 condition - Hypno: Positive
beta_15Hsub_2_5=readMat('beta_15Hsub_cond_2_5.mat')$beta.sub
beta_15Hsub_2_5=as.vector(beta_15Hsub_2_5)
# 3. condition - After
beta_15Hsub_3=readMat('beta_15Hsub_cond_3.mat')$beta.sub
beta_15Hsub_3=as.vector(beta_15Hsub_3)


# Create a single data frame from vectors
# to have a similar data frame to https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/#two-way-repeated-measures-anova

id <- sprintf("% d", 1:15)
c1_Before <- c(beta_15Hsub_1)
c2.1_Intro <- c(beta_15Hsub_2_1)
c2.2_Body <- c(beta_15Hsub_2_2)
c2.3_Mind <- c(beta_15Hsub_2_3)
c2.4_Negative <- c(beta_15Hsub_2_4)
c2.5_Positive <- c(beta_15Hsub_2_5)
c3_After <- c(beta_15Hsub_3)

df_beta_15H <- data.frame(id, c1_Before, c2.1_Intro, c2.2_Body, c2.3_Mind, c2.4_Negative, c2.5_Positive, c3_After)

# Gather the columns c1_Before - c3_After into long format.
# Convert id and cond into factor variables
df_beta_15H <- df_beta_15H %>%
  gather(key = "cond", value = "score", c1_Before, c2.1_Intro, c2.2_Body, c2.3_Mind, c2.4_Negative, c2.5_Positive, c3_After) %>%
  convert_as_factor(id, cond)

# Normalize the Score variable by rescaling it.
# Multiply its super small values by 1e13 to have the values >1.
# Thus we can avoid numeric issues due to floating point precision
df_beta_15H$score <- df_beta_15H$score*1e13

# Set color palette for visualizations
# Two contrast colors are chosen using this tool https://paletton.com/#uid=70u1i0ktToJkUvjq-ttxZkPHjfD 
# Hypnosis conditions are shown using different shades of one color
# Before and After conditions are marked with different shades of the other color
cond.colors <- c("c1_Before" ="#EBE128", "c2.1_Intro" = "#F9A057", "c2.2_Body" = "#EB8028", "c2.3_Mind" = "#C5600D", "c2.4_Negative" = "#A64B00", "c2.5_Positive" = "#7D3700", "c3_After" = "#A69E00")

# Summary statistics
df_beta_15H %>%
  group_by(cond) %>%
  get_summary_stats(score, type = "mean_sd")

# Technical Visualization
# Create box plots:
bxp <- ggboxplot(
  df_beta_15H, x = "cond", y = "score",
  add = 'median',
  color = cond.colors
)
bxp



########## ANOVA & PAIRWISE TESTS  --------------------------

##Check Assumptions

# Outliers
df_beta_15H %>%
  group_by(cond) %>%
  identify_outliers(score)
# Normality assumption - Shapiro-Wilk test (p-value > 0.05)
df_beta_15H %>%
  group_by(cond) %>%
  shapiro_test(score)
# QQ plot - correlation between the data and the normal distribution
ggqqplot(df_beta_15H, "score", facet.by = "cond")

# If the assumptions are met, we do One-way Anova & T-tests. 
# Otherwise - Kruskal-Wallis & Wilcoxon tests


## 1.1. Kruskal-Wallis & Wilcoxon test - paired ####

# Kruskal-Wallis test
one.way.np <- df_beta_15H %>%
  kruskal_test(score ~ cond)
one.way.np

# Effect size of the Cond on Beta score in Group 1
df_beta_15H %>%
  kruskal_effsize(score ~ cond)

# Wilcoxon test
Wlx.test <- df_beta_15H %>% 
  wilcox_test(score ~ cond, paired = TRUE,
              p.adjust.method = "BH") %>%
  add_significance()
Wlx.test

# Effect size
df_beta_15H %>% wilcox_effsize(score ~ cond)

# Report

# Bar-plot with p-values for 1 grouping variable - Cond
Wlx.test = Wlx.test %>% add_xy_position(x = "cond")
ggbarplot(df_beta_15H, x = "cond", y = "score", 
          ylab='median Beta score', fill=cond.colors,
          add = c("median"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test, label = "p = {p}", 
                     y.position = 550,
                     tip.length = 0.03, size=3,
                     hide.ns = TRUE,
                     step.increase = 0.5)+
  coord_cartesian(ylim = c(0, 700))+ # to set a limit to y value
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?????? ???????????? ???????????? - ????, ?????????? ?? ???? ?????????? 5 ???????????????????? ??????????-????????????",
    caption = paste("Kruskal-Wallis test: ", "p = ", round(one.way.np$p, digits = 3),
                    ". pwc: Wilcoxon test"))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 1.1. One-way Anova & T-test - paired ####

# One-Way Anova - simple main effect of Cond variable
one.way <- df_beta_15H %>%
  anova_test(dv = score, wid = id,  within = cond) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

# Compute Weltch T-test, which is the safer one
ttest <- df_beta_15H %>% 
  t_test(score ~ cond, paired = TRUE) %>%
  add_significance()
ttest

# Effect size
df_beta_15H %>% cohens_d(score ~ cond, paired = TRUE, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest = ttest %>% add_xy_position(x = "cond")

# Bar-plot for 1 grouping variable - Cond
ttest = ttest %>% add_xy_position(x = "cond")
ggbarplot(df_beta_15H, x = "cond", y = "score", 
          ylab='mean Beta score', fill=cond.colors,
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest, label = "p = {p}", 
                     y.position = 20,
                     tip.length = 0.03, size=3,
                     #hide.ns = TRUE,
                     step.increase = 0.1)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?????? ???????????? ???????????? - ????, ?????????? ?? ???? ?????????? 5 ???????????????????? ??????????-????????????",
    caption = get_test_label(one.way, detailed= FALSE))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))