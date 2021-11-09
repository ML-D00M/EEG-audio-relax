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
setwd("C:/Users/User/Desktop/EEG biolabs/audio - Relax/PSD EEG Estimation/C, CP, P, PO, O alpha sub16-10")

## Load the data saved from matlab

## 1 group - Hypnosis (16 subjects)
# 1st condition
beta_16sub_1=readMat('beta_16sub_cond_1.mat')$beta.sub
beta_16sub_1=as.vector(beta_16sub_1)
# 2nd condition
beta_16sub_2=readMat('beta_16sub_cond_2.mat')$beta.sub
beta_16sub_2=as.vector(beta_16sub_2)
# 3rd condition
beta_16sub_3=readMat('beta_16sub_cond_3.mat')$beta.sub
beta_16sub_3=as.vector(beta_16sub_3)

## 2 group - Control/Sea (10 subjects)
# 1st condition
beta_10Csub_1=readMat('beta_10Csub_cond_1.mat')$beta.sub 
beta_10Csub_1=as.vector(beta_10Csub_1)
# 2nd condition
beta_10Csub_2=readMat('beta_10Csub_cond_2.mat')$beta.sub
beta_10Csub_2=as.vector(beta_10Csub_2)
# 3rd condition
beta_10Csub_3=readMat('beta_10Csub_cond_3.mat')$beta.sub
beta_10Csub_3=as.vector(beta_10Csub_3)


# Create a single data frame from vectors
# to have a similar data frame to https://www.datanovia.com/en/lessons/mixed-anova-in-r/

id <- sprintf("% d", 1:26)
group <- c("1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis",
           "2.Control", "2.Control", "2.Control", "2.Control", "2.Control", "2.Control", "2.Control", "2.Control","2.Control", "2.Control")
c1_Before <- c(beta_16sub_1, beta_10Csub_1)
c2_During <- c(beta_16sub_2, beta_10Csub_2)
c3_After <- c(beta_16sub_3, beta_10Csub_3)

df_beta_16_10 <- data.frame(id, group, c1_Before, c2_During, c3_After)

# Gather the columns c1_Before and c2_After into long format.
# Convert id and cond into factor variables
df_beta_16_10 <- df_beta_16_10 %>%
  gather(key = "cond", value = "score", c1_Before, c2_During, c3_After) %>%
  convert_as_factor(id, cond)

# Normalize the Score variable by rescaling it.
# Multiply its super small values by 1e13 to have the values >1.
# Thus we can avoid numeric issues due to floating point precision
df_beta_16_10$score <- df_beta_16_10$score*1e13

# Inspect some random rows of the data by groups
set.seed(123)
df_beta_16_10 %>% sample_n_by(group, cond, size = 1)

# Summary statistics
# Group the data by cond and group, and then compute some summary statistics 
# of the score variable: mean and sd (standard deviation)
df_beta_16_10 %>%
  group_by(cond, group) %>%
  get_summary_stats(score, type = "mean_sd") 
# This function is hardcoded to return only 3 digits, so if we don't
# normalize the Score variable, it won't work. In this case
# we check the stats like this
#sum.table <- summary(df_beta_16_10)
#as.numeric(sub('.*:', '', sum.table[,4]))

# Set color palette for visualizations
# Two contrast colors are chosen using this tool https://paletton.com/#uid=20u1o0ktToJkUvjqzttxYkPHifD 
# The main distinction is between Groups. 
# Different Conditions in the same Group can be shown using different
# shades of a group color
group.colors <- c("1.Hypnosis" = "#C5610D", "2.Control" = "#087676")
cond.colors <- c("c1_Before" ="#C5610D", "c2_During" = "#C5BC0D", "c3_After" = "#087676")

# Technical Visualization
# Create box plots:
bxp <- ggboxplot(
  df_beta_16_10, x = "cond", y = "score",
  color = "group", palette = "jco"
)
bxp



### Two-way Mixed ANOVA analysis: within-subjects factor - condition, between-subjects factor - group ---------


### CHECK ASSUMPTIONS

## 1. Outliers
# Outliers can be easily identified using box plot methods, 
# implemented in the R function identify_outliers() [rstatix package]
out=df_beta_16_10 %>%
  group_by(cond, group) %>%
  identify_outliers(score)
out
# Remove just one outlier (id = 8 & 18) from Groups 1 & 2
#df_beta_16_10 <- df_beta_16_10[-c(8,18,34,44), ]
# to Remove all the outliers
#tbr=which(df_beta_16_10$score%in%out$score)
#df_beta_16_10=df_beta_16_10[-tbr,]

df_beta_16_10 %>%
  group_by(cond, group) %>%
  identify_outliers(score)

## 2. Normality
# This assumption can be checked by computing Shapiro-Wilk test 
# for each combinations of factor levels. 
# If the data is normally distributed, the p-value should be greater than 0.05
df_beta_16_10 %>%
  group_by(cond, group) %>%
  shapiro_test(score)

# QQ plot draws the correlation between a given data and the normal distribution
ggqqplot(df_beta_16_10, "score", ggtheme = theme_bw()) +
  facet_grid(cond ~ group)

## 3. The homogeneity of variance 
# assumption of the between-subject factor (group) can be checked using the Levene's test.
# The test is performed at each level of cond variable
df_beta_16_10 %>%
  group_by(cond) %>%
  levene_test(score ~ group)


### COMPUTATION

## 1. Two-way mixed ANOVA test
res.aov <- anova_test(
  data = df_beta_16_10, dv = score, wid = id,
  between = group, within = cond
)
get_anova_table(res.aov)


## 2. Nonparametric ANOVA - insensitive to outliers
library(nparLD)
np_anova_model<-nparLD(score ~ cond * group, data = df_beta_16_10,
                       subject = "id", description = FALSE)
summary(np_anova_model)

# Check the significance of the two-way interaction !



### POST-HOC TESTS -----


## 1. Procedure for SIGNIFICANT two-way interaction ====


## 1.1. T-tests & One-Way Anova - for parametric (regular) ANOVA ====


## 1.1.1. Paired T-test & One-Way Anova ####

# We investigate the effect of the within-subject factor "Cond"
# on the Beta score in each Group

# One-Way Anova - simple main effect of Cond variable
one.way2 <- df_beta_16_10 %>%
  group_by(group) %>%
  anova_test(dv = score, wid = id, within = cond) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

## Paired T-test - Cond 1 vs Cond 2 vs Cond 3 in Group 1 & 3
ttest.p <- df_beta_16_10 %>%
  group_by(group) %>%
  pairwise_t_test(
    score ~ cond, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )
ttest.p

# Effect size
df_beta_16_10 %>%  
  cohens_d(score ~ cond, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.p = ttest.p %>% add_xy_position(x = "cond")

# Bar-plot for 2 grouping variables
ggbarplot(df_beta_16_10, x = "cond", y = "score", 
          ylab='mean Beta score', fill='group',
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=group.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.p, label = "p = {p}", 
                     y.position = 15,
                     tip.length = 0.03, size=3,
                     step.increase = 0.1)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? ?????? ?????????? ??????????: ???????????? ?? ????????????????",
    caption = get_test_label(res.aov, detailed= FALSE))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 1.1.2. Non-paired T-test & One-Way Anova ####

# We investigate the effect of the between-subject factor "Group"
# on the Beta score in each ??ondition

# One-Way Anova - simple main effect of the Group variable
one.way <- df_beta_16_10 %>%
  group_by(cond) %>%
  anova_test(dv = score, wid = id, between = group) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way

## Non-paired T-test - Group 1 vs Group 2 in Cond 1, 2 & 3
ttest <- df_beta_16_10 %>%
  group_by(cond) %>%
  pairwise_t_test(
    score ~ group, paired = FALSE, 
    p.adjust.method = "bonferroni"
  )
ttest

# Effect size
df_beta_16_10 %>%  
  cohens_d(score ~ group, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest = ttest %>% add_xy_position(x = "group")

# Bar-plot for 2 grouping variables
ggbarplot(df_beta_16_10, x = "group", y = "score", 
          ylab='mean Beta score', fill='cond',
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest, label = "p = {p}", 
                     y.position = 15,
                     tip.length = 0.03, size=3,
                     step.increase = 0.1)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? ?????? ???????? ??????????????: ????, ???? ?????????? ?? ?????????? ??????????-????????????",
    caption = get_test_label(res.aov, detailed= FALSE))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 1.2. Wilcoxon & Kruskal-Wallis tests - for a non-parametric ANOVA ====

# Wilcoxon test is a non-parametric alternative to the T-test for comparing two means.
# Kruskal-Wallis is a non-parametric alternative to the one-way Anova.
# It extends the two-samples Wilcoxon test in the situation where there are more than two groups to compare.
# It's particularly recommended in a situation where the data are not normally distributed.


## 1.2.1. Paired Wilcoxon & Kruskal-Wallis tests ####

# We investigate the effect of the within-subject factor "Cond"
# on the Beta score in each Group

# Kruskal-Wallis test 1
one.way.np <- df_beta_16_10 %>%
  group_by(group) %>%
  kruskal_test(score ~ cond)
one.way.np

# Effect size of the Cond on Beta score
df_beta_16_10 %>%
  group_by(group) %>%
  kruskal_effsize(score ~ cond)

# Paired Wilcoxon test - Cond 1 vs Cond 2 vs Cond 3 in Group 1 & 2
Wlx.test.p <- df_beta_16_10 %>%
  group_by(group) %>%
  pairwise_wilcox_test(
    score ~ cond, paired = TRUE,
    p.adjust.method = "BH"
  )
Wlx.test.p

#Effect size
df_beta_16_10  %>%
  group_by(group) %>%
  wilcox_effsize(score ~ cond, paired = TRUE)

# Report
# Auto-compute p-value label positions
Wlx.test.p = Wlx.test.p %>% add_xy_position(x = "cond")
# Bar-plot for 2 grouping variables
ggbarplot(df_beta_16_10, x = "cond", y = "score", 
          ylab='median Beta score', fill='group',
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  scale_fill_manual(values=group.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.p, label = "p = {p}", 
                     y.position = 650,
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? ?????? ?????????? ??????????: ???????????? ?? ????????????????",
    caption = paste("Tests: non-parametric Anova, ", "p = ",  round(np_anova_model$ANOVA.test[3,3], digits = 3), 
                    ". One-way - Kruskal-Wallis test: ", "1. Hypno p = ", round(one.way.np$p[1], digits = 3),
                    ". 2. Ctrl p = ", round(one.way.np$p[2], digits = 3),
                    ". pwc: Wilcoxon test"))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 1.2.2. Non-paired Wilcoxon & Kruskal-Wallis tests ####

# We investigate the effect of the between-subject factor "Group"
# on the Beta score in each condition

# Kruskal-Wallis test 2
one.way.np2 <- df_beta_16_10 %>%
  group_by(cond) %>%
  kruskal.test(score ~ group, paired = FALSE)
one.way.np2

# Effect size of the Group on Beta score
df_beta_16_10 %>%
  #group_by(cond) %>%
  kruskal_effsize(score ~ group)

# Non-paired Wilcoxon test - Group 1 vs Group 2 in Cond 1, 2 & 3
Wlx.test <- df_beta_16_10 %>%
  group_by(cond) %>%
  wilcox_test(
    score ~ group, paired = FALSE,
    p.adjust.method = "BH"
  )
Wlx.test

#Effect size
df_beta_16_10  %>%
  wilcox_effsize(score ~ group, paired = FALSE)

# Report
# Auto-compute p-value label positions
Wlx.test = Wlx.test %>% add_xy_position(x = "group")
# Bar-plot for 2 grouping variables
ggbarplot(df_beta_16_10, x = "group", y = "score", 
          ylab='median Beta score', fill='cond',
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test, label = "p = {p}", 
                     y.position = 650,
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? ?????? ???????? ??????????????: ????, ???? ?????????? ?? ?????????? ??????????-????????????",
    caption = paste("Tests: non-parametric Anova, ", "p = ",  round(np_anova_model$ANOVA.test[3,3], digits = 3), 
                    ". One-way - Kruskal-Wallis test, ", "p = ", formatC(one.way.np2[["p.value"]], format = "e", digits = 2),
                    ". pwc: Wilcoxon test"))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


Wlx.test.p = Wlx.test.p %>% add_xy_position(x = "group")
ggbarplot(df_beta_16_10, x = "group", y = "score", 
          ylab='median Beta score', fill='cond',
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.p, label = "p = {p}", 
                     y.position = 450,
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  #coord_cartesian(ylim = c(0, 7.5))+ # to set a limit to y value
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? ?????? ???????? ??????????????: ????, ???? ?????????? ?? ?????????? ??????????-????????????",
    caption = paste("Tests: non-parametric Anova, ", "p = ",  round(np_anova_model$ANOVA.test[3,3], digits = 3), 
                    ". One-way - Kruskal-Wallis test: ", "1. Hypno p = ", round(one.way.np$p[1], digits = 3),
                    ". 2. Ctrl p = ", round(one.way.np$p[2], digits = 3),
                    ". pwc: Wilcoxon test"))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 2. Procedure for NON-SIGNIFICANT two-way interaction ====

# If the interaction is not significant, you need to interpret 
# the main effects for each of the two variables: Group and Cond.
# A significant main effect can be followed up with pairwise comparisons.


## 2.1. T-tests - for parametric (regular) ANOVA ====

## 2.1.1. Paired T-test - Cond 1 vs Cond 2 vs Cond 3 in Group 1 & 2 ####

# Perform multiple pairwise paired t-tests for the Cond variable, 
# ignoring Group. P-values are adjusted using the Bonferroni multiple testing correction method

ttest.p <- df_beta_16_10 %>%
  # here we don't group by Group
  pairwise_t_test(
    score ~ cond, paired = TRUE, 
    p.adjust.method = "bonferroni"
  )
ttest.p

# Effect size
df_beta_16_10 %>%  
  cohens_d(score ~ cond, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.p = ttest.p %>% add_xy_position(x = "cond")

# Bar-plot for 2 grouping variables
ggbarplot(df_beta_16_10, x = "cond", y = "score", 
          ylab='mean Beta score', fill='group',
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=group.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.p, label = "p = {p.adj}", 
                     y.position = 45,
                     tip.length = 0.03, size=3,
                     step.increase = 0.1)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ???????????????? ?????? ?????????? ???????????? ?? ????????????????",
    caption = "pwc: T-test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 2.1.2. Non-paired T-test - Group 1 vs Group 2 in Cond 1, 2 & 3 ####

# Perform multiple pairwise paired t-tests for the Group variable, 
# ignoring Cond. P-values are adjusted using the Bonferroni multiple testing correction method

ttest <- df_beta_16_10 %>%
  # here we don't group by cond
  pairwise_t_test(
    score ~ group, paired = FALSE, 
    p.adjust.method = "bonferroni"
  )
ttest

# Effect size
df_beta_16_10 %>%  
  cohens_d(score ~ group, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest = ttest %>% add_xy_position(x = "group")

# Bar-plot for 2 grouping variables
ggbarplot(df_beta_16_10, x = "group", y = "score", 
          ylab='mean Beta score', fill='cond',
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest, label = "p = {p.adj}", 
                     y.position = 15,
                     tip.length = 0.03, size=3,
                     step.increase = 0.1)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ???????????????? ?????? ???????? ??????????????: ????, ???? ?????????? ?? ?????????? ??????????-????????????",
    caption = "pwc: T-test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 2.2. Wilcoxon tests - for a non-parametric ANOVA ====

# It's a non-parametric alternative to the T-test for comparing two means.
# It's particularly recommended in a situation where the data are not normally distributed.

## 2.2.1. Paired Wilcoxon test - Cond 1 vs Cond 2 vs Cond 3 in Group 1 & 2 ####

Wlx.test.p <- df_beta_16_10 %>%
  # here we don't group by Group
  pairwise_wilcox_test(
    score ~ cond, paired = TRUE,
    p.adjust.method = "BH"
  )
Wlx.test.p

#Effect size
df_beta_16_10  %>%
  wilcox_effsize(score ~ cond, paired = TRUE)

# Report
# Auto-compute p-value label positions
Wlx.test.p = Wlx.test.p %>% add_xy_position(x = "cond")
# Bar-plot for 1 grouping variable - Cond
ggbarplot(df_beta_16_10, x = "cond", y = "score", 
          ylab='median Beta score', 
          fill=cond.colors,
          add = c("median"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=group.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.p, label = "p = {p.adj}", 
                     y.position = 5,
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ???????????????? ?????? ?????????? ???????????? ?? ????????????????",
    caption = "pwc: Wilcoxon test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 2.2.1. Non-paired Wilcoxon test - Group 1 vs Group 2 in Cond 1, 2 & 3 ####

Wlx.test <- df_beta_16_10 %>%
  # here we don't group by Cond
  wilcox_test(
    score ~ group, paired = FALSE,
    p.adjust.method = "BH"
  )
Wlx.test

#Effect size
df_beta_16_10  %>%
  wilcox_effsize(score ~ group, paired = FALSE)

# Report
# Auto-compute p-value label positions
Wlx.test = Wlx.test %>% add_xy_position(x = "group")
# Bar-plot for 1 grouping variable - Group
ggbarplot(df_beta_16_10, x = "group", y = "score", 
          ylab='median Beta score', 
          fill=group.colors,
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  geom_jitter(size=1, width=0.1)+ # add subjects
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test, label = "p = {p}", 
                     y.position = 5000, # put 5000 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ???????????????? ?????? ???????? ??????????????: ????, ???? ?????????? ?? ?????????? ??????????-????????????",
    caption = "pwc: Wilcoxon test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))




#________________________________________________________




########## PAIRWISE TESTS (NOT POST-HOC) --------------------------


# Some additional data for parwise tests
idg1 <- sprintf("% d", 1:16)
idg2 <- sprintf("% d", 1:10)
g1_Hypnosis <- c("1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis", "1.Hypnosis")
g2_Sea <- c("2.Control", "2.Control", "2.Control", "2.Control", "2.Control", "2.Control", "2.Control", "2.Control","2.Control", "2.Control")
c1_Before_Hypno <- c(beta_16sub_1)
c2_During_Hypno <- c(beta_16sub_2)
c3_After_Hypno <- c(beta_16sub_3)
c1_Before_Ctrl <- c(beta_10Csub_1)
c2_During_Ctrl <- c(beta_10Csub_2)
c3_After_Ctrl <- c(beta_10Csub_3)

# Additional color shades for visualizations
cond.colors_hyp <- c(c1_Before_Hypno ="#EB8028", c2_During_Hypno = "#C5600D", c3_After_Hypno = "#A64B00")
cond.colors_ctrl <- c(c1_Before_Ctrl ="#188D8D", c2_During_Ctrl = "#087676", c3_After_Ctrl = "#006464")


## 1. Cond 1 vs Cond 2 vs Cond 3 for Group 1 - paired comparison =====

# Data preparation
df_beta_16_10_g1 <- data.frame(idg1, c1_Before_Hypno, c2_During_Hypno, c3_After_Hypno) # create df
df_beta_16_10_g1 <- df_beta_16_10_g1 %>% # convert df
  gather(key = "cond", value = "score", c1_Before_Hypno, c2_During_Hypno, c3_After_Hypno) %>%
  convert_as_factor(idg1, cond)
#df_beta_16_10_g1 <- df_beta_16_10_g1[-c(18), ] # remove outlier-18
df_beta_16_10_g1$score <- df_beta_16_10_g1$score*1e13 # normalization

# Summary stats
df_beta_16_10_g1 %>%
  group_by(cond) %>%
  get_summary_stats(score, type = "mean_sd")

## Check assumptions

# Check outliers 
df_beta_16_10_g1 %>%
  group_by(cond) %>%
  identify_outliers(score) # should be no extreme outliers
# Check normality - data of 2 groups should be normally distributed
# Shapiro Wilk test
df_beta_16_10_g1 %>%
  group_by(cond) %>%
  shapiro_test(score)
# Draw a qq plot by cond
ggqqplot(df_beta_16_10_g1, x = "score", facet.by = "cond")
# Check the equality of variances
df_beta_16_10_g1 %>% levene_test(score ~ cond) # p should be >0.05

# If the assumptions are met, we do T-test. Otherwise - Wilcoxon test


## 1.1. One-way Anova & T-test - paired ####

# One-Way Anova - simple main effect of Cond variable
one.way.g1 <- df_beta_16_10_g1 %>%
  anova_test(dv = score, wid = idg1,  within = cond) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way.g1

# Compute Weltch T-test, which is the safer one
ttest.g1 <- df_beta_16_10_g1 %>% 
  t_test(score ~ cond, paired = TRUE) %>%
  add_significance()
ttest.g1

# Effect size
df_beta_16_10_g1 %>% cohens_d(score ~ cond, paired = TRUE, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.g1 = ttest.g1 %>% add_xy_position(x = "cond")

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_g1, x = "cond", y = "score", 
          ylab='mean Beta score', fill=cond.colors_hyp,
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=group.colors_hyp)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.g1, label = "p = {p}", 
                     y.position = 15,
                     tip.length = 0.03, size=3,
                     step.increase = 0.1, paired = TRUE)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????",
    caption = get_test_label(one.way.g1, detailed= FALSE))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))

# Visualization for paired samples
ggpaired(df_beta_16_10_g1, x = "cond", y = "score", 
         fill=cond.colors_hyp,
         point.size=1.5,line.color = "gray", line.size = 0.4)+
  ylab('Beta score')+
  stat_compare_means(paired = TRUE)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????",
    caption = "???????????? ?????????? ???? ?????????????? - ?????? ?????????????????? ????????????????????. 
    ?????????? ?????????????????? ???????????????????? ?????????? ?? ?????? ???? ???????????????????? ???? ?? ?????????? ??????????-????????????")+
  theme(text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 1.2. Kruskal-Wallis & Wilcoxon test - paired ####

# Kruskal-Wallis test 1
one.way.np.g1 <- df_beta_16_10_g1 %>%
  kruskal_test(score ~ cond)
one.way.np.g1

# Effect size of the Cond on Beta score in Group 1
df_beta_16_10_g1 %>%
  kruskal_effsize(score ~ cond)

# Wilcoxon test
Wlx.test.g1 <- df_beta_16_10_g1 %>% 
  wilcox_test(score ~ cond, paired = TRUE,
              p.adjust.method = "BH") %>%
  add_significance()
Wlx.test.g1

# Effect size
df_beta_16_10_g1 %>% wilcox_effsize(score ~ cond)

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_g1, x = "cond", y = "score", 
          ylab='median Beta score', fill=cond.colors_hyp,
          add = c("median"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  scale_fill_manual(values=group.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.g1, label = "p = {p}", 
                     y.position = 650, # put 1500 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????",
    caption = paste("Kruskal-Wallis test: ", "p = ", round(one.way.np.g1$p, digits = 3),
                    ". pwc: Wilcoxon test"))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))

# Visualization for paired samples
ggpaired(df_beta_16_10_g1, x = "cond", y = "score", 
         fill=cond.colors_hyp,
         point.size=1.5,line.color = "gray", line.size = 0.4)+
  ylab('Beta score')+
  stat_compare_means(paired = TRUE)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????",
    caption = "???????????? ?????????? ???? ?????????????? - ?????? ?????????????????? ????????????????????. 
    ?????????? ?????????????????? ???????????????????? ?????????? ?? ?????? ???? ???????????????????? ????, ???? ?????????? ?? ?????????? ??????????-????????????")+
  theme(text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 2. Cond 1 vs Cond 2 vs Cond 3 for Group 2 - paired ====

# Data preparation
df_beta_16_10_g2 <- data.frame(idg2, c1_Before_Ctrl, c2_During_Ctrl, c3_After_Ctrl) # create df
df_beta_16_10_g2 <- df_beta_16_10_g2 %>% # convert df
  gather(key = "cond", value = "score", c1_Before_Ctrl, c2_During_Ctrl, c3_After_Ctrl) %>%
  convert_as_factor(idg2, cond)
#df_beta_16_10_g1 <- df_beta_16_10_g1[-c(18), ] # remove outlier-18
df_beta_16_10_g2$score <- df_beta_16_10_g2$score*1e13 # normalization

# Summary stats
df_beta_16_10_g2 %>%
  group_by(cond) %>%
  get_summary_stats(score, type = "mean_sd")

## Check assumptions

# Check outliers 
df_beta_16_10_g2 %>%
  group_by(cond) %>%
  identify_outliers(score) # should be no extreme outliers
# Check normality - data of 2 groups should be normally distributed
# Shapiro Wilk test
df_beta_16_10_g2 %>%
  group_by(cond) %>%
  shapiro_test(score)
# Draw a qq plot by cond
ggqqplot(df_beta_16_10_g2, x = "score", facet.by = "cond")
# Check the equality of variances
df_beta_16_10_g2 %>% levene_test(score ~ cond) # p should be >0.05

# If the assumptions are met, we do T-test. Otherwise - Wilcoxon test

## 2.1. One-way Anova & T-test - paired ####

# One-Way Anova - simple main effect of Cond variable for Group 2
one.way.g2 <- df_beta_16_10_g2 %>%
  anova_test(dv = score, wid = idg2,  within = cond) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way.g2

# Compute Weltch T-test, which is the safer one
ttest.g2 <- df_beta_16_10_g2 %>% 
  t_test(score ~ cond, paired = TRUE) %>%
  add_significance()
ttest.g2

# Effect size
df_beta_16_10_g2 %>% cohens_d(score ~ cond, paired = TRUE, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.g2 = ttest.g2 %>% add_xy_position(x = "cond")
# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_g2, x = "cond", y = "score", 
          ylab='mean Beta score', fill=cond.colors_ctrl,
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=group.colors_ctrl)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.g2, label = "p = {p}", 
                     y.position = 25, #put 25 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.1, paired = TRUE)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????????",
    caption = get_test_label(one.way.g2, detailed= FALSE))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))

# Visualization for paired samples
ggpaired(df_beta_16_10_g2, x = "cond", y = "score", 
         fill=cond.colors_ctrl,
         point.size=1.5,line.color = "gray", line.size = 0.4)+
  ylab('Beta score')+
  stat_compare_means(paired = TRUE)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????????",
    caption = "???????????? ?????????? ???? ?????????????? - ?????? ?????????????????? ????????????????????. 
    ?????????? ?????????????????? ???????????????????? ?????????? ?? ?????? ???? ???????????????????? ????, ???? ?????????? ?? ?????????? ??????????-????????????")+
  theme(text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 2.2. Kruskal-Wallis & Wilcoxon test - paired ####

# Kruskal-Wallis test 2
one.way.np.g2 <- df_beta_16_10_g2 %>%
  kruskal_test(score ~ cond)
one.way.np.g2

# Effect size of the Cond on Beta score in Group 2
df_beta_16_10_g2 %>%
  kruskal_effsize(score ~ cond)

# Wilcoxon test
Wlx.test.g2 <- df_beta_16_10_g2 %>% 
  wilcox_test(score ~ cond, paired = TRUE,
              p.adjust.method = "BH") %>%
  add_significance()
Wlx.test.g2

# Effect size
df_beta_16_10_g2 %>% wilcox_effsize(score ~ cond)

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_g2, x = "cond", y = "score", 
          ylab='median Beta score', fill=cond.colors_ctrl,
          add = c("median"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=group.colors_ctrl)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.g2, label = "p = {p}", 
                     y.position = 650, # put 4000 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????????",
    caption = paste("Kruskal-Wallis test: ", "p = ", round(one.way.np.g2$p, digits = 3),
                    ". pwc: Wilcoxon test"))+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))

# Visualization for paired samples
ggpaired(df_beta_16_10_g2, x = "cond", y = "score", 
         fill=cond.colors_ctrl,
         point.size=1.5,line.color = "gray", line.size = 0.4)+
  ylab('Beta score')+
  stat_compare_means(paired = TRUE)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "????, ???? ?????????? ?? ?????????? ??????????-???????????? - ?????? ???????????? ????????????????",
    caption = "???????????? ?????????? ???? ?????????????? - ?????? ?????????????????? ????????????????????. 
    ?????????? ?????????????????? ???????????????????? ?????????? ?? ?????? ???? ???????????????????? ????, ???? ?????????? ?? ?????????? ??????????-????????????")+
  theme(text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))



### 3. Group 1 vs Group 2 for Cond 1 - non-paired ====

## Data preparation (the same as in the beginning)

df_beta_16_10_c1 <- data.frame(id, group, c1_Before) # create df
names(df_beta_16_10_c1)[names(df_beta_16_10_c1) == 
                                          'c1_Before'] <- 'score'
#df_beta_16_10_c1 <- df_beta_16_10_c1[-c(18), ] # remove outlier-18
df_beta_16_10_c1$score <- df_beta_16_10_c1$score*1e13 # normalization

# Summary stats
df_beta_16_10_c1 %>%
  group_by(group) %>%
  get_summary_stats(score, type = "mean_sd")

# Check assumptions

# Check outliers 
df_beta_16_10_c1 %>%
  group_by(group) %>%
  identify_outliers(score) # should be no extreme outliers
# Check normality - data of 2 groups should be normally distributed
# Shapiro Wilk test
df_beta_16_10_c1 %>%
  group_by(group) %>%
  shapiro_test(score)
# Draw a qq plot by group
ggqqplot(df_beta_16_10_c1, x = "score", facet.by = "group")
# Check the equality of variances
df_beta_16_10_c1 %>% levene_test(score ~ group) # p should be >0.05

# If the assumptions are met, we do T-test. Otherwise - Wilcoxon test


## 3.1. T-test - non-paired ####

# Compute Weltch T-test, which is the safer one
ttest.c1 <- df_beta_16_10_c1 %>% 
  t_test(score ~ group, paired = FALSE) %>%
  add_significance()
ttest.c1

# Effect size
df_beta_16_10_c1 %>% cohens_d(score ~ group, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.c1 = ttest.c1 %>% add_xy_position(x = "group")

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_c1, x = "group", y = "score", 
          ylab='mean Beta score',
          fill = group.colors,
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.c1, label = "p = {p}", 
                     y.position = 15,
                     tip.length = 0.03, size=3,
                     step.increase = 0)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ???? ??????????-????????????",
    caption = "pwc: T-test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 3.2. Wilcoxon test - non-paired ####

Wlx.test.c1 <- df_beta_16_10_c1 %>% 
  wilcox_test(score ~ group, paired = FALSE,
              p.adjust.method = "BH") %>%
  add_significance()
Wlx.test.c1

# Effect size
df_beta_16_10_c1 %>% wilcox_effsize(score ~ group)
ggbarplot(df_beta_16_10_c1, x = "group", y = "score", 
          ylab='median Beta score', fill=group.colors,
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.c1, label = "p = {p}", 
                     y.position = 650, # put 5000 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ???? ??????????-????????????",
    caption = "pwc: Wilcoxon test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))
# Bar-plot for 1 grouping variable




## 4. Group 1 vs Group 2 for Cond 2 - non-paired ====

## Data preparation

df_beta_16_10_c2 <- data.frame(id, group, c2_During) # create df
names(df_beta_16_10_c2)[names(df_beta_16_10_c2) == 
                                          'c2_During'] <- 'score'
#df_beta_16_10_c2 <- df_beta_16_10_c2[-c(18), ] # remove outlier-18
df_beta_16_10_c2$score <- df_beta_16_10_c2$score*1e13 # normalization

# Summary stats
df_beta_16_10_c2 %>%
  group_by(group) %>%
  get_summary_stats(score, type = "mean_sd")

## Check assumptions

# Check outliers 
df_beta_16_10_c2 %>%
  group_by(group) %>%
  identify_outliers(score) # should be no extreme outliers
# Check normality - data of 2 groups should be normally distributed
# Shapiro Wilk test
df_beta_16_10_c2 %>%
  group_by(group) %>%
  shapiro_test(score)
# Draw a qq plot by group
ggqqplot(df_beta_16_10_c2, x = "score", facet.by = "group")
# Check the equality of variances
df_beta_16_10_c2 %>% levene_test(score ~ group) # p should be >0.05

# If the assumptions are met, we do T-test. Otherwise - Wilcoxon test


## 4.1. T-test - non-paired ####

# Compute Weltch T-test, which is the safer one
ttest.c2 <- df_beta_16_10_c2 %>% 
  t_test(score ~ group, paired = FALSE) %>%
  add_significance()
ttest.c2

# Effect size
df_beta_16_10_c2 %>% cohens_d(score ~ group, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.c2 = ttest.c2 %>% add_xy_position(x = "group")
# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_c2, x = "group", y = "score", 
          ylab='mean Beta score',
          fill = group.colors,
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.c2, label = "p = {p}", 
                     y.position = 15, # set 15 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ???? ?????????? ??????????-????????????",
    caption = "pwc: T-test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 4.2. Wilcoxon test - non-paired #### 

Wlx.test.c2 <- df_beta_16_10_c2 %>% 
  wilcox_test(score ~ group, paired = FALSE,
              p.adjust.method = "BH") %>%
  add_significance()
Wlx.test.c2

# Effect size
df_beta_16_10_c2 %>% wilcox_effsize(score ~ group)

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_c2, x = "group", y = "score", 
          ylab='median Beta score', fill=group.colors,
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.c2, label = "p = {p}", 
                     y.position = 650, # set 2000 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ???? ?????????? ??????????-????????????",
    caption = "pwc: Wilcoxon test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 5. Group 1 vs Group 2 for Cond 3 - non-paired ==== 

## Data preparation

df_beta_16_10_c3 <- data.frame(id, group, c3_After) # create df
names(df_beta_16_10_c3)[names(df_beta_16_10_c3) == 
                                          'c3_After'] <- 'score'
#df_beta_16_10_c3 <- df_beta_16_10_c3[-c(18), ] # remove outlier-18
df_beta_16_10_c3$score <- df_beta_16_10_c3$score*1e13 # normalization

# Summary stats
df_beta_16_10_c3 %>%
  group_by(group) %>%
  get_summary_stats(score, type = "mean_sd")

## Check assumptions

# Check outliers 
df_beta_16_10_c3 %>%
  group_by(group) %>%
  identify_outliers(score) # should be no extreme outliers
# Check normality - data of 2 groups should be normally distributed
# Shapiro Wilk test
df_beta_16_10_c3 %>%
  group_by(group) %>%
  shapiro_test(score)
# Draw a qq plot by group
ggqqplot(df_beta_16_10_c3, x = "score", facet.by = "group")
# Check the equality of variances
df_beta_16_10_c3 %>% levene_test(score ~ group) # p should be >0.05

# If the assumptions are met, we do T-test. Otherwise - Wilcoxon test


## 5.1. T-test - non-paired ####

# Compute Weltch T-test, which is the safer one
ttest.c3 <- df_beta_16_10_c3 %>% 
  t_test(score ~ group, paired = FALSE) %>%
  add_significance()
ttest.c3

# Effect size
df_beta_16_10_c3 %>% cohens_d(score ~ group, var.equal = FALSE)

# Report
# Auto-compute p-value label positions
ttest.c3 = ttest.c3 %>% add_xy_position(x = "group")

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_c3, x = "group", y = "score", 
          ylab='mean Beta score',
          fill = group.colors,
          add = c("mean_se"),
          error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  #geom_jitter(size=1, width=0.1)+ # add subjects
  #scale_fill_manual(values=cond.colors)+ #set colors
  # add pairwise comparisons p-values
  stat_pvalue_manual(ttest.c3, label = "p = {p}", 
                     y.position = 15,
                     tip.length = 0.03, size=3,
                     step.increase = 0)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ?????????? ??????????-????????????",
    caption = "pwc: T-test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))


## 5.2. Wilcoxon test - non-paired #### 

Wlx.test.c3 <- df_beta_16_10_c3 %>% 
  wilcox_test(score ~ group, paired = FALSE,
              p.adjust.method = "BH") %>%
  add_significance()
Wlx.test.c3

# Effect size
df_beta_16_10_c3 %>% wilcox_effsize(score ~ group)

# Bar-plot for 1 grouping variable
ggbarplot(df_beta_16_10_c3, x = "group", y = "score", 
          ylab='median Beta score', fill=group.colors,
          add = c("median"),
          #error.plot = 'upper_errorbar',
          position = position_dodge(width=0.8))+
  geom_jitter(size=1, width=0.1)+ # add subjects
  # add pairwise comparisons p-values
  stat_pvalue_manual(Wlx.test.c3, label = "p = {p}", 
                     y.position = 4000, # set 4000 if geom_jitter
                     tip.length = 0.03, size=3,
                     step.increase = 0.2)+
  labs( # add labels
    title = "?????????????????? Beta-????????????????????",
    subtitle = "?? ?????????????? ???????????? ?? ???????????????? - ?????????? ??????????-????????????",
    caption = "pwc: Wilcoxon test")+
  theme(legend.title = element_blank(),
        text = element_text(size=15),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=10))