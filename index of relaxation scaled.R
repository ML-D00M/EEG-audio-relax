library(scales)
library(R.matlab)
library(tidyverse)
library(rstatix)
library(ggpubr)

# set working directory
setwd("E:/PAOLO/biolabs/audio - Relax/PSD EEG Estimation/C, CP, P, PO, O alpha sub16-10")


## DATA PREPARATION ====

# Load the data saved from matlab

# 1. Baseline data: 
# 1st group, Before
sum_alpha_theta_16sub_1=readMat('sum_alpha_theta_16sub_cond_1.mat')$sum.alpha.theta.sub
sum_alpha_theta_16sub_1=as.vector(sum_alpha_theta_16sub_1)
# 2nd group, Before
sum_alpha_theta_10Csub_1=readMat('sum_alpha_theta_10Csub_cond_1.mat')$sum.alpha.theta.sub 
sum_alpha_theta_10Csub_1=as.vector(sum_alpha_theta_10Csub_1)

# 2.1 Test data #1:
# 2nd group, after
sum_alpha_theta_10Csub_3=readMat('sum_alpha_theta_10Csub_cond_3.mat')$sum.alpha.theta.sub 
sum_alpha_theta_10Csub_3=as.vector(sum_alpha_theta_10Csub_3)

# 2.2 Test data #2:
# 1st group, Hypno: Intro
sum_alpha_theta_15Hsub_2_1=readMat('sum_alpha_theta_15Hsub_cond_2_1.mat')$sum.alpha.theta.sub
sum_alpha_theta_15Hsub_2_1=as.vector(sum_alpha_theta_15Hsub_2_1)

# Create a baseline data frame
id <- as.numeric(sprintf("% d", 1:26))
relax <- c(sum_alpha_theta_16sub_1, sum_alpha_theta_10Csub_1)
df <- data.frame(id,relax)

# Normalize the Relax variable by rescaling it.
# Multiply its super small values by 1e13 to have the values >1
# This will allow to calculate statistics
df$relax <- df$relax*1e13

# Create a new data frame #1 - for further testing
id <- sprintf("% d", 1:10)
relax <- sum_alpha_theta_10Csub_3
df_new <- data.frame(id,relax)

# Create a new data frame #2 - for further testing
id <- sprintf("% d", 1:15)
relax <- sum_alpha_theta_15Hsub_2_1
df_new2 <- data.frame(id,relax)

group.colors <- c("Hypnosis" = "#C5610D", "Control" = "#087676")

## CHECK STATISTICS ====

# Check outliers
out=df %>%
  identify_outliers(relax)
out

# Remove extreme outliers - Values above Q3 + 3xIQR or below Q1 - 3xIQR 
#df <- df[-c(8,18,20), ]
# to Remove all the outliers - Values above Q3 + 1.5xIQR or below Q1 - 1.5xIQR
#tbr=which(df$relax%in%out$relax)
#df=df[-tbr,]

# Check the outliers again and remove again if necessary
#out=df %>%
#  identify_outliers(relax)
#out

# Summary statistics
stats <- df %>%
  get_summary_stats(relax, type = "full")
stats

# Visualize the basic statistics of our baseline data frame
raw_plot = df
raw_plot$id <- as.numeric(df$id)
raw_plot$group = "Hypnosis"
raw_plot$group[raw_plot$id>16] <- "Control"
raw_plot$relax <- as.integer(df$relax)
ggbarplot(raw_plot, x='id', y='relax',
          xlab="Participant number",ylab='Relaxation index',
          fill="group",
          #ylim=c(0, 105),
          label=T)+
  scale_fill_manual(values=group.colors)+ #set colors
  labs( # add labels
    title = "Alpha+Theta-????????????????????",
    subtitle = "?????????????? ???????????????????? - ???? ?????????????????????????? ??????????-????????????")+
  theme(text = element_text(size=15),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=15))

# Take basic statistical parameters for further index calculations
q1 = as.numeric(stats$q1)
q3 = as.numeric(stats$q3)
iqr = as.numeric(stats$iqr)



## CREATE RELAX INDICES - BASIC ==== 

# 1. 0 min group indices
min = subset(df, relax<=0)
i_min <- min
i_min$relax <- rescale(i_min$relax, to = 0)
i_min$relax <- as.integer(i_min$relax)

# 2. 0-25 low group indices
low = subset(df, relax>0 & relax<=q1)
i_low <- low
i_low$relax <- rescale(i_low$relax, to = c(1, 25))
i_low$relax <- as.integer(i_low$relax)

# 3. 26-50 mid group indices
mid = subset(df, relax>q1 & relax<=q3)
i_mid <- mid
i_mid$relax <- rescale(i_mid$relax, to = c(26, 50))
i_mid$relax <- as.integer(i_mid$relax)

# 4. 56-75 mid_high group indices
mid_high = subset(df, relax>q3 & relax<=q3+1.5*iqr)
i_mid_high <- mid_high
i_mid_high$relax <- rescale(i_mid_high$relax, to = c(56, 75))
mid_high$relax <- as.integer(mid_high$relax)

# 5. 76-99 high group indices
high = subset(df, relax>q3+1.5*iqr & relax<=q3+3*iqr)
i_high <- high
i_high$relax <- rescale(i_high$relax, to = c(75, 99))
i_high$relax <- as.integer(i_high$relax)

# 6. 100 max group indices
max = subset(df, relax>q3+3*iqr)
i_max <- max
i_max$relax <- rescale(i_max$relax, to = 100)
i_max$relax <- as.integer(i_max$relax)

# Concatinate df with indices back into one df
i_relax <- rbind(i_min, i_low, i_mid, i_mid_high, i_high, i_max)
i_relax <- i_relax[order(i_relax$id),]

# Visualize the df with Relaxation indices
relax_plot = i_relax
relax_plot$id <- as.numeric(df$id)
relax_plot$group = "Hypnosis"
relax_plot$group[relax_plot$id>16] <- "Control"
ggbarplot(relax_plot, x='id', y='relax',
          xlab="Participant number",ylab='Relaxation index',
          fill="group",
          ylim=c(0, 105),
          label=T)+
  scale_fill_manual(values=group.colors)+ #set colors
  labs( # add labels
    title = "?????????????? ???????????????????????? ???? ???????? Alpha+Theta-????????????????????",
    subtitle = "?????????????? ???????????????????? - ???? ?????????????????????????? ??????????-????????????")+
  theme(text = element_text(size=15),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=15))



## TEST NEW INDICIES #1 ====

# Check the Test data #1: Control: After

# Create a test data frame
df_test1 <- df_new
df_test1$relax <- df_test1$relax*1e13

# Create an empty dataset for our indices
i_relax_test1 <- min 

id1 = 0

# Check each Alpha+Theta value in our new dataset
for (j in df_test1$relax)
{
  # Take one new Alpha+Theta value for testing
  at = j
  
  # Give it an id for further addition to the baseline data
  id1 = id1+1
  
  # Compare it with our baseline statistics and 
  # turn into the Relaxation index (i.e. map onto the 0-100 scale)
  if (at<=0) {
    min[nrow(min) + 1,] = c(id1,at)
    i_min <- min
    i_min$relax <- rescale(i_min$relax, to = 0)
    i_min$relax <- as.integer(i_min$relax)
    i = as.numeric(i_min[i_min$id==id1,]$relax)
  } else if (at>0 & at<=q1) {
    low[nrow(low) + 1,] = c(id1,at)
    i_low <- low
    i_low$relax <- rescale(i_low$relax, to = c(1, 25))
    i_low$relax <- as.integer(i_low$relax)
    i = as.numeric(i_low[i_low$id==id1,]$relax)
  } else if (at>q1 & at<=q3) {
    mid[nrow(mid) + 1,] = c(id1,at)
    i_mid <- mid
    i_mid$relax <- rescale(i_mid$relax, to = c(26, 50))
    i_mid$relax <- as.integer(i_mid$relax)
    i = as.numeric(i_mid[i_mid$id==id1,]$relax)
  } else if (at>q3 & at<=q3+1.5*iqr) {
    mid_high[nrow(mid_high) + 1,] = c(id1,at)
    i_mid_high <- mid_high
    i_mid_high$relax <- rescale(i_mid_high$relax, to = c(51, 75))
    i_mid_high$relax <- as.integer(i_mid_high$relax)
    i = as.numeric(i_mid_high[i_mid_high$id==id1,]$relax)
  } else if (at>q3+1.5*iqr & at<=q3+3*iqr) {
    high[nrow(high) + 1,] = c(id1,at)
    i_high <- high
    i_high$relax <- rescale(i_high$relax, to = c(76, 99))
    i_high$relax <- as.integer(i_high$relax)
    i = as.numeric(i_high[i_high$id==id1,]$relax)
  } else {
    max[nrow(max) + 1,] = c(id1,at)
    i_max <- max
    i_max$relax <- rescale(i_max$relax, to = 100)
    i_max$relax <- as.integer(i_max$relax)
    i = as.numeric(i_max[i_max$id==id1,]$relax)
  }
  
  print(i) # this is our Relaxation index
  
  # Concatinate df with indices back into one df
  i_relax_test1[nrow(i_relax_test1) + 1,] = c(id1,i)
  
}

# Visualize our new Relaxation indices
relax_plot1 = i_relax_test1
relax_plot1$id <- as.numeric(df_new1$id)
ggbarplot(relax_plot1, x='id', y='relax',
          xlab="Participant number",ylab='Relaxation index',
          fill="#087676",
          ylim=c(0, 105),
          label=T)+
  #geom_jitter(size=2, width=0.1, alpha=0.9)+
  labs( # add labels
    title = "?????????????? ???????????????????????? ???? ???????? Alpha+Theta-????????????????????",
    subtitle = "?????????? ?????????????????????? ??????????-???????????? ?? ?????????? ????????")+
  theme(text = element_text(size=15),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=15))



## TEST NEW INDICIES #2 ====

# Check the Test data #2 - Hypno: Intro

# Create a test data frame
df_test2 <- df_new2
df_test2$relax <- df_test2$relax*1e13

# Create an empty dataset for our indices
i_relax_test2 <- min 

id2 = 0

# Check each Alpha+Theta value in our new dataset
for (j in df_test2$relax)
{
  # Take one new Alpha+Theta value for testing
  at = j
  
  # Give it an id for further addition to the baseline data
  id2 = id2+1
  
  # Compare it with our baseline statistics and 
  # turn into the Relaxation index (i.e. map onto the 0-100 scale)
  if (at<=0) {
    min[nrow(min) + 1,] = c(id2,at)
    i_min <- min
    i_min$relax <- rescale(i_min$relax, to = 0)
    i_min$relax <- as.integer(i_min$relax)
    i = as.numeric(i_min[i_min$id==id2,]$relax)
  } else if (at>0 & at<=q1) {
    low[nrow(low) + 1,] = c(id2,at)
    i_low <- low
    i_low$relax <- rescale(i_low$relax, to = c(1, 25))
    i_low$relax <- as.integer(i_low$relax)
    i = as.numeric(i_low[i_low$id==id2,]$relax)
  } else if (at>q1 & at<=q3) {
    mid[nrow(mid) + 1,] = c(id2,at)
    i_mid <- mid
    i_mid$relax <- rescale(i_mid$relax, to = c(26, 50))
    i_mid$relax <- as.integer(i_mid$relax)
    i = as.numeric(i_mid[i_mid$id==id2,]$relax)
  } else if (at>q3 & at<=q3+1.5*iqr) {
    mid_high[nrow(mid_high) + 1,] = c(id2,at)
    i_mid_high <- mid_high
    i_mid_high$relax <- rescale(i_mid_high$relax, to = c(51, 75))
    i_mid_high$relax <- as.integer(i_mid_high$relax)
    i = as.numeric(i_mid_high[i_mid_high$id==id2,]$relax)
  } else if (at>q3+1.5*iqr & at<=q3+3*iqr) {
    high[nrow(high) + 1,] = c(id2,at)
    i_high <- high
    i_high$relax <- rescale(i_high$relax, to = c(76, 99))
    i_high$relax <- as.integer(i_high$relax)
    i = as.numeric(i_high[i_high$id==id2,]$relax)
  } else {
    max[nrow(max) + 1,] = c(id2,at)
    i_max <- max
    i_max$relax <- rescale(i_max$relax, to = 100)
    i_max$relax <- as.integer(i_max$relax)
    i = as.numeric(i_max[i_max$id==id2,]$relax)
  }
  
  print(i) # this is our Relaxation index
  
  # Concatinate df with indices back into one df
  i_relax_test2[nrow(i_relax_test2) + 1,] = c(id2,i)
  
}

# Visualize our new Relaxation indices
relax_plot2 = i_relax_test2
relax_plot2$id <- as.numeric(df_new2$id)
ggbarplot(relax_plot2, x='id', y='relax',
          xlab="Participant number",ylab='Relaxation index',
          fill="#F9A057",
          ylim=c(0, 105),
          label=T)+
  #geom_jitter(size=2, width=0.1, alpha=0.9)+
  labs( # add labels
    title = "?????????????? ???????????????????????? ???? ???????? Alpha+Theta-????????????????????",
    subtitle = "???? ?????????? ?????????????? 5 ?????????? ??????????-???????????? ?? ????????????????")+
  theme(text = element_text(size=15),
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=15))



## CALCULATE A NEW INDEX ====

# Here we take a new Alpha+Theta power value, turn it into a number 
# and map onto our 0-100 relaxation scale

# Create a test data frame - with new data we obtained
df_test <- df_new

# Multiply by the fixed value to calculate statistics
df_test$relax <- df_test$relax*1e13

# Take one new Alpha+Theta value for testing
at = df_test$relax[1]

# Give it an id for further addition to the baseline data
id_at = as.numeric(tail(as.vector(df$id), n=1))+1

# Compare it with our baseline statistics and 
# turn into the Relaxation index (i.e. map onto the 0-100 scale)
if (at<=0) {
  min[nrow(min) + 1,] = c(id_at,at)
  i_min <- min
  i_min$relax <- rescale(i_min$relax, to = 0)
  i_min$relax <- as.integer(i_min$relax)
  i = as.numeric(i_min[i_min$id==id_at,]$relax)
} else if (at>0 & at<=q1) {
  low[nrow(low) + 1,] = c(id_at,at)
  i_low <- low
  i_low$relax <- rescale(i_low$relax, to = c(1, 25))
  i_low$relax <- as.integer(i_low$relax)
  i = as.numeric(i_low[i_low$id==id_at,]$relax)
} else if (at>q1 & at<=q3) {
  mid[nrow(mid) + 1,] = c(id_at,at)
  i_mid <- mid
  i_mid$relax <- rescale(i_mid$relax, to = c(26, 50))
  i_mid$relax <- as.integer(i_mid$relax)
  i = as.numeric(i_mid[i_mid$id==id_at,]$relax)
} else if (at>q3 & at<=q3+1.5*iqr) {
  mid_high[nrow(mid_high) + 1,] = c(id_at,at)
  i_mid_high <- mid_high
  i_mid_high$relax <- rescale(i_mid_high$relax, to = c(51, 75))
  i_mid_high$relax <- as.integer(i_mid_high$relax)
  i = as.numeric(i_mid_high[i_mid_high$id==id_at,]$relax)
} else if (at>q3+1.5*iqr & at<=q3+3*iqr) {
  high[nrow(high) + 1,] = c(id_at,at)
  i_high <- high
  i_high$relax <- rescale(i_high$relax, to = c(76, 99))
  i_high$relax <- as.integer(i_high$relax)
  i = as.numeric(i_high[i_high$id==id_at,]$relax)
} else {
  max[nrow(max) + 1,] = c(id_at,at)
  i_max <- max
  i_max$relax <- rescale(i_max$relax, to = 100)
  i_max$relax <- as.integer(i_max$relax)
  i = as.numeric(i_max[i_max$id==id_at,]$relax)
}

iR = i
iR # this is our new Relaxation index



## UPDATE THE BASELINE DATA ====

# After we obtained a new baseline Alpha+Theta measure, 
# we shoul add it to our baseline dataframe - df
df[nrow(df) + 1,] = c(id_at,at)

# Save the df in a file
save(df,file="Relax_baseline_raw.Rda")

# Later we can load the file 
#load("Relax_baseline_raw.Rda")