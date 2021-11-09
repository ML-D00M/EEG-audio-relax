# EEG-audio-relax
Code for an EEG study of relaxation levels in response to hypnosis audio session.

We studied how a particular hypnotic audio session might affect a person's relaxation level measured from neural activity.
For this, we designed an experiment with 2 groups of participants: 1.Hypnosis group listened to the audio with hypnosis, 2.Control group listened to the sounds of ocean waves with no speech presented. In both groups EEG activity was recorded during audio sessions (~45 min) as well as before (6-7 min) and after (6-7 min) the sessions. The experiment design was similar to 2006 study "[Quantifying Mental Relaxation with EEG for use in Computer Games](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.512.5161&rep=rep1&type=pdf)" by Lin, T. A., and L. R. John. as we compared several EEG relaxation indices based on Theta-, Alpha- & Beta-activity (9 in total).

The EEG data (32 electrodes) were first preprocessed (filtering, artifacts removal etc.) and transformed (Welch PSD) using Matlab & Brainstorm software.

In this repository you can find:
1. A Matlab script for calculation of the EEG relaxation indices
2. R scripts for the statistical analysis and visualization of the experiment results for each of 9 indicies.
a. Here we first peformed two-way mixed anova analysis for both groups
b. Then we performed one-way anova for the 1.Hypnosis group only
