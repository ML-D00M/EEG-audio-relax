%% Compute average alpha, theta and beta power indicies for 2 conditions

% The data - Power Spectal Density (calculated via Welch's method) of the
% EEG signal of: 

% a. the Main group — 32 trials, which is 2 trials per each
%   of the 15H subjects, i.e. 32 .mat files in total (initial code)

% b. the Control group — 20 trials, which is 2 trials per each
%   of the 15H subjects, i.e. 20 .mat files in total (switch 15H to 15H in
%   the code)

% "change condition index" — a marker where we need to switch between "_1" & "_3"

% load the data 

cond_1 = dir('*_1.mat'); % lists all the files with '_1' in the end of the name which are the files of the 1st contition of every subject 
%cond_2 = dir('*_2.mat');
cond_2_1 = dir('*_2.1.mat');
cond_2_2 = dir('*_2.2.mat');
cond_2_3 = dir('*_2.3.mat');
cond_2_4 = dir('*_2.4.mat');
cond_2_5 = dir('*_2.5.mat');
cond_3 = dir('*_3.mat');

% Run 3 times to calculate average alpha, theta and beta power for every condition

alpha_sub_3 = []; % change condition index
theta_sub_3 = []; % change condition index
beta_sub_3 = []; % change condition index
    
% Calculate average alpha, theta and beta power for every subject

for sub=1:length(cond_3) % change condition index

   thisSubject=cond_3(sub).name; % change condition index

   % Read the data
   path=pwd;
   data=load(thisSubject,'TF','Freqs','RowNames');
   psd = data.TF;

   % Calculate average alpha for all eletrodes
   avg_alpha_sub = mean(psd(33:49)); % indexes for alpha-band: 33 column - 8.00 Hz, 49 column - 12.00 Hz
    
   % Calculate average theta for all eletrodes
   avg_theta_sub = mean(psd(17:33)); % indexes for theta-band: 17 column - 4.00 Hz, 33 column - 8.00 Hz
   
   % Calculate average beta for all eletrodes
   avg_beta_sub = mean(psd(53:121)); % indexes for theta-band: 53 column - 13.00 Hz, 121 column - 30.00 Hz
   
   % Add the alpha and theta values in an array
   alpha_sub_3(end+1) = avg_alpha_sub; % change condition index
   theta_sub_3(end+1) = avg_theta_sub;
   beta_sub_3(end+1) = avg_beta_sub;

end

alpha_cond_3 = mean(alpha_sub_3); % change condition index
theta_cond_3 = mean(theta_sub_3); % change condition index
beta_cond_3 = mean(beta_sub_3); % change condition index

%% Add all the values in a single table 

ATB = [alpha_cond_1 theta_cond_1 beta_cond_1; alpha_cond_3 theta_cond_3 beta_cond_3];
rowNames = {'cond_1', 'cond_3'};
colNames = {'Alpha','Theta','Beta'};
ATB_tab = array2table(ATB,'RowNames',rowNames,'VariableNames',colNames)
save(['ATB_avg_values_15Hsub.mat'],'ATB');

% Save each frequency band vector for further statistical significance check
save(['alpha_15Hsub_cond_1.mat'],'alpha_sub_1');
%save(['alpha_15Hsub_cond_2.mat'],'alpha_sub_2');
save(['alpha_15Hsub_cond_2_1.mat'],'alpha_sub_2_1');
save(['alpha_15Hsub_cond_2_2.mat'],'alpha_sub_2_2');
save(['alpha_15Hsub_cond_2_3.mat'],'alpha_sub_2_3');
save(['alpha_15Hsub_cond_2_4.mat'],'alpha_sub_2_4');
save(['alpha_15Hsub_cond_2_5.mat'],'alpha_sub_2_5');
save(['alpha_15Hsub_cond_3.mat'],'alpha_sub_3');

save(['theta_15Hsub_cond_1.mat'],'theta_sub_1');
%save(['theta_15Hsub_cond_2.mat'],'theta_sub_2');
save(['theta_15Hsub_cond_2_1.mat'],'theta_sub_2_1');
save(['theta_15Hsub_cond_2_2.mat'],'theta_sub_2_2');
save(['theta_15Hsub_cond_2_3.mat'],'theta_sub_2_3');
save(['theta_15Hsub_cond_2_4.mat'],'theta_sub_2_4');
save(['theta_15Hsub_cond_2_5.mat'],'theta_sub_2_5');
save(['theta_15Hsub_cond_3.mat'],'theta_sub_3');

save(['beta_15Hsub_cond_1.mat'],'beta_sub_1');
%save(['beta_15Hsub_cond_2.mat'],'beta_sub_2');
save(['beta_15Hsub_cond_2_1.mat'],'beta_sub_2_1');
save(['beta_15Hsub_cond_2_2.mat'],'beta_sub_2_2');
save(['beta_15Hsub_cond_2_3.mat'],'beta_sub_2_3');
save(['beta_15Hsub_cond_2_4.mat'],'beta_sub_2_4');
save(['beta_15Hsub_cond_2_5.mat'],'beta_sub_2_5');
save(['beta_15Hsub_cond_3.mat'],'beta_sub_3');

%% Calculate ratios and sums of Alpha, Theta and Beta for each subject and each condition

% all the results below are 1x15H / 1x15H vectors

% Alpha/Beta ratio
rat_alpha_beta_sub_cond_1 = alpha_sub_1 ./ beta_sub_1;
%rat_alpha_beta_sub_cond_2 = alpha_sub_2 ./ beta_sub_2;
rat_alpha_beta_sub_cond_2_1 = alpha_sub_2_1 ./ beta_sub_2_1;
rat_alpha_beta_sub_cond_2_2 = alpha_sub_2_2 ./ beta_sub_2_2;
rat_alpha_beta_sub_cond_2_3 = alpha_sub_2_3 ./ beta_sub_2_3;
rat_alpha_beta_sub_cond_2_4 = alpha_sub_2_4 ./ beta_sub_2_4;
rat_alpha_beta_sub_cond_2_5 = alpha_sub_2_5 ./ beta_sub_2_5;
rat_alpha_beta_sub_cond_3 = alpha_sub_3 ./ beta_sub_3;

save(['rat_alpha_beta_15Hsub_cond_1.mat'],'rat_alpha_beta_sub_cond_1');
%save(['rat_alpha_beta_15Hsub_cond_2.mat'],'rat_alpha_beta_sub_cond_2');
save(['rat_alpha_beta_15Hsub_cond_2_1.mat'],'rat_alpha_beta_sub_cond_2_1');
save(['rat_alpha_beta_15Hsub_cond_2_2.mat'],'rat_alpha_beta_sub_cond_2_2');
save(['rat_alpha_beta_15Hsub_cond_2_3.mat'],'rat_alpha_beta_sub_cond_2_3');
save(['rat_alpha_beta_15Hsub_cond_2_4.mat'],'rat_alpha_beta_sub_cond_2_4');
save(['rat_alpha_beta_15Hsub_cond_2_5.mat'],'rat_alpha_beta_sub_cond_2_5');
save(['rat_alpha_beta_15Hsub_cond_3.mat'],'rat_alpha_beta_sub_cond_3');


% Alpha/Theta ratio
rat_alpha_theta_sub_cond_1 = alpha_sub_1 ./ theta_sub_1;
%rat_alpha_theta_sub_cond_2 = alpha_sub_2 ./ theta_sub_2;
rat_alpha_theta_sub_cond_2_1 = alpha_sub_2_1 ./ theta_sub_2_1;
rat_alpha_theta_sub_cond_2_2 = alpha_sub_2_2 ./ theta_sub_2_2;
rat_alpha_theta_sub_cond_2_3 = alpha_sub_2_3 ./ theta_sub_2_3;
rat_alpha_theta_sub_cond_2_4 = alpha_sub_2_4 ./ theta_sub_2_4;
rat_alpha_theta_sub_cond_2_5 = alpha_sub_2_5 ./ theta_sub_2_5;
rat_alpha_theta_sub_cond_3 = alpha_sub_3 ./ theta_sub_3;

save(['rat_alpha_theta_15Hsub_cond_1.mat'],'rat_alpha_theta_sub_cond_1');
%save(['rat_alpha_theta_15Hsub_cond_2.mat'],'rat_alpha_theta_sub_cond_2');
save(['rat_alpha_theta_15Hsub_cond_2_1.mat'],'rat_alpha_theta_sub_cond_2_1');
save(['rat_alpha_theta_15Hsub_cond_2_2.mat'],'rat_alpha_theta_sub_cond_2_2');
save(['rat_alpha_theta_15Hsub_cond_2_3.mat'],'rat_alpha_theta_sub_cond_2_3');
save(['rat_alpha_theta_15Hsub_cond_2_4.mat'],'rat_alpha_theta_sub_cond_2_4');
save(['rat_alpha_theta_15Hsub_cond_2_5.mat'],'rat_alpha_theta_sub_cond_2_5');
save(['rat_alpha_theta_15Hsub_cond_3.mat'],'rat_alpha_theta_sub_cond_3');

% Beta/Theta ratio
rat_beta_theta_sub_cond_1 = beta_sub_1 ./ theta_sub_1;
%rat_beta_theta_sub_cond_2 = beta_sub_2 ./ theta_sub_2;
rat_beta_theta_sub_cond_2_1 = beta_sub_2_1 ./ theta_sub_2_1;
rat_beta_theta_sub_cond_2_2 = beta_sub_2_2 ./ theta_sub_2_2;
rat_beta_theta_sub_cond_2_3 = beta_sub_2_3 ./ theta_sub_2_3;
rat_beta_theta_sub_cond_2_4 = beta_sub_2_4 ./ theta_sub_2_4;
rat_beta_theta_sub_cond_2_5 = beta_sub_2_5 ./ theta_sub_2_5;
rat_beta_theta_sub_cond_3 = beta_sub_3 ./ theta_sub_3;

save(['rat_beta_theta_15Hsub_cond_1.mat'],'rat_beta_theta_sub_cond_1');
%save(['rat_beta_theta_15Hsub_cond_2.mat'],'rat_beta_theta_sub_cond_2');
save(['rat_beta_theta_15Hsub_cond_2_1.mat'],'rat_beta_theta_sub_cond_2_1');
save(['rat_beta_theta_15Hsub_cond_2_2.mat'],'rat_beta_theta_sub_cond_2_2');
save(['rat_beta_theta_15Hsub_cond_2_3.mat'],'rat_beta_theta_sub_cond_2_3');
save(['rat_beta_theta_15Hsub_cond_2_4.mat'],'rat_beta_theta_sub_cond_2_4');
save(['rat_beta_theta_15Hsub_cond_2_5.mat'],'rat_beta_theta_sub_cond_2_5');
save(['rat_beta_theta_15Hsub_cond_3.mat'],'rat_beta_theta_sub_cond_3');


% Alpha + Beta sum
sum_alpha_beta_sub_cond_1 = alpha_sub_1 + beta_sub_1;
%sum_alpha_beta_sub_cond_2 = alpha_sub_2 + beta_sub_2;
sum_alpha_beta_sub_cond_2_1 = alpha_sub_2_1 + beta_sub_2_1;
sum_alpha_beta_sub_cond_2_2 = alpha_sub_2_2 + beta_sub_2_2;
sum_alpha_beta_sub_cond_2_3 = alpha_sub_2_3 + beta_sub_2_3;
sum_alpha_beta_sub_cond_2_4 = alpha_sub_2_4 + beta_sub_2_4;
sum_alpha_beta_sub_cond_2_5 = alpha_sub_2_5 + beta_sub_2_5;
sum_alpha_beta_sub_cond_3 = alpha_sub_3 + beta_sub_3;

save(['sum_alpha_beta_15Hsub_cond_1.mat'],'sum_alpha_beta_sub_cond_1');
%save(['sum_alpha_beta_15Hsub_cond_2.mat'],'sum_alpha_beta_sub_cond_2');
save(['sum_alpha_beta_15Hsub_cond_2_1.mat'],'sum_alpha_beta_sub_cond_2_1');
save(['sum_alpha_beta_15Hsub_cond_2_2.mat'],'sum_alpha_beta_sub_cond_2_2');
save(['sum_alpha_beta_15Hsub_cond_2_3.mat'],'sum_alpha_beta_sub_cond_2_3');
save(['sum_alpha_beta_15Hsub_cond_2_4.mat'],'sum_alpha_beta_sub_cond_2_4');
save(['sum_alpha_beta_15Hsub_cond_2_5.mat'],'sum_alpha_beta_sub_cond_2_5');
save(['sum_alpha_beta_15Hsub_cond_3.mat'],'sum_alpha_beta_sub_cond_3');

% Alpha + Theta sum
sum_alpha_theta_sub_cond_1 = alpha_sub_1 + theta_sub_1;
%sum_alpha_theta_sub_cond_2 = alpha_sub_2 + theta_sub_2;
sum_alpha_theta_sub_cond_2_1 = alpha_sub_2_1 + theta_sub_2_1;
sum_alpha_theta_sub_cond_2_2 = alpha_sub_2_2 + theta_sub_2_2;
sum_alpha_theta_sub_cond_2_3 = alpha_sub_2_3 + theta_sub_2_3;
sum_alpha_theta_sub_cond_2_4 = alpha_sub_2_4 + theta_sub_2_4;
sum_alpha_theta_sub_cond_2_5 = alpha_sub_2_5 + theta_sub_2_5;
sum_alpha_theta_sub_cond_3 = alpha_sub_3 + theta_sub_3;

save(['sum_alpha_theta_15Hsub_cond_1.mat'],'sum_alpha_theta_sub_cond_1');
%save(['sum_alpha_theta_15Hsub_cond_2.mat'],'sum_alpha_theta_sub_cond_2');
save(['sum_alpha_theta_15Hsub_cond_2_1.mat'],'sum_alpha_theta_sub_cond_2_1');
save(['sum_alpha_theta_15Hsub_cond_2_2.mat'],'sum_alpha_theta_sub_cond_2_2');
save(['sum_alpha_theta_15Hsub_cond_2_3.mat'],'sum_alpha_theta_sub_cond_2_3');
save(['sum_alpha_theta_15Hsub_cond_2_4.mat'],'sum_alpha_theta_sub_cond_2_4');
save(['sum_alpha_theta_15Hsub_cond_2_5.mat'],'sum_alpha_theta_sub_cond_2_5');
save(['sum_alpha_theta_15Hsub_cond_3.mat'],'sum_alpha_theta_sub_cond_3');

% Theta + Beta sum
sum_theta_beta_sub_cond_1 = theta_sub_1 + beta_sub_1;
%sum_theta_beta_sub_cond_2 = theta_sub_2 + beta_sub_2;
sum_theta_beta_sub_cond_2_1 = theta_sub_2_1 + beta_sub_2_1;
sum_theta_beta_sub_cond_2_2 = theta_sub_2_2 + beta_sub_2_2;
sum_theta_beta_sub_cond_2_3 = theta_sub_2_3 + beta_sub_2_3;
sum_theta_beta_sub_cond_2_4 = theta_sub_2_4 + beta_sub_2_4;
sum_theta_beta_sub_cond_2_5 = theta_sub_2_5 + beta_sub_2_5;
sum_theta_beta_sub_cond_3 = theta_sub_3 + beta_sub_3;

save(['sum_theta_beta_15Hsub_cond_1.mat'],'sum_theta_beta_sub_cond_1');
%save(['sum_theta_beta_15Hsub_cond_2.mat'],'sum_theta_beta_sub_cond_2');
save(['sum_theta_beta_15Hsub_cond_2_1.mat'],'sum_theta_beta_sub_cond_2_1');
save(['sum_theta_beta_15Hsub_cond_2_2.mat'],'sum_theta_beta_sub_cond_2_2');
save(['sum_theta_beta_15Hsub_cond_2_3.mat'],'sum_theta_beta_sub_cond_2_3');
save(['sum_theta_beta_15Hsub_cond_2_4.mat'],'sum_theta_beta_sub_cond_2_4');
save(['sum_theta_beta_15Hsub_cond_2_5.mat'],'sum_theta_beta_sub_cond_2_5');
save(['sum_theta_beta_15Hsub_cond_3.mat'],'sum_theta_beta_sub_cond_3');

% Alpha + Theta + Beta sum
sum_alpha_theta_beta_sub_cond_1 = alpha_sub_1 + theta_sub_1 + beta_sub_1;
%sum_alpha_theta_beta_sub_cond_2 = alpha_sub_2 + theta_sub_2 + beta_sub_2;
sum_alpha_theta_beta_sub_cond_2_1 = alpha_sub_2_1 + theta_sub_2_1 + beta_sub_2_1;
sum_alpha_theta_beta_sub_cond_2_2 = alpha_sub_2_2 + theta_sub_2_2 + beta_sub_2_2;
sum_alpha_theta_beta_sub_cond_2_3 = alpha_sub_2_3 + theta_sub_2_3 + beta_sub_2_3;
sum_alpha_theta_beta_sub_cond_2_4 = alpha_sub_2_4 + theta_sub_2_4 + beta_sub_2_4;
sum_alpha_theta_beta_sub_cond_2_5 = alpha_sub_2_5 + theta_sub_2_5 + beta_sub_2_5;
sum_alpha_theta_beta_sub_cond_3 = alpha_sub_3 + theta_sub_3 + beta_sub_3;

save(['sum_alpha_theta_beta_15Hsub_cond_1.mat'],'sum_alpha_theta_beta_sub_cond_1');
%save(['sum_alpha_theta_beta_15Hsub_cond_2.mat'],'sum_alpha_theta_beta_sub_cond_2');
save(['sum_alpha_theta_beta_15Hsub_cond_2_1.mat'],'sum_alpha_theta_beta_sub_cond_2_1');
save(['sum_alpha_theta_beta_15Hsub_cond_2_2.mat'],'sum_alpha_theta_beta_sub_cond_2_2');
save(['sum_alpha_theta_beta_15Hsub_cond_2_3.mat'],'sum_alpha_theta_beta_sub_cond_2_3');
save(['sum_alpha_theta_beta_15Hsub_cond_2_4.mat'],'sum_alpha_theta_beta_sub_cond_2_4');
save(['sum_alpha_theta_beta_15Hsub_cond_2_5.mat'],'sum_alpha_theta_beta_sub_cond_2_5');
save(['sum_alpha_theta_beta_15Hsub_cond_3.mat'],'sum_alpha_theta_beta_sub_cond_3');


