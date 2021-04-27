%% treatment Data (EMG) Hometech project 
% code 1 : Pre-processed EMG
% author : Robin Macchi 

% clear workspace...
clear all
close all
clc 

%% STEP 1 - EMG raw data processing ----
% import data EMG
prompt = 'Select the patheway of the EMG data'
dirname=uigetdir % Where are your files (EMG) located? 
cd(dirname)
list=dir(dirname);
list=list(3:end);

for k=1:length(list)
cd(dirname)
filename=list(k).name;
opts = delimitedTextImportOptions("NumVariables", 15);
opts.DataLines = [3, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["time", "pectoralis_major_nd", "deltoidus_nd", "infraspinatus_nd", "biceps_brachialis_nd", "triceps_brachialis_nd", "flexor_ulnaris_carpi_nd", "extensor_ulnaris_carpi_nd", "pectoralis_major_d", "deltoidus_d", "infraspinatus_d", "biceps_brachialis_d", "triceps_brachialis_d", "flexor_ulnaris_carpi_d", "extensor_ulnaris_carpi_d"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
data = readtable(filename, opts);
data = table2array(data);
clear opts

% dependent of kinematic data (some files begin to 9 secondes) 
if k==2 || k==5
data=data(9000:end,:);
else
end

% Filtering

[B,A]=butter(4,[20 500].*0.802/1000); % regular bandpass filter, ex: 20 500 
for i=2:15
data_detrend(:,i)=detrend(data(:,i));
data_filtre(:,i)=filtfilt(B,A,data_detrend(:,i));
end
data_filtre(:,1)=data(:,1);
t=data_filtre(:,1);             

% Full-wave rectification
for i=2:15
    val_EMG_abs(:,i)=abs(data_filtre(:,i));
end
val_EMG_abs(:,1)=t;

  % Low-pass (creating the linear envelope of the signal) Butterworth zero-phase filter design
EmgFreq=1000;
[D,C]=butter(2,5.*0.802/EmgFreq);  %low pass filter, ex: 5 hz
for i=2:15
EmgEnv(:,i)=filtfilt(D,C,val_EMG_abs(:,i));
end
EmgEnv(:,1)=t;

%% step 2 - marker kinematic data processing ----

% load marker file (R_MC1)
if k==1
    prompt = 'Select the patheway of the R_MC1 marker data'
dirname_marker=uigetdir % Where are your files (kinematic marker) located? 
else
end
cd(dirname_marker)
list_marker=dir(dirname_marker);
list_marker=list_marker(3:end);

filename_marker=list_marker(k).name;
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["time", "x", "y", "z"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
data_marker = readtable(filename_marker, opts);
data_marker = table2array(data_marker);
clear opts

R_MC1=data_marker(:,2:4);
t2=data_marker(:,1); % time of the kinematic data

% if the data begin to 9s

if k==2 || k==5
t2=t2(2700:end,:);
R_MC1=R_MC1(2700:end,:);
else
end

% visualization

subplot(3,2,1)
plot(t,val_EMG_abs(:,2)); %  EMG file -> filtered and rectified from pectoralis_major_nd muscle

subplot(3,2,2)
plot(t2,R_MC1);  % kinematic data from R_MC1 marker 

subplot(3,2,3)
plot(t,data(:,2)); % raw EMG file from pectoralis_major_nd muscle

% test of synchronisation : signal superimposing (EMG + kinematic marker R_MC1)

R_MC1_bis=R_MC1(:,3); % vertical axis 
R_MC1_bis=R_MC1_bis*1500;   % amplified signal for a better visualization 

plot(t,val_EMG_abs(:,14)); % flexor_ulnaris_cerpi_d muscle
hold on 
plot(t2,R_MC1_bis(:,1));  % R_MC1 (vertical axis)

% to get the R_MC1 peaks

% Reduce file size

clear R_MC1_bis
clear deb1 fin1 x

prompt = 'Select the start (just before the first peak) and the end of the trial'
figure('name','R_MC1','NumberTitle','off');
p=plot(t2,R_MC1); 
[x,~] = ginput(2);

deb1=x(1);   % choose as the start, just before the movement begins (before the first peak)
deb1=find(t2>=deb1);
deb1=deb1(1);
    
fin1=x(2);   % choose as the end, just after the last peak
fin1=find(t2>=fin1); 
fin1=fin1(1);

R_MC1=R_MC1([deb1:fin1],:); 
t2=t2([deb1:fin1],:);

% move mean to smooth the peaks 

R_MC1=R_MC1(:,3);
window_size=6;
R_MC1_bis(1:window_size) = R_MC1(1:window_size);
for i=(window_size+1):length(R_MC1)-(window_size+1)
    R_MC1_bis(i) = mean(R_MC1(i-window_size:i+window_size));
end
R_MC1_bis(length(R_MC1)-window_size:length(R_MC1)) = R_MC1(length(R_MC1)-window_size:length(R_MC1));

% check the smooth peaks

figure;
plot(t2,R_MC1)
hold on 
plot(t2,R_MC1_bis)

%% step 3 - identification of the different cycle phases ---

% 3 phases : analysis phase ; up-swing phase ; down-swing phase identified from the marker kinematic
prompt = 'Select the threshold for the peak detection'
p=plot(t2,R_MC1);  % choose a threshold to indentify the peaks 
[x,y] = ginput(1);
seuil_kinematic=y(1);  

% get peak_max and peak_min of the down-swing phase, so juste after the
% peak max (based on R_MC1 kinematic from the vertical axis)
n_max = 1;
n_min = 1;
mode = 1; %1 = peak max, 2 = peak min
for i=2:length(R_MC1_bis)-1
    if(mode == 1)
        if(R_MC1_bis(i)> seuil_kinematic && R_MC1_bis(i)>R_MC1_bis(i+1) && R_MC1_bis(i)>R_MC1_bis(i-1))
            R_MC1_pic_max(n_max)=R_MC1_bis(i);
            t_R_MC1_pic_max(n_max)=t2(i);
            p_R_MC1_max(n_max)=i;
            n_max=n_max+1;
            mode = 2;
        end
    else
        if(R_MC1_bis(i)< seuil_kinematic && R_MC1_bis(i)<R_MC1_bis(i+1) && R_MC1_bis(i)<R_MC1_bis(i-1))
            R_MC1_pic_min(n_min)=R_MC1_bis(i);
            t_R_MC1_pic_min(n_min)=t2(i);
            p_R_MC1_min(n_min)=i;
            n_min=n_min+1;
            mode = 1;
        end
    end
end

% find peak min of the start of the up-swing phase (so just before the peak max)

n_max=1;

for n_max=1:length(p_R_MC1_min)
for i=2:length(R_MC1_bis)-1
    if (R_MC1_bis(i)< seuil_kinematic && R_MC1_bis(i)<R_MC1_bis(i+1) && R_MC1_bis(i)<R_MC1_bis(i-1) && i<p_R_MC1_max(n_max))
        R_MC1_pic_deb(:,n_max)=R_MC1_bis(i);
        R_MC1_pos_deb(:,n_max)=i;
        t_R_MC1_pos_deb(:,n_max)=t2(i);
    end
end
end

% check the peaks 

figure;
plot(t2,R_MC1);
hold on 
plot(t2,R_MC1_bis);
hold on 
scatter(t_R_MC1_pic_max,R_MC1_pic_max);
hold on 
scatter(t_R_MC1_pic_min,R_MC1_pic_min);
hold on 
scatter(t_R_MC1_pos_deb,R_MC1_pic_deb);

% get the peak of the last cycle

j=length(R_MC1);
r=length(p_R_MC1_min);
r2=p_R_MC1_min(r);

% get time of each cycle (peak_min(i) up to au pic min (i+1))

n=2;
for i=1:r-1;
    cycle.(['number' num2str(1)])=t2(1:p_R_MC1_min(1)); % first cycle
    cycle.(['number' num2str(n)])=t2(p_R_MC1_min(i):p_R_MC1_min(i+1)); % other cycles
    n=n+1;
end

% time of the analysis phase

n=2;
for i=1:r-1;
    cycle_preparation.(['number' num2str(1)])=t2(1:R_MC1_pos_deb(1)); % first analysis phase 
    cycle_preparation.(['number' num2str(n)])=t2(p_R_MC1_min(i):R_MC1_pos_deb(i+1)); % other analysis phases
    n=n+1;
end

% time of the up-swing phase

n=1;
for i=1:r;
    cycle_ascendant.(['number' num2str(n)])=t2(R_MC1_pos_deb(i):p_R_MC1_max(i));
    n=n+1;
end

% time of the dow-swing phase

n=1;
for i=1:r;
    cycle_descendant.(['number' num2str(n)])=t2(p_R_MC1_max(i):p_R_MC1_min(i));
    n=n+1;
end


% get frames of the global cycles of the EMG data

cycle=struct2cell(cycle);
for i=1:numel(cycle)
    P0_cycle=find(t>=cycle{i,1}(1));
    P0_cycle_(i,:)=P0_cycle(1);
    P1_cycle=find(t>=max(cycle{i,1}));
    P1_cycle_(i,:)=P1_cycle(1);
    clear P0_cycle P1_cycle
end

% get frames of analysis phases for each cycle (EMG data)

cycle_preparation=struct2cell(cycle_preparation);
for i=1:numel(cycle_preparation)
    P0_cycle_preparation=find(t>=cycle_preparation{i,1}(1));
    P0_cycle_preparation_(i,:)=P0_cycle_preparation(1);
    P1_cycle_preparation=find(t>=max(cycle_preparation{i,1}));
    P1_cycle_preparation_(i,:)=P1_cycle_preparation(1);
    clear P0_cycle_preparation P1_cycle_preparation
end

% get frames of the up-swing phases for each cycle (EMG data)

cycle_ascendant=struct2cell(cycle_ascendant);
for i=1:numel(cycle_ascendant)
    P0_cycle_ascendant=find(t>=cycle_ascendant{i,1}(1));
    P0_cycle_ascendant_(i,:)=P0_cycle_ascendant(1);
    P1_cycle_ascendant=find(t>=max(cycle_ascendant{i,1}));
    P1_cycle_ascendant_(i,:)=P1_cycle_ascendant(1);
    clear P0_cycle_ascendant P1_cycle_ascendant
end

% get frames of the down-swing phases for each cycle (EMG data)

cycle_descendant=struct2cell(cycle_descendant);
for i=1:numel(cycle_descendant);
    P0_cycle_descendant=find(t>=cycle_descendant{i,1}(1));
    P0_cycle_descendant_(i,:)=P0_cycle_descendant(1);
    P1_cycle_descendant=find(t>=max(cycle_descendant{i,1}));
    P1_cycle_descendant_(i,:)=P1_cycle_descendant(1);
    clear P0_cycle_descendant P1_cycle_descendant
end

%% step 4 - maximum voluntary contraction EMG files processing ---

% open max EMG file (choose the file corresponds to the day of the trial performing)
if k==1
    prompt = 'Select the patheway of the EMG maximal voluntary contraction data'
dirname_max_EMG=uigetdir  % Where are your files (max_EMG) located? 
else
end
cd(dirname_max_EMG)

if k==1 || k==4 || k==6
    filename_max_EMG=dirname_max_EMG+"\max_EMG_day_1_EMG.txt";
else
    filename_max_EMG=dirname_max_EMG+"\max_EMG_day_2_EMG.txt";
end
opts = delimitedTextImportOptions("NumVariables", 15);
opts.DataLines = [3, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["time", "pectoralis_major_nd", "deltoidus_nd", "infraspinatus_nd", "biceps_brachialis_nd", "triceps_brachialis_nd", "flexor_ulnaris_carpi_nd", "extensor_ulnaris_carpi_nd", "pectoralis_major_d", "deltoidus_d", "infraspinatus_d", "biceps_brachialis_d", "triceps_brachialis_d", "flexor_ulnaris_carpi_d", "extensor_ulnaris_carpi_d"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
maxEMG = readtable(filename_max_EMG, opts);
maxEMG = table2array(maxEMG);
clear opts

% Filter the max EMG files (Pass band filter)
n=1;
[B,A]=butter(4,[20 500].*0.802/1000); % regular bandpass filter, ex: 20 500 
for i=2:15;
maxEMG_detrend(:,i)=detrend(maxEMG(:,i));
maxEMG_filtre(:,i)=filtfilt(B,A,maxEMG_detrend(:,i));
end

% rectified signal 

i=1;
for i=2:15
    val_EMG_max(:,i)=abs(maxEMG_filtre(:,i));
end

% linear envelope of the max EMG file (low pass filter 5 Hz)

EmgFreq=1000;
[D,C]=butter(2,5.*0.802/EmgFreq);  %low pass filter, ex: 5 hz

for i=2:15
EmgEnv_max(:,i)=filtfilt(D,C,val_EMG_max(:,i));
end

% Find the max value of EMG_max 

i=1;
for i=2:15
    Max_val_EMG(:,i)=max(EmgEnv_max(:,i));
end

%% step 5 - EMG cycles processing ---

% normalized the EMG data EMG from the maximal voluntary contraction

i=1;
for i=2:15
val_EMG_norm(:,i)=EmgEnv(:,i)/Max_val_EMG(:,i);
end

% get normalized EMG data for each cycle (global cycle)

for m=2:length(val_EMG_norm(1,:))
for i=1:numel(cycle)
    cycle_EMG.(['muscle' num2str(m)]).(['cycle' num2str(i)])=val_EMG_norm(P0_cycle_(i):P1_cycle_(i),m);
end
end

% get normalized EMG data of the analysis phase for each cycle

for m=2:length(val_EMG_norm(1,:))
for i=1:numel(cycle_preparation)
    cycle_EMG_preparation.(['muscle' num2str(m)]).(['cycle' num2str(i)])=val_EMG_norm(P0_cycle_preparation_(i):P1_cycle_preparation_(i),m);
end
end

% get normalized EMG data of the up-swing phase for each cycle

for m=2:length(val_EMG_norm(1,:))
for i=1:numel(cycle_ascendant)
    cycle_EMG_ascendant.(['muscle' num2str(m)]).(['cycle' num2str(i)])=val_EMG_norm(P0_cycle_ascendant_(i):P1_cycle_ascendant_(i),m);
end
end

% get normalized EMG data of the down-swing phase for each cycle

for m=2:length(val_EMG_norm(1,:))
for i=1:numel(cycle_descendant)
    cycle_EMG_descendant.(['muscle' num2str(m)]).(['cycle' num2str(i)])=val_EMG_norm(P0_cycle_descendant_(i):P1_cycle_descendant_(i),m);
end
end

%% step 6 - compute iEMG (trapezoid method) for each cycle phase ---

% get trapezoids of each global cycle 

for m=2:length(val_EMG_norm(1,:))
    for i=1:numel(cycle)
        n=1;
        for j=1:length(cycle_EMG.(['muscle' num2str(m)]).(['cycle' num2str(i)]))-1
iEMG_cycle.(['muscle' num2str(m)]).(['cycle' num2str(i)])(n)=((cycle_EMG.(['muscle' num2str(m)]).(['cycle' num2str(i)])(j+1)+cycle_EMG.(['muscle' num2str(m)]).(['cycle' num2str(i)])(j))*(1/1000))/2;      
    n=n+1;
        end
    end
end

% sum trapezoids to get iEMG of global cycle 

for m=2:length(val_EMG_norm(1,:))
    for i=1:numel(cycle)
sum_iEMG_cycle.(['muscle' num2str(m)]).(['cycle' num2str(i)])=sum(iEMG_cycle.(['muscle' num2str(m)]).(['cycle' num2str(i)]));
    end
end

%% step 7 - save file iEMG ---

% Creating output file (iEMG)
resultats_total=sum_iEMG_cycle;

% save in .mat
if k==1
    prompt = 'Select the patheway of the backup files'
dirname_save=uigetdir   % Where do you save the output files ? 
else
end
cd(dirname_save)

save(string(filename(1:10))+'_results_iEMG.mat','resultats_total');

% save in .txt

j=1;
for m=2:length(val_EMG_norm(1,:))
    for i=1:numel(cycle)
for k=1:length(resultats_total.(['muscle' num2str(m)]).(['cycle' num2str(i)]))
    resultats_bis(k,j)=resultats_total.(['muscle' num2str(m)]).(['cycle' num2str(i)]);
    j=j+1;
end
    end
end


resultats_bis=table([resultats_bis(:,1)],[resultats_bis(:,2)],[resultats_bis(:,3)],[resultats_bis(:,4)]...
    ,[resultats_bis(:,5)],[resultats_bis(:,6)],[resultats_bis(:,7)],[resultats_bis(:,8)],[resultats_bis(:,9)],[resultats_bis(:,10)]...
    ,[resultats_bis(:,11)],[resultats_bis(:,12)],[resultats_bis(:,13)],[resultats_bis(:,14)]);
resultats_bis.Properties.VariableNames = {'pectoralis_major_nd' 'deltoidus_nd' 'infraspinatus_nd' 'biceps_brachialis_nd' 'triceps_brachialis_nd' 'flexor_ulnaris_cerpi_nd' 'extensor_ulnaris_cerpi_nd' 'pectoralis_major_d' 'deltoidus_d' 'infraspinatus_d' 'biceps_brachialis_d' 'triceps_brachialis_d' 'flexor_ulnaris_cerpi_d' 'extensor_ulnaris_cerpi_d'}

writetable(resultats_bis,string(filename(1:10))+'_results_iEMG');

%% step 8 - EMG envelope processing ---

% EMG-enveloppe (concatenating up-swing + down-swing)

for m=2:length(val_EMG_norm(1,:))
for i=1:numel(cycle_ascendant)
    cycle_EmgEnv_standard.(['muscle' num2str(m)]).(['cycle' num2str(i)])=val_EMG_norm(P0_cycle_preparation_(i):P1_cycle_descendant_(i),m);
end
end

% creation vector time

i=1;
n=1;
for m=2:length(val_EMG_norm(1,:))
for i=1:numel(cycle)
            Nb_pts=length(cycle_EmgEnv_standard.(['muscle' num2str(m)]).(['cycle' num2str(i)]));
            Nb_sec=Nb_pts/1000;
            t_EmgEnv_bis.(['cycle' num2str(i)])(:,1)=linspace(0,Nb_sec,Nb_pts);
    end
end

% normalize time in % of cycle 

i=1;
n=1;
for i=1:numel(cycle)
        for n=1:length(t_EmgEnv_bis.(['cycle' num2str(i)]))
            t_EmgEnv_norm.(['cycle' num2str(i)])(n,1)=(t_EmgEnv_bis.(['cycle' num2str(i)])(n,1)*100)/(max(t_EmgEnv_bis.(['cycle' num2str(i)])));
        end
end

% plot all cycle to visualize (remove cycle errors)
i=1;
m=1;
for m=2:length(val_EMG_norm(1,:))
    figure('name','EMG_env')
for i=1:numel(cycle)
plot(t_EmgEnv_norm.(['cycle' num2str(i)]),cycle_EmgEnv_standard.(['muscle' num2str(m)]).(['cycle' num2str(i)]))
hold on 
end
hold off
end

%% save all variables to analyze the kinematic files 

save(string(filename(1:10))+'all_variables.mat')

%% step 9 - save EMG envelope ---

cd(dirname_save)

% save in .mat
save(string(filename(1:10))+'_t_EMG_env_norm.mat','t_EmgEnv_norm');
save(string(filename(1:10))+'_cycle_EmgEnv.mat','cycle_EmgEnv_standard');

clearvars -except dirname_save dirname dirname_marker dirname_max_EMG k list list_marker

end
