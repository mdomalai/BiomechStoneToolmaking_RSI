%% treatment Data (kinematic) Hometech project 
% code 2 : Pre-processed kinematic
% author : Robin Macchi 

% clear workspace...
clear all
close all
clc 

%% import data 
prompt = 'Select the patheway of the kinematic data'
dirname_kinematic=uigetdir % Where are your files (kinematic) located? 
cd(dirname_kinematic)
list=dir(dirname_kinematic);
list=list(3:end);

for k=1:length(list)
cd(dirname_kinematic)
list=dir(dirname_kinematic);
list=list(3:end);
filename=list(k).name;
opts = delimitedTextImportOptions("NumVariables", 15);
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["time", "shoulderflexd", "shoulderaddd", "shoulderrotd", "elbowflexd", "forearmprond", "wristflexd", "wristaddd", "shoulnderflexnd", "shoulnderandndnd", "shoulnderrotnd", "elbowflexnd", "forearmpronnd", "wristflexnd", "wristandndnd"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
data_kinematic = readtable(filename, opts);
data_kinematic = table2array(data_kinematic);
clear opts

% open previous EMG treatment for the phase detection (up-swing and
% down-swing)
if k==1
    prompt = 'Select the patheway of the backup EMG files'
dirname_save=uigetdir % Where are your results (EMG) located? 
else
end
cd(dirname_save)
load(string(filename(1:10))+'all_variables.mat')

%% reduced file size to synchronize with EMG (load the variables of the previous EMG treatment)

kinematic=data_kinematic([deb1:fin1],:);
t2=data_kinematic(:,1);

% detection of up-swing and down-swing phases (vector time)
n=1;
for i=1:length(cycle)
    cycle_kinematic_time.(['number' num2str(n)])=t2(R_MC1_pos_deb(i):p_R_MC1_min(i));
    n=n+1;
end

% detection of up-swing and down-swing phases (joint-angles data)

n=1;
for l=1:length(kinematic(1,:))
for i=1:length(cycle)
    cycle_kinematic.(['angle' num2str(l)]).(['number' num2str(i)])=kinematic(R_MC1_pos_deb(i):p_R_MC1_min(i),l);
end
end

% new vector time (starting with 0)

i=1;
n=1;
for i=1:length(cycle)
            Nb_pts=length((cycle_kinematic.(['angle' num2str(2)]).(['number' num2str(i)])(:,1)));
            Nb_sec=Nb_pts/300;  % acquisition frequency = 300 Hz
            t_kinematic_bis.(['cycle' num2str(i)])(:,1)=linspace(0,Nb_sec,Nb_pts);
end

% normalization of time in %

i=1;
n=1;
for i=1:length(cycle)
        for n=1:length(t_kinematic_bis.(['cycle' num2str(i)]))
            t_kinematic_norm.(['cycle' num2str(i)])(n,1)=(t_kinematic_bis.(['cycle' num2str(i)])(n,1)/(max(t_kinematic_bis.(['cycle' num2str(i)]))))*100;
        end
end

% matrix with joint-angle names

name_angle.(['angle' num2str(2)])={'F-Es_r'};
name_angle.(['angle' num2str(3)])={'Ab-Ads_r'};
name_angle.(['angle' num2str(4)])={'Rots_r'};
name_angle.(['angle' num2str(5)])={'F-Ee_r'};
name_angle.(['angle' num2str(6)])={'P-Se_r'};
name_angle.(['angle' num2str(7)])={'F-Ew_r'};
name_angle.(['angle' num2str(8)])={'Ab-Adw_r'};
name_angle.(['angle' num2str(9)])={'F-Es_l'};
name_angle.(['angle' num2str(10)])={'Ab-Ads_l'};
name_angle.(['angle' num2str(11)])={'Rots_l'};
name_angle.(['angle' num2str(12)])={'F-Ee_l'};
name_angle.(['angle' num2str(13)])={'P-Se_l'};
name_angle.(['angle' num2str(14)])={'F-Ew_l'};
name_angle.(['angle' num2str(15)])={'Ab-Adw_l'};

% find the frame for each %

z=1:1:100;
z=z';

for i=1:length(cycle)
    n=1;
for c=1:100
   P_deb=find(t_kinematic_norm.(['cycle' num2str(i)])>=z(c));
   P_c.(['cycle' num2str(i)])(n,1)=P_deb(1);
   n=n+1;
   clear P_deb
end
end

% averaged kinematic data every 1%

for l=1:1:length(kinematic(1,:))
for i=1:length(cycle)
        n=2;
        for k=1:length(P_c.(['cycle' num2str(i)]))-1
cycle_kinematic_moy.(['angle' num2str(l)]).(['cycle' num2str(i)])(1,1)=mean(cycle_kinematic.(['angle' num2str(l)]).(['number' num2str(i)])(1:P_c.(['cycle' num2str(i)])(k),1));
cycle_kinematic_moy.(['angle' num2str(l)]).(['cycle' num2str(i)])(n,1)=mean(cycle_kinematic.(['angle' num2str(l)]).(['number' num2str(i)])(P_c.(['cycle' num2str(i)])(k):P_c.(['cycle' num2str(i)])(k+1),1));
    n=n+1;
        end
end
end
% Need to manually remove the cycle errors 

% save matlab file

if k==1
    prompt = 'Select the patheway of your backup files'
dirname_save=uigetdir
else
end
cd(dirname_save)
save(string(filename(1:10))+'_cycle_kinematic.mat','cycle_kinematic_moy');


% computation of mean cycle and standard deviation

for l=1:1:length(kinematic(1,:))
for y=1:100
                 n=1;
    for i=1:23
        passage_visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(n,y)=cycle_kinematic_moy.(['angle' num2str(l)]).(['cycle' num2str(i)])(y,1);
    n=n+1;
        end
             end
end

for l=1:1:length(kinematic(1,:))
y=1;
    for n=1:100
visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1)=mean(passage_visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(:,n));
error_visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1)=std(passage_visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(:,n));
y=y+1;
    end
end

% create vector of standard deviation

for l=1:1:length(kinematic(1,:))
 y=1
    for n=1:100
       visu_sd_sup_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1)=visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1)+error_visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1);
       visu_sd_inf_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1)=visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1)-error_visu_1.(['angle' num2str(l)]).(['cycle' num2str(l)])(y,1);
       y=y+1; 
    end
end

% plot averaged cycles (after manually remove cycle errors)
for j=2:length(kinematic(1,:))
n=100;
t_visu=1:1:100;
t_visu=t_visu';
figure(j)
plot(t_visu,visu_sd_sup_1.(['angle' num2str(j)]).(['cycle' num2str(j)]))
hold on
plot(t_visu,visu_sd_inf_1.(['angle' num2str(j)]).(['cycle' num2str(j)]))
h=area(z,[visu_sd_sup_1.(['angle' num2str(j)]).(['cycle' num2str(j)]) visu_sd_inf_1.(['angle' num2str(j)]).(['cycle' num2str(j)])-visu_sd_sup_1.(['angle' num2str(j)]).(['cycle' num2str(j)])],'LineStyle',':')
z=1:n;
h(1).FaceColor = [1 1 1];
h(2).FaceColor = [0.8 0.8 0.8];
hold on
plot(t_visu,visu_1.(['angle' num2str(j)]).(['cycle' num2str(j)]),'linewidth',1,'color','k')
ylabel(string(name_angle.(['angle' num2str(j)]))+' (°)','FontSize',10,...              
       'FontName','Times New Roman','Color','k')
end

clearvars -except dirname_kinematic dirname_save

pause 
prompt = 'Click enter to continue the code when you have finished checking the curves'

close all

end









