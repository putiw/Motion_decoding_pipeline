clear all
close all
clc
restoredefaultpath
addpath(genpath('~/Documents/GitHub/gifti')); % https://github.com/gllmflndn/gifti
addpath(genpath('~/Documents/GitHub/fMRI-Matlab')); % https://github.com/HaiyangJin/fMRI-Matlab
BASE = '/Users/pw1246/Desktop/MRI/Decoding/';

addpath(genpath(BASE));

fs_setup('/Applications/freesurfer/7.2.0',[],1);
fs_subjdir([BASE 'derivatives/freesurfer']);
%%
filename_L_V = fullfile([pwd '/maps/l.fsaverage6_vertical.mgz']);
filename_R_V = fullfile([pwd '/maps/r.fsaverage6_vertical.mgz']);
filename_L_H = fullfile([pwd '/maps/l.fsaverage6_horizontal.mgz']);
filename_R_H = fullfile([pwd '/maps/r.fsaverage6_horizontal.mgz']);
aaRH = fs_readmgh(filename_R_H);aaRV = fs_readmgh(filename_R_V);
aaLH = fs_readmgh(filename_L_H);aaLV = fs_readmgh(filename_L_V);
aaRH=aaRH.vol;aaRV=aaRV.vol;aaLH=aaLH.vol;aaLV=aaLV.vol;
diffL = aaLH-aaLV;
diffR = aaRH-aaRV;
%%
idx = kmeans(diffL, 3); % 3 number of classes
scatter3(diffL(1,:), diffL(2,:), diffL(3,:), 15, idx, 'filled');
%% t test
Map_H_R = zeros(40962,9);
Map_H_L = zeros(40962,9);

subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
for sub_idx = 1:9
sub = subject{sub_idx};
ses = {'01','02'};
filename_R = fullfile([pwd '/maps/r.' sub '_' ses{:} '.mgz']);
filename_L = fullfile([pwd '/maps/l.' sub '_' ses{:} '.mgz']);
aaR = fs_readmgh(filename_R);
aaL = fs_readmgh(filename_L);
Map_H_R(:,sub_idx) = aaR.vol;
Map_H_L(:,sub_idx) = aaL.vol;
end
%%
%fv_surf
Map_R = zeros(40962,9);
Map_L = zeros(40962,9);
Map_R_H = zeros(40962,9);
Map_L_H = zeros(40962,9);
Map_R_V = zeros(40962,9);
Map_L_V = zeros(40962,9);
subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
for sub_idx = 1:9
sub = subject{sub_idx};
ses = {'01','02'};
ses_V = {'03','04'};
filename_R = fullfile([pwd '/maps/r.' sub '_' ses{:} '.mgz']);
filename_L = fullfile([pwd '/maps/l.' sub '_' ses{:} '.mgz']);
filename_R_V = fullfile([pwd '/maps/r.' sub '_' ses_V{:} '.mgz']);
filename_L_V = fullfile([pwd '/maps/l.' sub '_' ses_V{:} '.mgz']);
aaR = fs_readmgh(filename_R);
aaL = fs_readmgh(filename_L);
aaR_V = fs_readmgh(filename_R_V);
aaL_V = fs_readmgh(filename_L_V);
Map_R_H(:,sub_idx) = aaR.vol;
Map_L_H(:,sub_idx) = aaL.vol;
Map_R_V(:,sub_idx) = aaR_V.vol;
Map_L_V(:,sub_idx) = aaL_V.vol;
Map_R(:,sub_idx) = aaR.vol-aaR_V.vol; 
Map_L(:,sub_idx) = aaL.vol-aaL_V.vol; 

end
Map_R = mean(Map_R,2);Map_L = mean(Map_L,2);

filename_L_t = ['l.H_V.mgz'];
fs_savemgz(sub, Map_L, filename_L_t, [pwd,'/maps/'], 'lh');
filename_R_t = ['r.H_V.mgz'];
fs_savemgz(sub, Map_R, filename_R_t, [pwd,'/maps/'], 'rh');
%% get chance diff
Map_R_H = mean(Map_R_H,2);Map_L_H = mean(Map_L_H,2);
Map_R_V = mean(Map_R_V,2);Map_L_V = mean(Map_L_V,2);
%%
chance = zeros(2000,1);
all = [Map_R_H;Map_R_V;Map_R_H;Map_R_V];
for ii = 1:5000
 chance(ii,1) = all(randi(40962*4))-all(randi(40962*4));
end
   
prctile(chance(:),95)

bar_H = [Map_R_H(Map_R>=0.0383);Map_L_H(Map_L>=0.0383)];
bar_V = [Map_R_V(Map_R>=0.0383);Map_L_V(Map_L>=0.0383)];
%%
[h,p,ci,stats_R] = ttest2(Map_R_H',Map_R_V','Tail','right');
[h,p,ci,stats_L] = ttest2(Map_L_H',Map_L_V','Tail','right');
Map_R_t = zeros(40962,1);
Map_L_t = zeros(40962,1);
inflate = fullfile('/Users/pw1246/Desktop/MRI/Decoding/derivatives/freesurfer/fsaverage6/surf/lh.inflated');
inflate_file = fs_read_surf(inflate);
inflate_R = fullfile('/Users/pw1246/Desktop/MRI/Decoding/derivatives/freesurfer/fsaverage6/surf/rh.inflated');
inflate_file_R = fs_read_surf(inflate_R);
SearchLight = zeros(40962,100);
SearchLight_R = zeros(40962,100);

for ii = 1:40962
center = inflate_file.coord(:,ii);
distance = sqrt(sum((inflate_file.coord - center).^ 2));
[~,SearchLight(ii,:)] = mink(distance,100);
centerR = inflate_file_R.coord(:,ii);
distanceR = sqrt(sum((inflate_file_R.coord - centerR).^ 2));
[~,SearchLight_R(ii,:)] = mink(distanceR,100);
Map_R_t(ii,:) = sum(stats_R.tstat(SearchLight_R(ii,1:30))>=1.86);
Map_L_t(ii,:) = sum(stats_L.tstat(SearchLight(ii,1:30))>=1.86);
end
%% random data to find t cluster threshold

kkn = 100;
cluster_R = zeros(kkn,40962);
cluster_L = zeros(kkn,40962);
for kk = 1:kkn
ALL_R = [Map_R_H;Map_R_V]'; 
ALL_R = ALL_R(:,randperm(81924));
[~,~,~,temp_R] = ttest2(ALL_R(:,1:40962),ALL_R(:,40963:end));
ALL_L = [Map_L_H;Map_L_V]'; 
ALL_L = ALL_L(:,randperm(81924));
[~,~,~,temp_L] = ttest2(ALL_L(:,1:40962),ALL_L(:,40963:end));
nanidx = isnan(temp_R.tstat);
temp_R.tstat(nanidx) = 0;
nanidx = isnan(temp_L.tstat);
temp_L.tstat(nanidx) = 0;
for ii = 1:40962
cluster_R(kk,ii)=sum(temp_R.tstat(SearchLight_R(ii,1:30))>=1.86);
cluster_L(kk,ii)=sum(temp_L.tstat(SearchLight(ii,1:30))>=1.86);
end
end

%%
all_cluster = [cluster_R,cluster_L];
finalsize = prctile(all_cluster(:),95);

stats_R.tstat = stats_R.tstat';
stats_L.tstat = stats_L.tstat';
stats_R.tstat(Map_R_t<=finalsize)=0;
stats_L.tstat(Map_L_t<=finalsize)=0;

filename_L_t = ['l.t17.mgz'];
fs_savemgz(sub, stats_L.tstat, filename_L_t, [pwd,'/maps/'], 'lh');
filename_R_t = ['r.t17.mgz'];
fs_savemgz(sub, stats_R.tstat, filename_R_t, [pwd,'/maps/'], 'rh');

%%


%%
subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
for sub_idx = 9
sub = subject{sub_idx};
ses = {'01','02'};
load(fullfile([pwd '/maps/' sub '-ses-' ses{:} '-L.mat']));
% load(fullfile([pwd '/maps/' sub '-ses-' ses{:} '-R.mat']));
% filename_R = ['r.' sub '_' ses{:} '.mgz'];
filename_L = ['l.' sub '_' ses{:} '.mgz'];
fs_savemgz(sub, L_result, filename_L, [pwd,'/maps/'], 'lh');
% fs_savemgz(sub, R_result, filename_R, [pwd,'/maps/'], 'rh');
end

