clear all;clc;close all;

% load 20 runs of one subject's surface data and then average across all
% stimuli repeats

% for each vertex pick the direction with the largest response
% plot each vertex's response to each direction in separate maps 
% generate 18 .mgz files

restoredefaultpath
addpath(genpath('~/Documents/GitHub/gifti')); % https://github.com/gllmflndn/gifti
addpath(genpath('~/Documents/GitHub/fMRI-Matlab')); % https://github.com/HaiyangJin/fMRI-Matlab

BASE = '/Users/pw1246/Desktop/MRI/Decoding/';
addpath(genpath(BASE));

fs_setup('/Applications/freesurfer/7.2.0',[],1);
fs_subjdir([BASE 'derivatives/freesurfer']);

%%
sub = 'sub-0228';
ses = {'01','02'};
run = [1:10]';

[L_samples R_samples stim_label] = load_surf(BASE,sub,ses,run);
%%

L_hemi = zeros(8,size(L_samples,2));
R_hemi = zeros(8,size(R_samples,2));
for dir = 1:8
L_hemi(dir,:) = abs(mean(L_samples(stim_label == dir,:)));
R_hemi(dir,:) = abs(mean(R_samples(stim_label == dir,:)));
end

[L_m L_i] = max(L_hemi);
nanidx = isnan(L_m);
L_i(nanidx) = 0;
[R_m R_i] = max(R_hemi);
nanidx = isnan(R_m);
R_i(nanidx) = 0;

%% max direction map
fs_savemgz(sub, L_i','l.max.mgz', [pwd,'/maps'], 'lh');
fs_savemgz(sub, R_i', 'r.max.mgz', [pwd,'/maps'], 'rh');

%% response map to each direction
fs_savemgz(sub, L_hemi(1,:)','l.1.mgz', [pwd,'/maps'], 'lh');
fs_savemgz(sub, R_hemi(1,:)', 'r.1.mgz', [pwd,'/maps'], 'rh');

%%
fv_surf