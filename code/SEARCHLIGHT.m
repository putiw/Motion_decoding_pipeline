clear all;clc;close all;

restoredefaultpath
addpath(genpath('~/Documents/GitHub/gifti')); % https://github.com/gllmflndn/gifti
addpath(genpath('~/Documents/GitHub/fMRI-Matlab')); % https://github.com/HaiyangJin/fMRI-Matlab

% BASE = '/Users/pw1246/Desktop/MRI/Decoding/';
BASE =  '/Users/pw1246/Desktop/MRI/Decoding/';
addpath(genpath(BASE));

fs_setup('/Applications/freesurfer/7.2.0',[],1);
fs_subjdir([BASE 'derivatives/freesurfer']);
subjdir = [BASE 'derivatives/freesurfer'];

%%
subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
for sub_idx = 1:9

sub = subject{sub_idx};
ses = {'01','02'}
datafile = [pwd '/dalma_data/' sub '-ses-' ses{:} '.mat'];
load(datafile)

%%

% subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
% for sub_idx = 1
% 
% sub = subject{sub_idx};
% 
% 
% ses =   {'01','02'};
% run = [1:10]';
% 
% 
% inflate = fullfile('/Users/pw1246/Desktop/MRI/Decoding/derivatives/freesurfer/fsaverage6/surf/lh.inflated');
% inflate_file = fs_read_surf(inflate);
% [L_samples R_samples stim_label] = load_surf(BASE,sub,ses,run,'fsaverage6');
% L_label = stim_label; R_label = stim_label;
% R_label(find(stim_label<=4.5)) = stim_label(find(stim_label<=4.5))+4;
% R_label(find(stim_label>=4.5)) = stim_label(find(stim_label>=4.5))-4;
% 
% v_n = size(L_samples,2);
% % find what is inside for each searchlight
% 
% SearchLight = zeros(round(v_n),100);
% for ii = 1:round(v_n)
% center = inflate_file.coord(:,ii);
% distance = sqrt(sum((inflate_file.coord - center).^ 2));
% [~,SearchLight(ii,:)] = mink(distance,100);
% end
%  R_inflate = fullfile('/Users/pw1246/Desktop/MRI/Decoding/derivatives/freesurfer/fsaverage6/surf/rh.inflated');
% 
% R_inflate_file = fs_read_surf(R_inflate);
% R_v_n = size(R_samples,2);
% % find what is inside for each searchlight
% 
% R_SearchLight = zeros(round(v_n),100);
% for ii = 1:round(R_v_n)
% center = R_inflate_file.coord(:,ii);
% distance = sqrt(sum((R_inflate_file.coord - center).^ 2));
% [~,R_SearchLight(ii,:)] = mink(distance,100);
% end
% savedata = fullfile(pwd,'dalma_data',[sub '-ses-' ses{:} '.mat']);
% save(savedata,'R_label','L_label','SearchLight','L_samples','R_samples','stim_label','R_SearchLight','R_v_n');
% end
%%

% C = intersect(SearchLight(:,1),SearchLight(5000,:));
% %% visual
% close all
% downsample = 20;
% idx = find(inflate_file.coord(2,:)<0);
% idx = idx(:,1:round(end/downsample));
% f1 = figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
% s1 = scatter3(inflate_file.coord(1,idx), ...
%     inflate_file.coord(2,idx),inflate_file.coord(3,idx),30,'.');
% axis equal
% view(-90,0)
% xlabel('x'); xlim([min(inflate_file.coord(1,:)) max(inflate_file.coord(1,:))]);
% ylabel('y'); ylim([min(inflate_file.coord(2,:)) max(inflate_file.coord(2,:))]);
% zlabel('z'); zlim([min(inflate_file.coord(3,:)) max(inflate_file.coord(3,:))]);
% %%
% hold on
% for kk = 1:10
%    randidx = randi(numel(idx));
% scatter3(inflate_file.coord(1,SearchLight(idx(randidx),1)), ...
%     inflate_file.coord(2,SearchLight(idx(randidx),1)),...
%     inflate_file.coord(3,SearchLight(idx(randidx),1)),200,'r.');
% 
% scatter3(inflate_file.coord(1,SearchLight(idx(randidx),:)), ...
%     inflate_file.coord(2,SearchLight(idx(randidx),:)),...
%     inflate_file.coord(3,SearchLight(idx(randidx),:)),80,'.','MarkerEdgeAlpha',.5);
% end
%%

%
%

%% MATLAB Classify
tic
code_rep = 300; % bootstraps repeats
test_n = 1; % percetange of data set as testing trials

trial_n = 160;  % number of trials 
dir_n = 8; % number of directions
nScans = trial_n/dir_n; % number of scans

Training_group = sort(repmat(repelem(1:8,1)',nScans-test_n,1)); % training dataset stimulus label
Testing_group = repmat(repelem(1:8,1)',test_n,1); % testing dataset stimulus label 

pc = zeros(code_rep,size(SearchLight,1)); % percentage correct
parfor ii = 1:size(SearchLight,1) 
    data = L_samples(:,SearchLight(ii,:));

    [n,voxel_n] = size(data); % trial count and voxel size
        
    for kk = 1:code_rep
        
        Training_data = [];
        Testing_data = [];     
        
        for mm = 1:dir_n    
            train_run = 1:nScans;
            test_run = randperm(nScans,test_n); %random select (test_n) trial for each direction as testing trial index (without replacement)
            train_run(test_run) = []; % set the rest as the training trial index          
            dir_idx = find(R_label==mm); % find trial index for this direction          
            Testing_data = [Testing_data; data(dir_idx(test_run),:)];
            Training_data = [Training_data; data(dir_idx(train_run),:)];            
        end
     
        decoded_group = classify(Testing_data,Training_data,Training_group,'diaglinear');                
        pc(kk,ii) = length(find(decoded_group == Testing_group))/length(Testing_group); % percentage correct          
    end
end

L_result = mean(pc)';
filesavenamel = ['r.' sub '_' ses{:} '.mgz'];

 fs_savemgz(sub, L_result',filesavenamel, [pwd,'/maps'], 'rh');

%% Right



R_pc = zeros(code_rep,R_v_n); % percentage correct
parfor ii = 1:R_v_n
    data = R_samples(:,R_SearchLight(ii,:));

    [n,voxel_n] = size(data); % trial count and voxel size
        
    for kk = 1:code_rep
        
        Training_data = [];
        Testing_data = [];     
        
        for mm = 1:dir_n    
            train_run = 1:nScans;
            test_run = randperm(nScans,test_n); %random select (test_n) trial for each direction as testing trial index (without replacement)
            train_run(test_run) = []; % set the rest as the training trial index          
            dir_idx = find(R_label==mm); % find trial index for this direction          
            Testing_data = [Testing_data; data(dir_idx(test_run),:)];
            Training_data = [Training_data; data(dir_idx(train_run),:)];            
        end
     
        decoded_group = classify(Testing_data,Training_data,Training_group,'diaglinear');                
        R_pc(kk,ii) = length(find(decoded_group == Testing_group))/length(Testing_group); % percentage correct          
    end
end

toc

R_result = mean(R_pc)';
filesavenamer = ['r.' sub '_' ses{:} '.mgz'];

 fs_savemgz(sub, R_result',filesavenamer, [pwd,'/maps'], 'rh');
end