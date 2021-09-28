% Setup and run classification and/or decoding algorithm

clear all
close all
clc
tic

% Dependencies
restoredefaultpath
addpath(genpath('~/Documents/GitHub/TAFKAP')); % https://github.com/jeheelab/TAFKAP
% Toolboxes for single trial BOLD amplitude estimation
addpath(genpath('~/Documents/GitHub/GLMsingle')); % https://github.com/kendrickkay/GLMsingle
addpath(genpath('~/Documents/GitHub/fracridge')); % https://github.com/nrdg/fracridge

addpath(genpath('~/Documents/GitHub/rokers_mri_lab/code/invChol'));

%BASE = '~/Desktop/motion/';
BASE = '~/Dropbox (RVL)/MRI/Decoding/';
% BASE = '/Volumes/Macintosh HD/Decoding/';
addpath(genpath(BASE));

% Figure defaults
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize', 14)

% Set up parameters
sub = 'sub-0204';
ses =  {'01','02'}; % {'03','04'};
run = [1:10]';

%roi =  {'V1','hMT'};
% roi =  {'V1','V2','hMT','IPS0'};
%roi =  {'V1','V2','V3','V3A','hV4','LO','hMT','MST','IPS'};
roi = {'V1','V2','V3','V3A','V3B','hV4','LO1','LO2','hMT','MST','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','VO1','VO2','SPL1','PHC1','PHC2','FEF'};
params = SetupTAFKAP(); 
[samples,stim_label] = loadmydata(BASE,sub,ses,run,roi,params);

%
voxelsize = numel(roi,1);
for i = 1:numel(roi)
    voxelsize(i,1)=size(samples{i},2);
end
savedata = fullfile(pwd,'data',[sub '-ses-' ses{:} '.mat']);
save(savedata,'samples','stim_label','roi','voxelsize');


%% Estimate single trial BOLD amplitude 

switch params.sample_unit
    case {'raw'}
        % Requires data in cell format: subject x run
        % Not by ROI since you need a noise pool
        design = {};
        for p=1:10
            design{p} = zeros(250,7);
            if mod(p,2)==1
                cnt = 1; curcond = 8;
                while cnt <= 250
                    if curcond ~= 1
                        design{p}(cnt,curcond-1) = 1;
                    end
                    cnt = cnt + 2;
                    curcond = mod2(curcond-1,8);
                end
            else
                cnt = 1; curcond = 1;
                while cnt <= 250
                    if curcond ~= 1
                        design{p}(cnt,curcond-1) = 1;
                    end
                    cnt = cnt + 2;
                    curcond = mod2(curcond+1,8);
                end
            end
        end
        
        % Copy design nSessions times
        
        % Drop first 10 frames
        design = cellfun(@(x) x(11:end,:),design,'UniformOutput',0);
        samples = cellfun(@(x) x(:,:,:,11:end),samples,'UniformOutput',0);
        results = GLMestimatesingletrial(design,samples,3,1.5,'testglmsingle',struct());
end

%% Run TAFKAP

% Setup design (parameters)
nScans = length(ses)*length(run); % scans per subject
nFolds = 8; % Set to multiple of number of processing cores
nDirs = 8; % motion directions
% params.stimval = 22.5.*reshape(repmat([5:-1:1 8:-1:6 4:8 1:3],1,nScans/2),1,[])'; % stimulus labels
%params.stimval = reshape(repmat([5:-1:1 8:-1:6 4:8 1:3],1,nScans/2),1,[])'; % categorical stimulus labels
params.stimval = stim_label;

% TODO: Move to SetupTAFKAP
switch params.sample_unit
    case {'scan'}
        params.runNs = reshape(repmat(1:nScans,nDirs,1),1,[])'; %trainruns; % stimulus block/run
    case {'trial'}
        nRepeats = 15;
        params.runNs = reshape(repmat(1:nScans,nRepeats*nDirs,1),1,[])'; %trainruns; % stimulus block/run
    case {'raw'}
        % do nothing
    otherwise
        error('Unknown sample unit')
end
        
pre = cell(nFolds,1);  % Preallocate
for ii = 1:nFolds
    p{ii} = params;
    c = cvpartition(params.stimval, 'Holdout', 0.1); % stratify by motion direction, but not scan
    p{ii}.train_trials = c.training;
    p{ii}.test_trials = c.test;
    pre{ii} = params.stimval(c.test);
end

% Run TAFKAP
ests = cell(numel(roi),1); % Preallocate
uncs = cell(numel(roi),1); % Preallocate
pres = cell(numel(roi),1); % Preallocate
for whichRoi = 1:numel(roi)
    parfor ii = 1:nFolds
        rng(ii);% To counter the effects of TAFKAP_Decode setting the system rand seed to const. This was an EXTRAORDINARLY hard bug to find.
        [est{ii}, unc{ii}, liks{ii}, hypers{ii}] = TAFKAP_Decode(samples{whichRoi}, p{ii});
    end
    ests{whichRoi} = cell2mat(est');
    uncs{whichRoi} = cell2mat(unc');
    pres{whichRoi} = cell2mat(pre);
end

%% Calculate performance

params.subjects = sub;
% plot_results(user, params, rois, acc) % Plot results
saveresult = cell(numel(roi),1);
for whichRoi = 1:numel(roi)
    % Plot confusion matrix
    % conmat = confusionmat(categorical(ceil(pres{rr}/22.5)),categorical(ceil(ests{rr}/22.5))); % continuous data
    conmat = confusionmat(categorical(pres{whichRoi}),categorical(ests{whichRoi})); % categorical data
    conmat = 100.*conmat./sum(conmat,2); % convert to % accuracy
    
    figure(1); hold on;
    sqrtRois = ceil(sqrt(numel(roi)));
    subplot(sqrtRois,sqrtRois,whichRoi)
    conmat = [conmat; conmat(1,:)]; % wrap matrix
    conmat = [conmat, conmat(:,1)];
    % imagesc(conmat');
    clim = [0 100]; % upper, lower limits
    imagesc(conmat', clim);
    saveresult{whichRoi} = conmat';
    title(roi{whichRoi})
    xlabel('Presented direction')
    ylabel('Decoded direction')
    
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    yticks([1:9])
    yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))    

    set(gca,'YDir','normal')
    
    if whichRoi == length(roi)
        cb = colorbar;
        cb.Label.String = 'Classification performance (%)';
    end
    
    disp(['Classification performance ' roi{whichRoi} ': ' num2str(100.*mean(pres{whichRoi}==ests{whichRoi})) '%'])
end


 f = fullfile(pwd,'result',[sub '-ses-' ses{:} '-TAFKAP.mat']);
 save(f,'saveresult','ests','uncs','pres','roi','voxelsize');


% hp4 = get(subplot(2,2,4),'Position');
% h = colorbar('Position', [hp4(1)+hp4(3)+0.02  hp4(2)  0.03  hp4(2)+hp4(3)*2.1]);
% h.Label.String = 'Classification performance (%)';
% axis tight
    
matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual')
set(groot,'defaultFigurePaperPositionMode','manual')
saveas(gcf, ['../figures/Confusion_matrix-' datestr(now,30) '.pdf'])


% Plot uncertainty as a function of motion direction
figure
for whichRoi = 1:numel(roi)
    subplot(sqrtRois,sqrtRois,whichRoi)
    hold on
    title(roi{whichRoi})
    scatter(pres{whichRoi},uncs{whichRoi})
    y = splitapply(@mean,uncs{whichRoi},pres{whichRoi});
    plot(1:8,y)
    
    xlim([.5 8.5])
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    
    xlabel('Presented motion direction')
    ylabel('Entropy (nats)')
end
saveas(gcf, ['../figures/Uncertainty-' datestr(now,30) '.pdf'])


%% MATLAB Classify
% 
% tic
% code_rep = 1000; % bootstraps repeats
% voxelmax = 1e5; % maximum voxels to use per ROI
% test_n = 1; % percetange of data set as testing trials
% 
% trial_n = size(samples{1,1},1);  % number of trials 
% dir_n = 8; % number of directions
% nScans = trial_n/dir_n; % number of scans
% result = zeros(numel(roi),4); % allocate results
% cnfm = cell(numel(roi),1); % each runs confusion matrix is restored separately
% 
% Training_group = sort(repmat(repelem(1:8,1)',nScans-test_n,1)); % training dataset stimulus label
% Testing_group = repmat(repelem(1:8,1)',test_n,1); % testing dataset stimulus label 
% random_Training_group = Training_group(randperm(length(Training_group))); % random shuffle stimulus label for decoding chance
% 
% pval = zeros(numel(roi),1); % p value (significance from chance)
% pc = zeros(code_rep,numel(roi)); % percentage correct
% chancepc = zeros(code_rep,numel(roi)); % percentage correct from chance
% 
% 
% for ii = 1:numel(roi) % loop over ROIs
%     
%     data = samples{1,ii}; % get data for this ROI    
%     [n,voxel_n] = size(data); % trial count and voxel size
%     
%     cnfm{ii,1} = zeros(dir_n,dir_n); % setting up empty confusion matrix before bootstraps
%     
%     for kk = 1:code_rep
%         
%         if voxel_n > voxelmax
%             data = data(:,randperm(voxel_n,voxelmax),:); % randomly select (voxelmax) amount of voxels if there is a upper limit for voxel counts
%         end
% 
%         Training_data = [];
%         Testing_data = [];     
%         
%         for mm = 1:dir_n    
%             train_run = 1:nScans;
%             test_run = randperm(nScans,test_n); %random select (test_n) trial for each direction as testing trial index (without replacement)
%             train_run(test_run) = []; % set the rest as the training trial index          
%             dir_idx = find(stim_label==mm); % find trial index for this direction          
%             Testing_data = [Testing_data; data(dir_idx(test_run),:)];
%             Training_data = [Training_data; data(dir_idx(train_run),:)];            
%         end
%     
%         decoded_group = classify(Testing_data,Training_data,Training_group,'diaglinear');
%         random_decoded_group = classify(Testing_data,Training_data,random_Training_group,'diaglinear');
%                 
%         pc(kk,ii) = length(find(decoded_group == Testing_group))/length(Testing_group); % percentage correct
%         chancepc(kk,ii) = length(find(random_decoded_group == Testing_group))/length(Testing_group); % percentage correct from chance
%         
%         temp_cnfm = zeros(dir_n,dir_n); % confusion matrix for this bootstrap       
%         for de = 1:dir_n
%             for te = 1:dir_n
%                 temp_cnfm(de,te) = length(find(decoded_group == de & Testing_group == te))/length(find(Testing_group == te));
%             end
%         end
%         cnfm{ii,1} = cnfm{ii,1} + temp_cnfm; 
%     end
%     
%     cnfm{ii,1} = cnfm{ii,1}./code_rep; % average confusion matrix after bootstraps
%     
%     [~,pval(ii),~,~] = ztest(pc(:,ii),mean(chancepc(:,ii)),std(chancepc(:,ii)),'Tail','Right'); % p value from chance
%        
%     SEM = std(pc(:,ii))/sqrt(length(pc(:,ii)));    
%     result(ii,1) = mean(pc(:,ii)).*100;
%     result(ii,2) = SEM.*196*2;   
%     SEM_chance = std(chancepc(:,ii))/sqrt(length(chancepc(:,ii)));
%     result(ii,3) = mean(chancepc(:,ii)).*100;
%     result(ii,4) = SEM_chance.*196*2;
%     voxelsize(ii)=size(data,2);
% end
% 
% toc
% %roi comparison stats (p value, confidence interval)
% CIroi_mean = zeros(numel(roi),numel(roi));
% CIroi = zeros(numel(roi),numel(roi));
% pvalueroi = zeros(numel(roi),numel(roi));
% for ii = 1:numel(roi)
%     for jj = 1:numel(roi)        
%         [~,pvalueroi(ii,jj)] = ttest(pc(:,ii),pc(:,jj));
%         CIroi_mean(ii,jj) = mean(pc(:,ii)-pc(:,jj));
%         CIroi(ii,jj) = (std(pc(:,ii)-pc(:,jj))/sqrt(code_rep))*1.96;      
%     end
% end
% 
% %% save data in the result folder
% f = fullfile(pwd,'result',[sub '-ses-' ses{:} '-classify.mat']);
% save(f,'cnfm','result','pvalueroi','pval','CIroi_mean','CIroi','roi','voxelsize');
% 
% %% plot confusion matrix
% %close all
% %whichroi = [1 4 9 11];
% whichroi = 1:size(cnfm,1);
% plot_confusionM(cnfm,roi,voxelsize,whichroi);
