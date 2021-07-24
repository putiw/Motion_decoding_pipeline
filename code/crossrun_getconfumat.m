% Setup and run classification and/or decoding algorithm

clear all
close all
clc
tic

% Dependencies
restoredefaultpath
addpath(genpath('~/Documents/GitHub/TAFKAP')); % https://github.com/jeheelab/TAFKAP
% Toolboxes for single trial BOLD amplitude estimation
addpath(genpath('~/Documents/GitHub/GLMsingle'));
addpath(genpath('~/Documents/GitHub/fracridge'));

addpath(genpath('~/Documents/GitHub/rokers_mri_lab/code/invChol'));

 BASE = '/Users/pw1246/Desktop/motion/';
%BASE = '~/Dropbox (RVL)/MRI/Decoding/';
addpath(genpath(BASE));

% Figure defaults
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize', 14)

% Set up parameters
sub = 'sub-0204';
ses = {'03'};
run = [1:10]';

roi =  {'V1','hMT'};
%roi =  {'V1','V2','hMT','IPS0'};
%roi =  {'V1','V2','V3','V3A','hV4','LO','hMT','MST','IPS'};
%roi = {'V1','V2','V3','V3A','V3B','hV4','LO1','LO2','hMT','MST','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','VO1','VO2','SPL1','PHC1','PHC2','FEF'};
[dataset,stim_label] = loadmydata(BASE,sub,ses,run,roi);


%% Estimate single trial BOLD amplitude 
% 
% Requires data in cell format: subject by run
% design = {};
% for p=1:10
%     design{p} = zeros(250,7);
%     if mod(p,2)==1
%         cnt = 1; curcond = 8;
%         while cnt <= 250
%             if curcond ~= 1
%              design{p}(cnt,curcond-1) = 1;
%             end
%             cnt = cnt + 2;
%             curcond = mod2(curcond-1,8);
%         end
%     else
%         cnt = 1; curcond = 1;
%         while cnt <= 250
%             if curcond ~= 1
%              design{p}(cnt,curcond-1) = 1;
%             end
%             cnt = cnt + 2;
%             curcond = mod2(curcond+1,8);
%         end
%     end
% end
% 
% design = cellfun(@(x) x(11:end,:),design,'UniformOutput',0);
% Func2 = cellfun(@(x) x(:,:,:,11:end),Func2,'UniformOutput',0);
% results = GLMestimatesingletrial(design,Func2,3,1.5,'testglmsingle',struct());

%% Run TAFKAP

nScans = length(ses)*length(run); % scans per subject
% % samples = cell(numel(roi),1);
% for whichRoi = 1:numel(roi)
%     samples{whichRoi} = cell2mat(DATA(:,whichRoi));
% end

samples = dataset;

% Setup design (parameters)
nFolds = 8; % Set to multiple of number of processing cores
params = SetupTAFKAP();
nDirs = 8; % motion directions
% params.stimval = 22.5.*reshape(repmat([5:-1:1 8:-1:6 4:8 1:3],1,nScans/2),1,[])'; % stimulus labels
%params.stimval = reshape(repmat([5:-1:1 8:-1:6 4:8 1:3],1,nScans/2),1,[])'; % categorical stimulus labels
params.stimval = stim_label;
params.runNs = reshape(repmat(1:nScans,nDirs,1),1,[])'; %trainruns; % stimulus block/run

pre = cell(nFolds,1);
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
for rr = 1:numel(roi)
    parfor ii = 1:nFolds
        rng(ii);% To counter the effects of TAFKAP_Decode setting the system rand seed to const. This was an EXTRAORDINARLY hard bug to find.
        [est{ii}, unc{ii}, liks{ii}, hypers{ii}] = TAFKAP_Decode(samples{rr}, p{ii});
    end
    ests{rr} = cell2mat(est');
    uncs{rr} = cell2mat(unc');
    pres{rr} = cell2mat(pre);
end

%% Calculate performance

params.subjects = sub;
% plot_results(user, params, rois, acc) % Plot results

for rr = 1:numel(roi)
    % Plot confusion matrix
    % conmat = confusionmat(categorical(ceil(pres{rr}/22.5)),categorical(ceil(ests{rr}/22.5))); % continuous data
    conmat = confusionmat(categorical(pres{rr}),categorical(ests{rr})); % categorical data
    % conmat = conmat.*nDirs./length(est);
    %conmat = 100.*conmat./nFolds; % convert to % accuracy
    conmat = 100.*conmat./sum(conmat,2); % convert to % accuracy
    
    figure(3); hold on;
    subplot(2,2,rr)
    conmat = [conmat; conmat(1,:)]; % wrap matrix
    conmat = [conmat, conmat(:,1)];
    % imagesc(conmat');
    clim = [0 100]; % upper, lower limits
    imagesc(conmat', clim);
    
    title(roi{rr})
    xlabel('Presented direction')
    ylabel('Decoded direction')
    
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    yticks([1:9])
    yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))    
    axis tight
    set(gca,'YDir','normal')
    
    if whichRoi == length(roi)
        cb = colorbar;
    end
    cb.Label.String = 'Classification performance (%)';
    
    disp(['Classification performance ' roi{rr} ': ' num2str(100.*mean(pres{rr}==ests{rr})) '%'])
end

matlab.graphics.internal.setPrintPreferences('DefaultPaperPositionMode','manual')
set(groot,'defaultFigurePaperPositionMode','manual')
saveas(gcf, ['../figures/Confusion_matrix-' datestr(now,30) '.pdf'])


% Plot uncertainty as a function of motion direction
figure
for rr = 1:numel(roi)
    subplot(2,2,rr)
    hold on
    title(roi{rr})
    scatter(pres{rr},uncs{rr})
    y = splitapply(@mean,uncs{rr},pres{rr});
    plot(1:8,y)
    
    xlim([.5 8.5])
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    
    xlabel('Presented motion direction')
    ylabel('Entropy (nats)')
end



%% MATLAB Classify

tic
code_rep = 1000; % bootstraps repeats
voxelmax = 1e5; % maximum voxels to use per ROI
test_n = 1; % percetange of data set as testing trials

trial_n = size(dataset{1,1},1);  % number of trials 
dir_n = 8; % number of directions
nScans = trial_n/dir_n; % number of scans
result = zeros(numel(roi),4); % allocate results
cnfm = cell(numel(roi),1); % each runs confusion matrix is restored separately

Training_group = sort(repmat(repelem(1:8,1)',nScans-test_n,1)); % training dataset stimulus label
Testing_group = repmat(repelem(1:8,1)',test_n,1); % testing dataset stimulus label 
random_Training_group = Training_group(randperm(length(Training_group))); % random shuffle stimulus label for decoding chance

pval = zeros(numel(roi),1); % p value (significance from chance)
pc = zeros(code_rep,numel(roi)); % percentage correct
chancepc = zeros(code_rep,numel(roi)); % percentage correct from chance


for ii = 1:numel(roi) % loop over ROIs

    cnfm{ii,1} = zeros(dir_n,dir_n);
    data = dataset{1,ii};
    
    [n,voxel_n] = size(data);
    
    for kk = 1:code_rep
        
        if voxel_n > voxelmax
            data = data(:,randperm(voxel_n,voxelmax),:); % randomly select (voxelmax) amount of voxels if there is a upper limit for voxel counts
        end

        Training_data = [];
        Testing_data = [];     
        
        for mm = 1:dir_n    
            train_run = 1:nScans;
            test_run = randperm(nScans,test_n); %random select (test_n) trial for each direction as testing dataset without replacement
            train_run(test_run) = [];           
            dir_idx = find(stim_label==mm);           
            Testing_data = [Testing_data; data(dir_idx(test_run),:)];
            Training_data = [Training_data; data(dir_idx(train_run),:)];            
        end
        
        % does this code only ever do 1 session at a time? (it will do 2 sessions if there are 2 sessions)
        
        decoded_group = classify(Testing_data,Training_data,Training_group,'diaglinear');
        random_decoded_group = classify(Testing_data,Training_data,random_Training_group,'diaglinear');
                
        pc(kk,ii) = length(find(decoded_group == Testing_group))/length(Testing_group);
        temp_chance = length(find(random_decoded_group == Testing_group))/length(Testing_group);
        chancepc(kk,ii) = temp_chance;
        temp_cnfm = zeros(dir_n,dir_n);
        for de = 1:dir_n
            for te = 1:dir_n
                temp_cnfm(de,te) = length(find(decoded_group == de & Testing_group == te))/length(find(Testing_group == te));
            end
        end
        cnfm{ii,1} = cnfm{ii,1} + temp_cnfm;
    end
    
    cnfm{ii,1} = cnfm{ii,1}./code_rep;
    
    [~,pval(ii),~,~] = ztest(pc(:,ii),mean(chancepc(:,ii)),std(chancepc(:,ii)),'Tail','Right');
       
    SEM = std(pc(:,ii))/sqrt(length(pc(:,ii)));    
    result(ii,1) = mean(pc(:,ii)).*100;
    result(ii,2) = SEM.*196*2;   
    SEM_chance = std(chancepc(:,ii))/sqrt(length(chancepc(:,ii)));
    result(ii,3) = mean(chancepc(:,ii)).*100;
    result(ii,4) = SEM_chance.*196*2;
    voxelsize(ii)=size(data,2);
end

toc
CIroi_mean = zeros(numel(roi),numel(roi));
CIroi = zeros(numel(roi),numel(roi));
pvalueroi = zeros(numel(roi),numel(roi));
for ii = 1:numel(roi);
    for jj = 1:numel(roi)        
        [~,pvalueroi(ii,jj)] = ttest(pc(:,ii),pc(:,jj));
        CIroi_mean(ii,jj) = mean(pc(:,ii)-pc(:,jj));
        CIroi(ii,jj) = (std(pc(:,ii)-pc(:,jj))/sqrt(code_rep))*1.96;      
    end
end

%% save data in the result folder
f = fullfile(pwd,'result',[sub '-ses-' ses{:} '-classify.mat']);
save(f,'cnfm','result','pvalueroi','pval','CIroi_mean','CIroi','roi','voxelsize');

%% plot confusion matrix
%close all
%whichroi = [1 4 9 11];
whichroi = 1:size(cnfm,1);
plot_confusionM(cnfm,roi,voxelsize,whichroi);
