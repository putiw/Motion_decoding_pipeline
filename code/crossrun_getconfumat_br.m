clear all
close all
clc
tic

% Dependencies
restoredefaultpath
addpath(genpath('~/Documents/GitHub/TAFKAP')); % https://github.com/jeheelab/TAFKAP

% BASE = '/Users/pw1246/Desktop/motion/';
BASE = '~/Dropbox (RVL)/MRI/Decoding/';

addpath(genpath(BASE));

sub = 'sub-0201';
ses = {'01'}; %,'02'};
run = [1:10]';%,[1:10]'];
roi =  {'V1','V2','hMT','IPS0'};
%roiname =  {'V1','V2','V3','V3A','hV4','LO','hMT','MST','IPS'};
%roiname = {'V1','V2','V3','V3A','V3B','hV4','LO1','LO2','hMT','MST','IPS0','IPS1','IPS2','IPS3','IPS4','IPS5','VO1','VO2','SPL1','PHC1','PHC2','FEF'};
params = SetupTAFKAP(); 
DATA = loadmydata(BASE,sub,ses,run,roi,params);


%% TAFKAP

nScans = length(ses)*length(run); % scans per subject
for whichRoi = 1:numel(roi)
    samples{whichRoi} = [];
    for whichSession = 1:nScans
        tSeries = DATA{whichSession,whichRoi}; % extract tSeries
        
        % downsample from TRs (240) x voxels to average response per run per direction (8) x voxels
        tSeries = (tSeries(1:2:end-1,:) + tSeries(2:2:end,:)) ./2; % take average of every 2 TRs
        % probably needs a mean([],2)
        tSeries = squeeze(mean(reshape(tSeries,8,15,[]),2)); % average every 8th datapoint
        samples{whichRoi} = [samples{whichRoi}; tSeries]; % build samples across runs
    end
end

% run TAFKAP decoder
params = SetupTAFKAP();
nDirs = 8; % motion directions
% params.stimval = 22.5.*reshape(repmat([5:-1:1 8:-1:6 4:8 1:3],1,nScans/2),1,[])'; % stimulus labels
params.stimval = reshape(repmat([5:-1:1 8:-1:6 4:8 1:3],1,nScans/2),1,[])'; % categorical stimulus labels
%params.stimval = 22.5.*reshape(repmat([1:8 8:-1:1],1,nScans/2),1,[])'; % stimulus labels (wrong but placeholder)
params.runNs = reshape(repmat(1:nScans,nDirs,1),1,[])'; %trainruns; % stimulus block/run

nFolds = 10;
for rr = 1:numel(roi)
    temp_est = []; % estimated direction
    temp_pre = []; % presented direction
    for ii = 1:nFolds
        rng(ii);% To counter the effects of TAFKAP_Decode setting the system rand seed to const. This was an EXTRAORDINARLY hard bug to find.
        c = cvpartition(params.stimval, 'Holdout', 0.1); % stratify by motion direction, but not scan
        params.train_trials = c.training;
        params.test_trials = c.test;
        [est, unc, liks, hypers] = TAFKAP_Decode(samples{rr}, params);
        temp_est = [temp_est; est];
        pre = params.stimval(c.test);
        temp_pre = [temp_pre; pre];
    end
    ests{rr} = temp_est;
    pres{rr} = temp_pre;
end

%% Calculate performance
% i = 1; % nFolds
% for mm = 1 % ROIs
%     acc(i,mm) = mean(get_acc(ests, params)); % accuracy: boot x roi
% end
user = 'br';
rois = 'V1';
params.subjects = sub;
% plot_results(user, params, rois, acc) % Plot results

for rr = 1:numel(roi)
    % Plot confusion matrix
    % conmat = confusionmat(categorical(ceil(pres{rr}/22.5)),categorical(ceil(ests{rr}/22.5))); % continuous data
    conmat = confusionmat(categorical(pres{rr}),categorical(ests{rr})); % categorical data
    % conmat = conmat.*nDirs./length(est);
    conmat = 100.*conmat./nFolds;
    
    figure(3); hold on;
    subplot(2,2,rr)
    conmat = [conmat; conmat(1,:)]; % wrap matrix
    conmat = [conmat, conmat(:,1)];
    imagesc(conmat');
    % clim = [0 .5]; % upper, lower limits
    % imagesc(conmat, clim);
    
    title(roi{rr})
    xlabel('Presented direction')
    ylabel('Decoded direction')
    
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    
    yticks([1:9])
    yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    
    axis tight
    % title(roi(whichRoi))
    
    if whichRoi == length(roi)
        cb = colorbar;
    end
    cb.Label.String = 'Classification performance (%)';
    
    disp(['Classification performance ' roi{rr} ':' num2str(mean(pres{rr}==ests{rr}))])
end

%%
code_rep = 1000; % bootstraps??
voxelmax = 10000;
trial_n = size(DATA{1,1},1)/2;
dir_n = 8;
dir_rep = trial_n/dir_n;
runn = numel(run);
gg = [[5:-1:1 8:-1:6]' [4:8 1:3]']; % true trial design matrix across all runs %lh



result = zeros(size(DATA,2),12);
CNFM = cell(runn,1); % each runs confusion matrix is restored separately
CNFM_random = cell(runn,1);

CNFM_std = cell(runn,1);
CNFM_std_random = cell(runn,1);

test_n = 1;
Training_group = repmat(repelem(1:8,1)',runn-test_n,1);
Testing_group = repmat(repelem(1:8,1)',test_n,1);

pval = zeros(size(DATA,2),1);
pc = zeros(code_rep,size(DATA,2));pc2 = zeros(code_rep,size(DATA,2));pc3 = zeros(code_rep,size(DATA,2));pc4 = zeros(code_rep,size(DATA,2));
chancepc = zeros(code_rep,size(DATA,2));

for ii = 1:size(DATA,2) % loop over ROIs?
    clear Resp_dir temp_Res
    cnfm{ii,1} = zeros(dir_n,dir_n);
    
    for RUN = 1:runn
        run_group = repmat(gg(:,mod(RUN,2)+1),dir_rep,1);
        data = DATA{RUN,ii};
        [n,voxel_n] = size(data);
        voxelsize(ii)=voxel_n;
        index = 1:n;
        elem = ones(1,n/2).*2;
        endv = n-sum(elem);
        if(~endv)
            endv = [];
        end
        index = mat2cell(index,1,[elem,endv])';
        
        clear B temp_B temp_dir
        temp_B = cell2mat(cellfun(@(x) sum(data(x,:),1),index,'un',0)); % average 2 TRs
        
        variance8 = zeros(dir_n,voxel_n);
        for dir = 1:dir_n
            variance8(dir,:)= var(temp_B(find(run_group==dir),:));
        end
        sigma = sqrt((sum(variance8).*(dir_rep-test_n))/(trial_n-dir_n));
        B = bsxfun(@rdivide,temp_B,sigma);
        
        for dir = 1:dir_n % # directions
            idd = find(run_group==dir); % index of trials of this direction
            temp_dir = B(idd,:);
            temp_Res(dir,:,RUN) = mean(temp_dir,1);
        end
    end
    
    
    for kk = 1:code_rep
        clear Training_data Testing_data
        
        if voxel_n <= voxelmax;
            Resp_dir = temp_Res;
        else
            Resp_dir = temp_Res(:,randperm(voxel_n,100),:);
        end
        
        testrunn = randperm(size(temp_Res,3),test_n);
        run_idx = 1:size(temp_Res,3);
        run_idx(testrunn) = [];
        Training_group = repmat(repelem(1:8,1)',size(temp_Res,3)-test_n,1);
        
        temp_train = Resp_dir(:,:,run_idx);
        temp_test = Resp_dir(:,:,testrunn);
        split_train = num2cell(temp_train,[1,2]);
        split_test = num2cell(temp_test,[1,2]);
        Training_data = vertcat(split_train{:});
        Testing_data = vertcat(split_test{:});
        random_Training_group = Training_group(randperm(numel(Training_group),numel(Training_group)),:);
        
        %random select 1 trial for each direction as testing dataset
        
        Training_data = zeros(runn*dir_n-dir_n,size(Resp_dir,2));
        Testing_data = zeros(dir_n,size(Resp_dir,2));
        for mm = 1:dir_n    
            train_run = 1:runn;
            test_run = randi(runn);
            train_run(test_run) = [];
            Testing_data(mm,:) = Resp_dir(mm,:,test_run);
            trainruns = num2cell(Resp_dir(mm,:,train_run),[1,2]);
            Training_data((1:8:72)+mm-1,:)= vertcat(trainruns{:});
        end
        
        % does this code only ever do 1 session at a time?
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
        pc2(kk,ii) = mean(diag(flip(temp_cnfm)));
        pc3(kk,ii) = 0.25*mean(diag(temp_cnfm,-1))+0.25*mean(diag(temp_cnfm,1))+0.5*mean(diag(temp_cnfm));
        pc4(kk,ii) = 0.25*mean(diag(flip(temp_cnfm)))+0.25*mean(diag(flip(temp_cnfm),1))+0.5*mean(diag(flip(temp_cnfm),2));
    end
    
    
    cnfm{ii,1} = cnfm{ii,1}./code_rep;
    
    [~,pval(ii),~,~] = ztest(pc(:,ii),mean(chancepc(:,ii)),std(chancepc(:,ii)),'Tail','Right');
    
    
    SEM = std(pc(:,ii))/sqrt(length(pc(:,ii)));%SEM2 = std(pc2)/sqrt(length(pc));
    SEM3 = std(pc3(:,ii))/sqrt(length(pc));SEM4 = std(pc4(:,ii))/sqrt(length(pc));
    %CI = mean(pc) + ts*SEM;
    result(ii,1) = mean(pc(:,ii)).*100;
    result(ii,2) = SEM.*196*2;%(max(CI)-min(CI)).*100;
    %     result(ii,3) = mean(pc2(:,ii)).*100;
    %     result(ii,4) = SEM2.*196*2;
    result(ii,5) = mean(pc3(:,ii),1).*100;
    result(ii,6) = mean(pc4(:,ii),1).*100;
    result(ii,7) = SEM3.*196*2;
    result(ii,8) = SEM4.*196*2;
    result(ii,9) = mean(pc3(:,ii)-pc4(:,ii)).*100;
    result(ii,10) =(std(pc3(:,ii)-pc4(:,ii))/sqrt(code_rep))*196;
    result(ii,11) = mean(pc3(:,ii)-chancepc(:,ii)).*100;
    result(ii,12) = (std(pc3(:,ii)-chancepc(:,ii))/sqrt(code_rep))*196;
    result(ii,13) = mean(pc4(:,ii)-chancepc(:,ii)).*100;
    result(ii,14) = (std(pc4(:,ii)-chancepc(:,ii))/sqrt(code_rep))*196;
    
    SEM_chance = std(chancepc(:,ii))/sqrt(length(chancepc(:,ii)));
    result(ii,3) = mean(chancepc(:,ii)).*100;
    result(ii,4) = SEM_chance.*196*2;
end

toc
CIroi_mean = zeros(size(DATA,2),size(DATA,2));
CIroi = zeros(size(DATA,2),size(DATA,2));
pvalueroi = zeros(size(DATA,2),size(DATA,2));
for ii = 1:size(DATA,2)
    for jj = 1:size(DATA,2)
        
        [~,pvalueroi(ii,jj)] = ttest(pc(:,ii),pc(:,jj));
        CIroi_mean(ii,jj) = mean(pc(:,ii)-pc(:,jj));
        CIroi(ii,jj) = (std(pc(:,ii)-pc(:,jj))/sqrt(code_rep))*1.96;
        
        
    end
end

%%
f = fullfile(pwd,'result',[sub '-ses-' ses{:} '-classify.mat']);
save(f,'cnfm','result','pvalueroi','pval','CIroi_mean','CIroi','roi','voxelsize');


%%

close all
titleroi = roi;
lett = {'A','B','C','D','E','F'};
figure('Renderer', 'painters', 'Position', [10 10 1150 700]);
%     ha1 = tight_subplot(1,3,[-0.13 0.05],[-0.05 -0.068],[0.05 .05]);
%ha1 = tight_subplot(3,3,[0.09 -0.55],[0.05 0.05],[-0.25 -0.3]);

%     figure('Renderer', 'painters', 'Position', [10 10 1150 350]);
ha1 = tight_subplot(1,4,[-0.13 0.05],[-0.05 -0.068],[0.05 .05]);

RoI = [1 2 9 11];
RoI = 1:4;
for k1 = 1:length(RoI)
    axes(ha1(k1));
    kk = RoI(k1);
    combine_cfmx = cnfm{kk,1}.*100;
    
    combine_cfmx(9,:) = combine_cfmx(1,:);
    combine_cfmx(:,9) = combine_cfmx(:,1);
    
    imagesc(combine_cfmx)
    
    c1 = 100/dir_n; %chance
    c2 = 50; %max
    cb = colorbar(); caxis([0 c2]);
    cb.Ruler.TickLabelFormat = '%d%%';
    %colorgroup = [100 180 255; 255 255 255; 255 122 122]./255;
    colorgroup = [200 230 255; 255 255 255; 255 122 122]./255;
    ratio = (c2-c1)/c1;
    cell_len = 10;
    value1 = linspace(0, 1, cell_len);
    mymap1 = value1'*colorgroup(2,:)+(1 - value1)'*colorgroup(1,:);
    value2 = linspace(0, 1, round(cell_len.*ratio));
    mymap2 = value2'*colorgroup(3,:)+(1 - value2)'*colorgroup(2,:);
    mymap = [mymap1;mymap2];
    
    colormap(mymap);
    
    
    hold on
    yticks([1:9]);
    yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    xticks([1:9]);
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    
    text(0.42:8.42,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
    %text(0,10,lett(k1),'FontSize',15,'Fontweight','bold');
    set(gca,'YDir','normal')
    ylabel('decoded direction')
    xlabel('presented direction')
    axis square;
    %Title = [titleroi{kk} ' (' num2str(size(DATA{1,kk},2)) ' voxels)'];
    
    if voxelmax <=100
        if size(DATA{1,kk},2) <= 100
            voxelcount(kk) = size(DATA{1,kk},2);
        else
            voxelcount(kk) = 100;
        end
    else
        voxelcount(kk) = size(DATA{1,kk},2);
    end
    titleroi{kk}
    Title = [titleroi{kk} ' (' num2str(voxelcount(kk)) ' voxels)'];
    %Title =[{'EVC'},{'hMT+'},{'DVC'},{'VCC'}];
    title(Title);
    set(gca,'FontSize',15)
    
    hold off
end
%%
%     figure('Renderer', 'painters', 'Position', [10 10 800 700]);
%     %ha = tight_subplot(2,3,[-0.13 0.05],[-0.05 -0.068],[0.05 .05]);
%     ha2 = tight_subplot(3,3,[0.09 -0.55],[0.05 0.05],[-0.25 -0.3]);
%     %ha = tight_subplot(2,roicount/2,[0.2/roicount 0.4/roicount]);
%
%     for k = 1:(size(DATA,2)/2);
%         axes(ha2(k));
%         kk = k+9;
%         combine_cfmx = cnfm{kk,1}.*100;
%
%         combine_cfmx(9,:) = combine_cfmx(1,:);
%         combine_cfmx(:,9) = combine_cfmx(:,1);
%
%         imagesc(combine_cfmx)
%
%         c1 = 100/dir_n; %chance
%         c2 = 30; %max
%         cb = colorbar(); caxis([0 c2]);
%         cb.Ruler.TickLabelFormat = '%d%%';
%         %colorgroup = [100 180 255; 255 255 255; 255 122 122]./255;
%         colorgroup = [200 230 255; 255 255 255; 255 122 122]./255;
%         ratio = (c2-c1)/c1;
%         cell_len = 10;
%         value1 = linspace(0, 1, cell_len);
%         mymap1 = value1'*colorgroup(2,:)+(1 - value1)'*colorgroup(1,:);
%         value2 = linspace(0, 1, round(cell_len.*ratio));
%         mymap2 = value2'*colorgroup(3,:)+(1 - value2)'*colorgroup(2,:);
%         mymap = [mymap1;mymap2];
%
%         colormap(mymap);
%
%
%         hold on
%         yticks([1:9]);
%         yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%         xticks([1:9]);
%         xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%         %text(8.6*ones(8,1),1:8,num2str(mean(combine_cfmx,2),3),'FontSize',15);
%         text(0.42:8.42,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
%         set(gca,'YDir','normal')
%         ylabel('decoded direction')
%         xlabel('presented direction')
%         axis square;
%         %Title = [titleroi{kk} ' (' num2str(size(DATA{1,kk},2)) ' voxels)'];
%
%         if voxelmax <=100
%         if size(DATA{1,kk},2) <= 100
%             voxelcount(kk) = size(DATA{1,kk},2);
%         else
%             voxelcount(kk) = 100;
%         end
%         else
%             voxelcount(kk) = size(DATA{1,kk},2);
%         end
%
%         Title = [titleroi{kk} ' (' num2str(voxelcount(kk)) ' voxels)'];
%
%         title(Title)
%         set(gca,'FontSize',15)
%
%         hold off
%     end

%
%     f3 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
%
%         hold on
%         ba = bar(1:length(roiname),result(:,1),0.4,'k','Linewidth',2);
%
%
%         errorbar(1:length(roiname),result(:,1),result(:,2),'k.','Linewidth',3)
%         ylim([10 30])
%         plot([0 length(roiname)+1],[12.5 12.5],'b--','linewidth',1.5);
%
%
%         ba.FaceAlpha = 0.05; %ba2.FaceAlpha = 0.3;
%         f3.CurrentAxes.XTick = 1:length(roiname);
%         xticklabels(roiname)
%         title([{'Eight-way decoding accuracy'}])
%         ylabel('Decoding accuracy')
%         myRange = [10:5:35];
%         yticks(myRange)
%         yticklabels(strcat(cellstr(string(num2cell(myRange))),'%'))
%         legend('Decoding 3D','Errors in depth but not 2D')
%         set(gca,'FontSize',15)
%         xtickangle(45)
%                 hold on
%         ba2 = bar(1:length(roiname),(result(:,3)+result(:,4)./2),0.4,'k','Linewidth',2);
%         ba3 = bar(1:length(roiname),(result(:,3)-result(:,4)./2),0.4,'w','Linewidth',2);
%         ba4 = bar(1:length(roiname),(result(:,3)-result(:,4)./2),0.4,'k','Linewidth',2);
%         ba4.FaceAlpha = 0.05;
%         ba2.FaceAlpha = 0.3;
%  %%
% close all
%     f3 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
%     ha1 = tight_subplot(1,3,[0 0.05],[0.15 0.05],[0.08 .08]);
%      %ha1 = tight_subplot(1,4,[-0.13 0.05],[-0.05 -0.068],[0.05 .05]);
% axes(ha1(1));
%         hold on
%
%         ba = bar(1:length(roiname),result(:,5),0.4,'k','Linewidth',2);
%         ba2 = bar((1:length(roiname))+0.4,result(:,6),0.4,'k','Linewidth',2);
%         errorbar(1:length(roiname),result(:,5),result(:,7),'k.','Linewidth',3)
%         errorbar((1:length(roiname))+0.4,result(:,6),result(:,8),'k.','Linewidth',3)
%         ylim([10 30])
%         plot([0 length(roiname)+1],[12.5 12.5],'b--','linewidth',1.5);
%
%         ba.FaceAlpha = 0.05; ba2.FaceAlpha = 0.3;
%         f3.CurrentAxes.XTick = 1:length(roiname);
%         roin = {'EVC','hMT+','DVC','VVC'}
%         xticklabels(roin)
%         title([{'Eight-way decoding accuracy'}])
%         ylabel('Decoding accuracy')
%         myRange = [10:5:35];
%         yticks(myRange)
%         yticklabels(strcat(cellstr(string(num2cell(myRange))),'%'))
%         legend('Decoding 3D','Errors in depth but not 2D')
%         set(gca,'FontSize',15)
%         %%
%         xtickangle(45)
%
%         text(-1,31.2,'E','FontSize',15,'Fontweight','bold');
%
%        axis square;
%
%
%         axes(ha1(2));
%
%         cfmx_2d = diag(ones(1,9))+diag(0.5*ones(1,8),-1)+diag(0.5*ones(1,8),1);
%         imagesc(cfmx_2d)
%         hold on
%         text(-1,0,'F','FontSize',15,'Fontweight','bold');
%         yticks([1:9]);
%         yticklabels(cellstr([{char(8594)} {char(8600)} {char(8595)} {char(8601)} {char(8592)} {char(8598)} {char(8593)} {char(8599)} {char(8594)}]))
%         xticks([1:9]);
%         xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%         title('Errors in depth but not 2D');
%                 ylabel('decoded direction')
%         xlabel('presented direction')
%         set(gca,'FontSize',15)
%         axis square;
%         axes(ha1(3));
%
%         cfmx_3d = flip(cfmx_2d);
%         imagesc(cfmx_3d)
%         hold on
%         text(-1,0,'G','FontSize',15,'Fontweight','bold');
%         cmap = gray(256);colormap(cmap);colorbar;
%         yticks([1:9]);
%         yticklabels(cellstr([{char(8594)} {char(8600)} {char(8595)} {char(8601)} {char(8592)} {char(8598)} {char(8593)} {char(8599)} {char(8594)}]))
%         xticks([1:9]);
%         xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%         title('Decoding 3D');        ylabel('decoded direction')
%         xlabel('presented direction')
%         set(gca,'FontSize',15)
% axis square;
% ha1(1).Position(1:4)
% ha1(2).Position(1:4)
% ha1(3).Position(3:4) = ha1(1).Position(3:4)