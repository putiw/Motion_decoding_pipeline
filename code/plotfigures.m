clear all
close all
clc

% matlab classify (detrend vs. fft vs. detrend+norm)
BASE = '/Users/Puti/Documents/GitHub/classify_pipeline-/code/result/';
ROI = load(['/Users/Puti/Documents/GitHub/classify_pipeline-/code/result/sub-0204-ses-0102-TAFKAP.mat']);

BASE = '/Users/Puti/Documents/GitHub/plotf/';
BASE = '/Users/Puti/Documents/GitHub/classify_pipeline-/code/result/';

sub = {'01','03ny','04','05','06','29'};
sub = {'03ny'}
roi_n = 22;
ses_h = 'ses-0102';
ses_v = 'ses-0304';

acc = zeros(roi_n,numel(sub),2);
CNFM_c_t = cell(roi_n,numel(sub));
CNFM_i_t = cell(roi_n,numel(sub));
for subid = 1:length(sub)       
Coh = load([BASE 'sub-02' sub{subid} '-' ses_h '-TAFKAP.mat']);
Inc = load([BASE 'sub-02' sub{subid} '-' ses_v '-TAFKAP.mat']);

for kk = 1:roi_n  
    CNFM_c_t{kk,subid} = Coh.saveresult{kk,1}(1:8,1:8)./100;
    CNFM_i_t{kk,subid} = Inc.saveresult{kk,1}(1:8,1:8)./100;
    acc(kk,subid,1) = mean(diag(CNFM_c_t{kk,subid}))*100;
    acc(kk,subid,2) = mean(diag(CNFM_i_t{kk,subid}))*100;
end

end


roi_n = 22;
meanacc = mean(acc,2);
se = zeros(22,4);
for ii = 1:2
    for kk = 1:roi_n
se(kk,ii) = std(acc(kk,:,ii))/sqrt(length(sub));
    end
end

sum_fm = zeros(8,8);

cnfm_c_t = cell(roi_n,1);
cnfm_i_t = cell(roi_n,1);
for kk = 1:roi_n;
    sum_ct = sum_fm;
    sum_it = sum_fm;
    for subid = 1:length(sub)
        temp_c_t =  CNFM_c_t{kk,subid};
        temp_i_t =  CNFM_i_t{kk,subid}  ;
        sum_ct = sum_ct + temp_c_t;
        sum_it = sum_it + temp_i_t;
    end
    cnfm_c_t{kk,1} = sum_ct./length(sub);
    cnfm_i_t{kk,1} = sum_it./length(sub);
end


%% figure 1 example one subject V1 hMT accuracy bar graph

close all
clc

plotwhichroi = [1 9 11];
r_n = numel(plotwhichroi);
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],meanacc(plotwhichroi,:,1),0.6,'linewidth',2);
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',2);

bar1.FaceColor = 'flat';
bar1.CData(1,:) = [0.4 0.6 1];
bar1.CData(2,:) = [1 0.2 0.2];
bar1.CData(3,:) = [0.9290, 0.6940, 0.1250];
bar1.FaceAlpha = 0.9;
 plot([0 3.5],[12.5 12.5],'k--','LineWidth',2);
%  text(2.2,13.5,'chance')
 %text(1.5,52,'subject 0201','Fontsize',15,'HorizontalAlignment','center');
f1.CurrentAxes.XTick = [1.2 2 2.8];
f1.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
title({'Decoding Accuracy','(subject 0202)'})
xlabel('ROI')
box on 
ax = gca;
ax.LineWidth = 2;

% xtickangle(45)
ylim([0 70])
xlim([0.5 3.5])
 
% figure
% sqrtRois = ceil(sqrt(numel(Coh.roi)));
% for whichRoi = 1:numel(Coh.roi)
%     subplot(sqrtRois,sqrtRois,whichRoi)
%     hold on
%     title(Coh.roi{whichRoi})
%     scatter(Coh.pres{whichRoi},Coh.uncs{whichRoi})
%     y = splitapply(@mean,Coh.uncs{whichRoi},Coh.pres{whichRoi});
%     plot(1:8,y)
%     
%     xlim([.5 8.5])
%     xticks([1:9])
%     xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%     
%     xlabel('Presented motion direction')
%     ylabel('Entropy (nats)')
% end
%% figure 5 - sub 1 , H & V

close all
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],meanacc(plotwhichroi,:,1),0.3,'linewidth',2);
bar2 = bar([1.2 2 2.8]+0.25,meanacc(plotwhichroi,:,2),0.3,'linewidth',2);
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',2);
errorbar(bar2.XData,bar2.YData,se(plotwhichroi,2),'k.','linewidth',2);
bar1.FaceColor = 'flat';bar2.FaceColor = 'flat';
bar1.CData(1,:) = [0.4 0.6 1];bar2.CData(1,:) = [0.4 0.6 1];
bar1.CData(2,:) = [1 0.2 0.2];bar2.CData(2,:) = [1 0.2 0.2];
bar1.CData(3,:) = [0.9290, 0.6940, 0.1250];bar2.CData(3,:) = [0.9290, 0.6940, 0.1250];
bar1.FaceAlpha = 0.9;bar2.FaceAlpha = 0.3;
 plot([0 4],[12.5 12.5],'k--','LineWidth',2);
%  text(2.2,13.5,'chance')
 %text(1.5,52,'subject 0201','Fontsize',15,'HorizontalAlignment','center');
f1.CurrentAxes.XTick = [1.2 2 2.8]+0.13;
f1.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
title({'Decoding Accuracy','Horizontal vs. Vertical'})

title({'Decoding Accuracy','(subject 0202)'})

xlabel('ROI')
box on 
ax = gca;
ax.LineWidth = 2;
ylim([0 70])
xlim([0.7 3.5])
 


 %% figure 2 - example one subject V1 hMT confusion matrix

cnfm = cnfm_c_t; roi = ROI.roi; 
voxelsize = ROI.voxelsize;
whichroi = [1 9 11];
figure('Renderer', 'painters', 'Position', [10 10 1000 310]);
ha1 = tight_subplot(1,3,[0.05 0.05],[0.2 0.15],[0.07 .05]);

titleroi = roi(whichroi);
for mm = 1:length(whichroi);
    axes(ha1(mm));
    kk = whichroi(mm);
    
    combine_cfmx = cnfm{kk,1}.*100;
    
    combine_cfmx(9,:) = combine_cfmx(1,:);
    combine_cfmx(:,9) = combine_cfmx(:,1);
    
    imagesc(combine_cfmx)
    
    c1 = 100/8; %chance
    c2 = 50; %max
    cb = colorbar(); caxis([0 c2]);
    cb.Ruler.TickLabelFormat = '%d%%';
    
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
    xticklabels(cellstr([{char(8600)} {char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} ]))
    text(0.6:8.6,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
    set(gca,'YDir','normal')
    ylabel('Decoded Direction')
    xlabel('Presented Direction')
    axis square;
   
    if voxelsize == 0
        title(roi{kk});
    else
    Title = [roi{kk} ' (' num2str(voxelsize(kk)) ' voxels)'];
    %Title = [roi{kk}];
    title({'Horizontal',Title});
    end
      box on 
    set(gca,'FontSize',15,'LineWidth',2)
    hold off
end
%% figure 6 - example one subject V1 hMT confusion matrix (H & V)

cnfm = cnfm_i_t; roi = Coh.roi; 
voxelsize = ROI.voxelsize;
whichroi = [1 9 11];
figure('Renderer', 'painters', 'Position', [10 10 1000 310]);
ha1 = tight_subplot(1,3,[0.05 0.05],[0.2 0.15],[0.07 .05]);
titleroi = roi(whichroi);
for mm = 1:length(whichroi);
    axes(ha1(mm));
    kk = whichroi(mm);
    
    combine_cfmx = cnfm{kk,1}.*100;
    
    combine_cfmx(9,:) = combine_cfmx(1,:);
    combine_cfmx(:,9) = combine_cfmx(:,1);
    
    imagesc(combine_cfmx)
    
    c1 = 100/8; %chance
    c2 = 50; %max
    cb = colorbar(); caxis([0 c2]);
    cb.Ruler.TickLabelFormat = '%d%%';
    
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
    xticklabels(cellstr([{char(8600)} {char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} ]))
    text(0.6:8.6,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
    set(gca,'YDir','normal')
    ylabel('Decoded Direction')
    xlabel('Presented Direction')
    axis square;
    h=gca; h.XAxis.TickLength = [0 0];
    if voxelsize == 0
        title(roi{kk});
    else
    Title = [roi{kk} ' (' num2str(voxelsize(kk)) ' voxels)'];
    title({'Vertical',Title});
    end
      box on 
    set(gca,'FontSize',15,'LineWidth',2)
    hold off
end



%% figure 3 -  all subjects accuracy bar graph


clear all
close all
clc

% matlab classify (detrend vs. fft vs. detrend+norm)
BASE = '/Users/Puti/Documents/GitHub/classify_pipeline-/code/result/';
ROI = load([BASE 'sub-0201-ses-0102-TAFKAP.mat']);

BASE = '/Users/Puti/Documents/GitHub/plotf/';
BASE = '/Users/Puti/Documents/GitHub/classify_pipeline-/code/result/';

sub = {'01','05','06','29'}; %sub = {'01','03ny','04','05','06','29'};
% sub = {'02'}
roi_n = 22;
ses_h = 'ses-0102';
ses_v = 'ses-0304';

acc = zeros(roi_n,numel(sub),2);
CNFM_c_t = cell(roi_n,numel(sub));
CNFM_i_t = cell(roi_n,numel(sub));
for subid = 1:length(sub)       
Coh = load([BASE 'sub-02' sub{subid} '-' ses_h '-TAFKAP.mat']);
Inc = load([BASE 'sub-02' sub{subid} '-' ses_v '-TAFKAP.mat']);

for kk = 1:roi_n  
    CNFM_c_t{kk,subid} = Coh.saveresult{kk,1}(1:8,1:8)./100;
    CNFM_i_t{kk,subid} = Inc.saveresult{kk,1}(1:8,1:8)./100;
    acc(kk,subid,1) = mean(diag(CNFM_c_t{kk,subid}))*100;
    acc(kk,subid,2) = mean(diag(CNFM_i_t{kk,subid}))*100;
end

end


roi_n = 22;
meanacc = mean(acc,2);
se = zeros(22,1);
for ii = 1:2
    for kk = 1:roi_n
se(kk,ii) = std(acc(kk,:,ii))/sqrt(length(sub));
    end
end

sum_fm = zeros(8,8);
cnfm_c_c = cell(roi_n,1);
cnfm_i_c = cell(roi_n,1);
cnfm_c_t = cell(roi_n,1);
cnfm_i_t = cell(roi_n,1);
for kk = 1:roi_n;
    sum_ct = sum_fm;
    sum_it = sum_fm;
    for subid = 1:length(sub)
        temp_c_t =  CNFM_c_t{kk,subid};
        temp_i_t =  CNFM_i_t{kk,subid}  ;
        sum_ct = sum_ct + temp_c_t;
        sum_it = sum_it + temp_i_t;
    end
    cnfm_c_t{kk,1} = sum_ct./length(sub);
    cnfm_i_t{kk,1} = sum_it./length(sub);
end
%% figure 3
clc 
close all
f3 = figure('Renderer', 'painters', 'Position', [10 10 650 310]);
hold on
bar1 = bar(1:roi_n,meanacc(:,:,1),0.8,'linewidth',2);
errorbar(bar1.XData,bar1.YData,se(:,1),'k.','linewidth',2);

bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
 plot([0 23],[12.5 12.5],'k--','LineWidth',2);
f3.CurrentAxes.XTick = (1:roi_n);
f3.CurrentAxes.XTickLabel = Coh.roi;
box on 
ylim([0 59])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Decoding accuracy averaged across six subjects'})
xlabel('ROI')
set(gca,'FontSize',15,'LineWidth',2)

%% figure 4 - all confusion matrix
voxelsize = 0;
whichroi = 1:22;
plot_confusionM(cnfm_c_t,Coh.roi,voxelsize,whichroi)


%% figure 7
close all
f3 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,meanacc(:,:,1),0.4,'k','linewidth',2);
bar2 = bar((1:roi_n)+0.4,meanacc(:,:,2),0.4,'k','linewidth',2);
errorbar(bar1.XData,bar1.YData,se(:,1),'k.','linewidth',2);

errorbar(bar2.XData,bar2.YData,se(:,2),'k.','linewidth',2);
 plot([0 23],[12.5 12.5],'k--','LineWidth',2);
bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
bar2.FaceColor = [0 0.1 0];
bar2.FaceAlpha = 0.1;
legend('Horizontal','Vertical')
f3.CurrentAxes.XTick = (1:roi_n)+0.4;
f3.CurrentAxes.XTickLabel = Coh.roi;
box on 
ylim([0 59])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Decoding accuracy averaged across six subjects'})
xlabel('ROI')
set(gca,'FontSize',15,'LineWidth',2)


%%
close all
voxelsize = 0;
whichroi = 1:22;

% plot_confusionM(cnfm_c_c,Coh.roi,voxelsize,whichroi)
% text(18,8,'Horizontal (classify)','FontSize',25)
% plot_confusionM(cnfm_i_c,Coh.roi,voxelsize,whichroi)
% text(18,8,'Vertical (classify)','FontSize',25)
plot_confusionM(cnfm_c_t,Coh.roi,voxelsize,whichroi)
text(18,8,'Horizontal (TAFKAP)','FontSize',25)
 plot_confusionM(cnfm_i_t,Coh.roi,voxelsize,whichroi)
% text(18,8,'Vertical (TAFKAP)','FontSize',25)

%% motion index
close all
clc
indexacc = acc-12.5;
indexacc(indexacc<=0)=0;
a_b = indexacc(:,:,1)-indexacc(:,:,2);
ab = indexacc(:,:,1)+indexacc(:,:,2);
index(:,:) = a_b./ab;

meanidx = mean(index,2);
see = zeros(22,1);

    for kk = 1:roi_n
see(kk,1) = std(index(kk,:,1))/sqrt(length(sub));
    end

close all
f3 = figure('Renderer', 'painters', 'Position', [10 10 650 310]);
hold on

f3.CurrentAxes.XTick = (1:roi_n);
f3.CurrentAxes.XTickLabel = Coh.roi;
title('Motion index')
xtickangle(45)
% ylim([-1 1.5])
 set(gca,'FontSize',15,'LineWidth',2)
s = scatter(1:22,index,80,'o','filled','MarkerEdgeColor',[1 1 1],'LineWidth',2,'MarkerFaceAlpha',0.5);
errorbar(1:roi_n,meanidx(:,:),see(:,1),'ko','linewidth',2);
scatter(1:roi_n,meanidx(:,:),55,'ko','MarkerFaceColor',[1 1 1],'LineWidth',2);
xlim([0 23])
plot([0 23],[0 0],'k--')
ylabel('Motion index')
xlabel('ROI')
h=gca; h.XAxis.TickLength = [0 0];
box on 


%% motion index only v1 hMT IPS

close all
f11 = figure('Renderer', 'painters', 'Position', [10 10 310 310]);
hold on
plotwhichroi = [1 9 11];
f11.CurrentAxes.XTick = (1:numel(plotwhichroi));
f11.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
title('Motion Index')

% ylim([-1 1.5])
 set(gca,'FontSize',15,'LineWidth',2)
s = scatter(1:numel(plotwhichroi),index(plotwhichroi,:,:),100,'o','filled','MarkerEdgeColor',[1 1 1],'LineWidth',2,'MarkerFaceAlpha',0.5);
errorbar(1:numel(plotwhichroi),meanidx(plotwhichroi,:),see(plotwhichroi,1),'ko','linewidth',2);
scatter(1:numel(plotwhichroi),meanidx(plotwhichroi,:),65,'ko','MarkerFaceColor',[1 1 1],'LineWidth',2);
xlim([0 4])
plot([0 4],[0 0],'k--')
ylabel('Motion Index')
xlabel('ROI')
h=gca; h.XAxis.TickLength = [0 0];
box on 




%% ACC_H - ACC_V, all Ss
close all
f10 = figure('Renderer', 'painters', 'Position', [10 10 650 310]);
hold on

diff = acc(:,:,1)-acc(:,:,2);
bar1 = bar(1:roi_n,mean(diff,2),0.8,'linewidth',2);
errorbar(bar1.XData,bar1.YData,std(diff')./sqrt(6),'k.','linewidth',2);

bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
f10.CurrentAxes.XTick = 1:22;
f10.CurrentAxes.XTickLabel = Coh.roi;
box on 
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Differences in decoding accuracy: horizontal vs. vertical'})
xlabel('ROI')
set(gca,'FontSize',15,'LineWidth',2)



%% individual data with error bar 
clear all

BASE = '/Users/Puti/Documents/GitHub/classify_pipeline-/code/result/';

sub = {'01','04','05','06','29'};
sub = {'02'};
roi_n = 22;
ses_h = 'ses-0102';
ses_v = 'ses-0304';

acc = zeros(roi_n,length(sub),2);
se = zeros(roi_n,length(sub),2);
for subid = 1:length(sub)       
Coh = load([BASE 'sub-02' sub{subid} '-' ses_h '-TAFKAP.mat']);
Inc = load([BASE 'sub-02' sub{subid} '-' ses_v '-TAFKAP.mat']);

for kk = 1:roi_n  
    aa_h = squeeze(mean(reshape((Coh.pres{kk}==Coh.ests{kk}),16,8,[])));
    aa_v = squeeze(mean(reshape((Inc.pres{kk}==Inc.ests{kk}),16,8,[])));
    
    acc(kk,subid,1) = 100.*mean(aa_h);
    acc(kk,subid,2) = 100.*mean(aa_v);
    
    se(kk,subid,1) = 100.*std(aa_h)/sqrt(8);
    se(kk,subid,2) = 100.*std(aa_v)/sqrt(8);
end

end



roi_n = 22;
meanacc = mean(acc,2);
see = zeros(22,1);
for ii = 1:2
    for kk = 1:roi_n
see(kk,ii) = std(acc(kk,:,ii))/sqrt(length(sub));
    end
end


plotwhichroi = [1 9 11];
r_n = numel(plotwhichroi);
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],meanacc(plotwhichroi,:,1),0.6,'linewidth',2);
if numel(sub) ==1
    se = se;
else
    se = see;
end
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',2);

bar1.FaceColor = [0.4 0.6 1];
bar1.FaceAlpha = 0.9;
 plot([0 3.5],[12.5 12.5],'k--','LineWidth',2);
f1.CurrentAxes.XTick = [1.2 2 2.8];
f1.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
    Title = ['(subject 02' sub{1} ')'];
    title({'Decoding Accuracy',Title});
xlabel('ROI')
box on 
ax = gca;
ax.LineWidth = 2;

% xtickangle(45)
ylim([0 70])
xlim([0.5 3.5])
 
% figure
% sqrtRois = ceil(sqrt(numel(Coh.roi)));
% for whichRoi = 1:numel(Coh.roi)
%     subplot(sqrtRois,sqrtRois,whichRoi)
%     hold on
%     title(Coh.roi{whichRoi})
%     scatter(Coh.pres{whichRoi},Coh.uncs{whichRoi})
%     y = splitapply(@mean,Coh.uncs{whichRoi},Coh.pres{whichRoi});
%     plot(1:8,y)
%     
%     xlim([.5 8.5])
%     xticks([1:9])
%     xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%     
%     xlabel('Presented motion direction')
%     ylabel('Entropy (nats)')
% end
% figure 5 - sub 1 , H & V


f2 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],meanacc(plotwhichroi,:,1),0.3,'linewidth',2);
bar2 = bar([1.2 2 2.8]+0.248,meanacc(plotwhichroi,:,2),0.3,'linewidth',2);
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',2);
errorbar(bar2.XData,bar2.YData,se(plotwhichroi,2),'k.','linewidth',2);
bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.8;bar2.FaceAlpha = 0.75;
 plot([0 4],[12.5 12.5],'k--','LineWidth',2);
%  text(2.2,13.5,'chance')
 %text(1.5,52,'subject 0201','Fontsize',15,'HorizontalAlignment','center');
f2.CurrentAxes.XTick = [1.2 2 2.8]+0.13;
f2.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
title({'Decoding Accuracy','Horizontal vs. Vertical'})

    Title = ['(subject 02' sub{1} ')'];
    title({'Decoding Accuracy',Title});
    
l = legend('Horizontal','Vertical');
l.Position = l.Position + [0.02 0.025 0 0];
    
    
xlabel('ROI')
box on 
ax = gca;
ax.LineWidth = 2;
ylim([0 70])
xlim([0.7 3.5])
 
