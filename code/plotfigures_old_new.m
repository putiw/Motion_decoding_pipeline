
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2],[mean(bar_H)*100,mean(bar_V)*100],0.6,'linewidth',1);

% errorbar(bar1.XData,bar1.YData,'k.','linewidth',1);

bar1.FaceColor = [0.4 0.6 1];
bar1.FaceAlpha = 0.9;
 plot([0 3.5],[12.96 12.96],'k--','linewidth',1);
f1.CurrentAxes.XTick = [1.2 2];
% f1.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
Title = ['(subject 02' sub{1} ')'];
title({'Decoding Accuracy',Title});
xlabel('ROI')
box on
ax = gca;
ax.LineWidth = 2;

% xtickangle(45)
ylim([0 50])
xlim([0.5 2.7])
%% xxxxxxx
%% XXXXXXX  individual data with error bar
clear all
close all

BASE = [pwd,'/result/'];
BASEnew = [pwd,'/braimcore_result/'];

sub = {'0201','0202','0204','0205','0206','0228','0229','0203ny'};
subnew = {'0201','0202','0204','0205','0206','0228','0229','0903'};
% sub = {'01'};
roi_n = 22;
ses = 'ses-0304';
ses_h = 'ses-0102';
ses_v = 'ses-0304';
sub_n = length(sub);
CNFM_c_t = cell(roi_n,sub_n);
CNFM_i_t = cell(roi_n,sub_n);

acc = zeros(roi_n,length(sub),2);
se = zeros(roi_n,length(sub),2);
for subid = 1:length(sub)
    Coh = load([BASEnew 'sub-' subnew{subid} '-' ses_h '-TAFKAP.mat']);
    Inc = load([BASEnew 'sub-' subnew{subid} '-' ses_v '-TAFKAP.mat']);
    
    for kk = 1:roi_n
        CNFM_c_t{kk,subid} = Coh.saveresult{kk,1}(1:8,1:8)./100;
        CNFM_i_t{kk,subid} = Inc.saveresult{kk,1}(1:8,1:8)./100;
        
        aa_h = squeeze(mean(reshape((Coh.pres{kk}==Coh.ests{kk}),16,32,[])));
        aa_v = squeeze(mean(reshape((Inc.pres{kk}==Inc.ests{kk}),16,32,[])));
        
        acc(kk,subid,1) = 100.*mean(aa_h);
        acc(kk,subid,2) = 100.*mean(aa_v);
        
        se(kk,subid,1) = 100.*std(aa_h)/sqrt(32);
        se(kk,subid,2) = 100.*std(aa_v)/sqrt(32);
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


sum_fm = zeros(8,8);

cnfm_c_t = cell(roi_n,1);
cnfm_i_t = cell(roi_n,1);

miss2 = zeros(roi_n,sub_n);
miss3 = zeros(roi_n,sub_n);

offdiag = zeros(roi_n,sub_n);

idx2x = [1 7 8 1 7 8 3 4 5 3 4 5 3 4 5 1 7 8];
idx2y = [1 1 1 2 2 2 4 4 4 5 5 5 6 6 6 8 8 8];
idx3x = [1 2 3 1 2 3 1 2 3 5 6 7 5 6 7 5 6 7];
idx3y = [2 2 2 3 3 3 4 4 4 6 6 6 7 7 7 8 8 8];

for kk = 1:roi_n;
    sum_ct = sum_fm;
    sum_it = sum_fm;
    for subid = 1:length(sub)
        temp_c_t =  CNFM_c_t{kk,subid};
        temp_i_t =  CNFM_i_t{kk,subid};
        sum_ct = sum_ct + temp_c_t;
        sum_it = sum_it + temp_i_t;

        % miss classification
        miss2(kk,subid) = (CNFM_c_t{kk,subid}(1,2) + CNFM_c_t{kk,subid}(5,4) + CNFM_c_t{kk,subid}(5,6) + CNFM_c_t{kk,subid}(1,8))/4;
        miss3(kk,subid) = (CNFM_c_t{kk,subid}(3,2) + CNFM_c_t{kk,subid}(3,4) + CNFM_c_t{kk,subid}(7,6) + CNFM_c_t{kk,subid}(7,8))/4;

        % off diag
        offdiag(kk,subid) = (CNFM_c_t{kk,subid}(1,1)+CNFM_c_t{kk,subid}(8,2)+CNFM_c_t{kk,subid}(6,4)+CNFM_c_t{kk,subid}(5,5)+CNFM_c_t{kk,subid}(4,6)+CNFM_c_t{kk,subid}(2,8))/6;
temp2 = zeros(18,1);
temp3 = zeros(18,1);
        % correct 2/3d
for ii = 1:18
        temp2(ii) = CNFM_c_t{kk,subid}(idx2x(ii),idx2y(ii));
        temp3(ii) = CNFM_c_t{kk,subid}(idx3x(ii),idx3y(ii));
end
cor2(kk,subid) = mean(temp2);
        cor3(kk,subid) = mean(temp3);
    end
    cnfm_c_t{kk,1} = sum_ct./length(sub);
    cnfm_i_t{kk,1} = sum_it./length(sub);
end

%% figure 1 V1 hMT IPS
plotwhichroi = [1 9 11];
r_n = numel(plotwhichroi);
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],meanacc(plotwhichroi,:,1),0.6,'linewidth',1);
if sub_n ==1
    se = se;
else
    se = see;
end
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',1);

bar1.FaceColor = [0.4 0.6 1];
bar1.FaceAlpha = 0.9;
plot([0 3.5],[12.5 12.5],'k--','linewidth',1);
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
ylim([0 50])
xlim([0.5 3.5])

%% figure 1 V1 hMT
close all
plotwhichroi = [1 9];
r_n = numel(plotwhichroi);
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2],meanacc(plotwhichroi,:,1),0.6,'linewidth',1);
if sub_n ==1
    se = se;
else
    se = see;
end
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',1);

bar1.FaceColor = [0.4 0.6 1];
bar1.FaceAlpha = 0.9;
plot([0 3.5],[12.5 12.5],'k--','linewidth',1);
f1.CurrentAxes.XTick = [1.2 2];
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
ylim([0 50])
xlim([0.5 2.7])
%% figure 6 V1 hMT IPS H/V
plotwhichroi = [1 9 11];
f2 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],meanacc(plotwhichroi,:,1),0.3,'linewidth',1);
bar2 = bar([1.2 2 2.8]+0.248,meanacc(plotwhichroi,:,2),0.3,'linewidth',1);
errorbar(bar1.XData,bar1.YData,se(plotwhichroi,1),'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,se(plotwhichroi,2),'k.','linewidth',1);
bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.8;bar2.FaceAlpha = 0.75;
plot([0 4],[12.5 12.5],'k--','linewidth',1);
%  text(2.2,13.5,'chance')
%text(1.5,52,'subject 0201','Fontsize',15,'HorizontalAlignment','center');
f2.CurrentAxes.XTick = [1.2 2 2.8]+0.13;
f2.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
title({'Decoding Accuracy','Horizontal vs. Vertical'})


if sub_n ==1
    Title = ['(subject 02' sub{1} ')'];
    title({'Decoding Accuracy',Title});
else
    title({'Decoding Accuracy','(Nine subjects)'});
end

l = legend('original','braimcore');
l.Position = l.Position + [0.02 0.025 0 0];


xlabel('ROI')
box on
ax = gca;
ax.LineWidth = 2;
ylim([0 50])
xlim([0.7 3.5])
%% FIGURE 4 bar all ROIs
clc
close all
f3 = figure('Renderer', 'painters', 'Position', [10 10 650 310]);
hold on
bar1 = bar(1:roi_n,meanacc(:,:,1),0.8,'linewidth',1);
errorbar(bar1.XData,bar1.YData,se(:,1),'k.','linewidth',1);

bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
plot([0 23],[12.5 12.5],'k--','linewidth',1);
f3.CurrentAxes.XTick = (1:roi_n);
f3.CurrentAxes.XTickLabel = Coh.roi;
box on
ylim([0 50])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Decoding accuracy averaged across nine subjects'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)

%%

%% figure 7 - example one subject V1 hMT confusion matrix
close all
cnfm = cnfm_c_t; roi = Coh.roi;
if sub_n == 1
    voxelsize = Coh.voxelsize;
else
    voxelsize = 0;
end

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
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    %text(0.6:8.6,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
    set(gca,'YDir','normal')
    ylabel('Decoded Direction')
    xlabel('Presented Direction')
    axis square;
    
    if voxelsize == 0
        Title = [roi{kk}];
        title({'Horizontal',Title});
        
    else
        Title = [roi{kk} ' (' num2str(voxelsize(kk)) ' voxels)'];
        title({'Horizontal',Title});
    end
    box on
    set(gca,'FontSize',15,'linewidth',1)
    hold off
end
% figure 8 - example one subject V1 hMT confusion matrix (H & V)


cnfm = cnfm_i_t; roi = Coh.roi;

whichroi = [1 3 9];
figure('Renderer', 'painters', 'Position', [10 10 1000 310]);
ha2 = tight_subplot(1,3,[0.05 0.05],[0.2 0.15],[0.07 .05]);
titleroi = roi(whichroi);
for mm = 1:length(whichroi);
    axes(ha2(mm));
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
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    %text(0.6:8.6,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
    set(gca,'YDir','normal')
    ylabel('Decoded Direction')
    xlabel('Presented Direction')
    axis square;
    h=gca; h.XAxis.TickLength = [0 0];
    if voxelsize == 0
        Title = [roi{kk}];
        title({'Vertical',Title});
        
    else
        Title = [roi{kk} ' (' num2str(voxelsize(kk)) ' voxels)'];
        title({'Vertical',Title});
    end
    box on
    set(gca,'FontSize',15,'linewidth',1)
    hold off
end

%% Misclassification 
plotwhichroi = [1 9 11];

mean_miss_2 = mean(miss2,2).*100;
sd_miss_2 = std(miss2')'/sqrt(sub_n-1).*100;
mean_miss_3 = mean(miss3,2).*100;
sd_miss_3 = std(miss3')'/sqrt(sub_n-1).*100;

f2 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],mean_miss_2(plotwhichroi),0.3,'linewidth',1);
bar2 = bar([1.2 2 2.8]+0.248,mean_miss_3(plotwhichroi),0.3,'linewidth',1);
errorbar(bar1.XData,bar1.YData,sd_miss_2(plotwhichroi),'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,sd_miss_3(plotwhichroi),'k.','linewidth',1);
bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.8;bar2.FaceAlpha = 0.75;
plot([0 4],[12.5 12.5],'k--','linewidth',1);
%  text(2.2,13.5,'chance')
%text(1.5,52,'subject 0201','Fontsize',15,'HorizontalAlignment','center');
f2.CurrentAxes.XTick = [1.2 2 2.8]+0.13;
f2.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
title({'Decoding Accuracy','Horizontal vs. Vertical'})


if sub_n ==1
    Title = ['(subject 02' sub{1} ')'];
    title({'Decoding Accuracy',Title});
else
    title({'Decoding Misclassification','(Nine subjects)'});
end

l = legend('2D','3D');
l.Position = l.Position + [0.02 0.025 0 0];


xlabel('ROI')
box on
ax = gca;
ax.LineWidth = 1;
ylim([0 30])
xlim([0.7 3.5])

%% off diag
close all
mean_off_diag = mean(offdiag,2).*100;
sd_off_diag = std(offdiag')'./(sub_n-1).*100;

plotwhichroi = [1 9 11];
r_n = numel(plotwhichroi);
f1 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],mean_off_diag(plotwhichroi),0.6,'linewidth',1);
if sub_n ==1
    se = se;
else
    se = see;
end
errorbar(bar1.XData,bar1.YData,sd_off_diag(plotwhichroi),'k.','linewidth',1);

bar1.FaceColor = [0.4 0.6 1];
bar1.FaceAlpha = 0.9;
plot([0 3.5],[12.5 12.5],'k--','linewidth',1);
f1.CurrentAxes.XTick = [1.2 2 2.8];
f1.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Decoding Percentage (%)')
set(gca,'FontSize',15)
Title = ['(subject 02' sub{1} ')'];
title({'Off-diagnal Decoding',Title});
xlabel('ROI')
box on
ax = gca;
ax.LineWidth = 1;

% xtickangle(45)
ylim([0 30])
xlim([0.5 3.5])


f3 = figure('Renderer', 'painters', 'Position', [10 10 650 310]);
hold on
bar1 = bar(1:roi_n,mean_off_diag,0.8,'linewidth',1);
errorbar(bar1.XData,bar1.YData,sd_off_diag,'k.','linewidth',1);

bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
plot([0 23],[12.5 12.5],'k--','linewidth',1);
f3.CurrentAxes.XTick = (1:roi_n);
f3.CurrentAxes.XTickLabel = Coh.roi;
box on
ylim([0 30])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Decoding Percentage (%)')
title({'Off-diagnal Decoding Percentage'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)


%% decode 2d


close all
mean_cor2 = mean(cor2,2).*100;
sd_cor2 = std(cor2')'./(sub_n-1).*100;
mean_cor3 = mean(cor3,2).*100;
sd_cor3 = std(cor3')'./(sub_n-1).*100;

plotwhichroi = [1 9 11];
f2 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on
bar1 = bar([1.2 2 2.8],mean_cor2(plotwhichroi),0.3,'linewidth',1);
bar2 = bar([1.2 2 2.8]+0.248,mean_cor3(plotwhichroi),0.3,'linewidth',1);
errorbar(bar1.XData,bar1.YData,sd_cor2(plotwhichroi),'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,sd_cor3(plotwhichroi),'k.','linewidth',1);
bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.8;bar2.FaceAlpha = 0.75;
plot([0 4],[12.5 12.5],'k--','linewidth',1);
%  text(2.2,13.5,'chance')
%text(1.5,52,'subject 0201','Fontsize',15,'HorizontalAlignment','center');
f2.CurrentAxes.XTick = [1.2 2 2.8]+0.13;
f2.CurrentAxes.XTickLabel = Coh.roi(plotwhichroi);
ylabel('Percentage Correct (%)')
set(gca,'FontSize',15)
title({'Decoding Accuracy','2D vs. 3D'})


if sub_n ==1
    Title = ['(subject 02' sub{1} ')'];
    title({'Decoding Accuracy',Title});
else
    title({'Decoding 2D vs. 3D','(Nine subjects)'});
end

l = legend('2D','3D');
l.Position = l.Position + [0.02 0.025 0 0];


xlabel('ROI')
box on
ax = gca;
ax.LineWidth = 1;
ylim([0 30])
xlim([0.7 3.5])







%% figure 3 - all confusion matrix
close all
voxelsize = 0;
whichroi = 1:22;
plot_confusionM(cnfm_c_t,Coh.roi,voxelsize,whichroi)
plot_confusionM(cnfm_i_t,Coh.roi,voxelsize,whichroi)

%% figure 7 bar all ROIs H/V
close all
f3 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,meanacc(:,:,1),0.4,'k','linewidth',1);
bar2 = bar((1:roi_n)+0.4,meanacc(:,:,2),0.4,'k','linewidth',1);




plot([bar1.XData;bar1.XData], [bar1.YData-se(:,1)',;bar1.YData+se(:,1)'], '-k','linewidth',1.5);
plot([bar2.XData;bar2.XData], [bar2.YData-se(:,2)',;bar2.YData+se(:,2)'], '-k','linewidth',1.5);

 plot([0 23],[12.5 12.5],'k--','linewidth',1);
bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
bar2.FaceColor = [0 0.1 0];
bar2.FaceAlpha = 0.1;
legend('Horizontal','Vertical')
f3.CurrentAxes.XTick = (1:roi_n)+0.25;
f3.CurrentAxes.XTickLabel = Coh.roi;
box on 
ylim([0 50])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Decoding accuracy averaged across nine subjects'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)


%%

close all
clc
indexacc = acc-12.5;
indexacc(indexacc<=0)=0;
a_b = indexacc(:,:,1)-indexacc(:,:,2);
ab = indexacc(:,:,1)+indexacc(:,:,2);
index(:,:) = a_b./ab;

meanidx = nanmean(index,2);
see = zeros(22,1);

    for kk = 1:roi_n
see(kk,1) = nanstd(index(kk,:,1))/sqrt(length(sub));
    end


f3 = figure('Renderer', 'painters', 'Position', [10 10 650 310]);
hold on

f3.CurrentAxes.XTick = (1:roi_n);
f3.CurrentAxes.XTickLabel = Coh.roi;
title('Motion index')
xtickangle(45)
% ylim([-1 1.5])
 set(gca,'FontSize',15,'linewidth',1)
 nanindex = ~isnan(index);
%  for ii = 1:22
%   roinanindex = index(ii,nanindex(ii,:));
% scatter(repmat(ii,numel(roinanindex),1),roinanindex,80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.5);
%  end
 roiidx = 1:22;
  for ii = 1:sub_n
  roinanindex = index(nanindex(:,ii),ii);
scatter(roiidx(nanindex(:,ii)),roinanindex,80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.5);
  end
 
  
errorbar(1:roi_n,meanidx(:,:),see(:,1),'ko','linewidth',1);
scatter(1:roi_n,meanidx(:,:),55,'ko','MarkerFaceColor',[1 1 1],'linewidth',1);
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
 set(gca,'FontSize',15,'linewidth',1)
%  nanindex = ~isnan(index);
%  for ii = 1:numel(plotwhichroi)
%      kk = plotwhichroi(ii);
%   roinanindex = index(kk,nanindex(kk,:));
% scatter(repmat(ii,numel(roinanindex),1),roinanindex,80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.5);
%  end
 
whichindex = index(plotwhichroi,:);
 nanindex = ~isnan(whichindex);
  for ii = 1:sub_n
  roinanindex = whichindex(nanindex(:,ii),ii);
scatter(1:3,roinanindex,80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.5);
  end
 
  
 errorbar(1:numel(plotwhichroi),meanidx(plotwhichroi,:),see(plotwhichroi,1),'ko','linewidth',1);
scatter(1:numel(plotwhichroi),meanidx(plotwhichroi,:),65,'ko','MarkerFaceColor',[1 1 1],'linewidth',1);
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
bar1 = bar(1:roi_n,mean(diff,2),0.8,'linewidth',1);
errorbar(bar1.XData,bar1.YData,std(diff')./sqrt(6),'k.','linewidth',1);

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
set(gca,'FontSize',15,'linewidth',1)

%%  entropy
% close all
% clear all
% clc

% BASE = '/Users/pw1246/Box Sync/plotf (pw1246@nyu.edu)/';
% sub = {'01','02','03ny','04','05','06','28','29','48'};
% sub = {'01'};
roi_n = 22;
ses_h = 'ses-0102';
ses_v = 'ses-0304';
sub_n = 9;
CNFM_c_t = cell(roi_n,sub_n);
CNFM_i_t = cell(roi_n,sub_n);

acc = zeros(roi_n,length(sub),2);
se = zeros(roi_n,length(sub),2);
y_h = zeros(roi_n,8,length(sub));
y_v = zeros(roi_n,8,length(sub));
std_y_h = zeros(roi_n,8,length(sub));
std_y_v = zeros(roi_n,8,length(sub));
for subid = 1:length(sub)
    Coh = load([BASEnew 'sub-' subnew{subid} '-' ses_h '-TAFKAP.mat']);
    Inc = load([BASEnew 'sub-' subnew{subid} '-' ses_v '-TAFKAP.mat']);
      
    for whichRoi = 1:numel(Coh.roi)
    temp_y_h  = splitapply(@mean,Coh.uncs{whichRoi},Coh.pres{whichRoi})';
    temp_y_v  = splitapply(@mean,Inc.uncs{whichRoi},Inc.pres{whichRoi})';
    std_temp_y_h  = splitapply(@std,Coh.uncs{whichRoi},Coh.pres{whichRoi})'./sqrt(512);
    std_temp_y_v  = splitapply(@std,Inc.uncs{whichRoi},Inc.pres{whichRoi})'./sqrt(512);
 
    y_h(whichRoi,1:8,subid) = temp_y_h; y_h(whichRoi,9,subid) = y_h(whichRoi,1,subid);
    y_v(whichRoi,1:8,subid) = temp_y_v; y_v(whichRoi,9,subid) = y_v(whichRoi,1,subid);
    std_y_h(whichRoi,1:8,subid) = std_temp_y_h; std_y_h(whichRoi,9,subid) = std_y_h(whichRoi,1,subid);
    std_y_v(whichRoi,1:8,subid) = std_temp_y_v; std_y_v(whichRoi,9,subid) = std_y_v(whichRoi,1,subid);

    end
end

mean_y_h = mean(y_h,3);
figure
hold on 
for whichRoi = 1:numel(Coh.roi)
    subplot(4,6,whichRoi)
    hold on
    title(Coh.roi{whichRoi})
scatter(1:9,mean_y_h(whichRoi,:),80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.5);

    xlim([.5 9.5])
    ylim([1.2 2])
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
    
    xlabel('Presented motion direction')
    ylabel('Entropy (nats)')
end


%%

gap = 0.25;

mean_sub_h = mean(y_h,3);
mean_sub_v = mean(y_v,3);

se_en = zeros(22,9,2);
if sub_n ==1
    se_en(:,:,1)=std_y_h;
    se_en(:,:,2)=std_y_v;
else
for kk = 1:22;
    se_en(kk,:,1)= std(reshape(y_h(kk,:,:),9,sub_n)')/sqrt(sub_n);
    se_en(kk,:,2)= std(reshape(y_v(kk,:,:),9,sub_n)')/sqrt(sub_n);
end 
end

close all
figure('Renderer', 'painters', 'Position', [10 10 1800 900]);

for whichRoi = 1:numel(Coh.roi)
    subplot(4,6,whichRoi)
    hold on
for subid = 1:length(sub)       
    title(Coh.roi{whichRoi})
    scatter((1:9)-gap,y_h(whichRoi,:,subid),70,'o','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.5 0.5 0.5],'linewidth',1.5,'MarkerFaceAlpha',0.3);

    xlim([.5 9.5])
     ylim([1.5 2.1])
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))    
    xlabel('Presented motion direction')
    ylabel('Entropy (nats)')
end

 plot((1:9)-gap,mean_sub_h(whichRoi,:),'k-','linewidth',1);
  errorbar((1:9)-gap,mean_sub_h(whichRoi,:),se_en(whichRoi,:,1),'ro','MarkerFaceColor',[1 1 1],'linewidth',1);
%scatter((1:9),mean_sub_h(whichRoi,:),60,'ro','MarkerFaceColor',[1 1 1],'linewidth',1);
 
end

%%

figure('Renderer', 'painters', 'Position', [10 10 1800 900]);

for whichRoi = 1:numel(Coh.roi)
    subplot(4,6,whichRoi)
    hold on
for subid = 1:length(sub)       
    title(Coh.roi{whichRoi})
    scatter((1:9)+gap,y_v(whichRoi,:,subid),70,'o','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.5 0.5 0.5],'linewidth',1.5,'MarkerFaceAlpha',0.3);

    xlim([.5 9.5])
    ylim([1.5 2.1])
    xticks([1:9])
    xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))    
    %xlabel('Presented motion direction')
    %ylabel('Entropy (nats)')
end
 plot((1:9)+gap,mean_sub_v(whichRoi,:),'k-','linewidth',1);

  errorbar((1:9)+gap,mean_sub_v(whichRoi,:),se_en(whichRoi,:,2),'ko','MarkerFaceColor',[1 1 1],'linewidth',1);
%  scatter((1:9)+gap,mean_sub_v(whichRoi,:),60,'ko','MarkerFaceColor',[1 1 1],'linewidth',1);
set(gca,'xticklabel',{[]})
set(gcf, 'color', 'none');   
set(gca, 'color', 'none');
end

%% lateral versus depth
y_h_4 = y_h(:,[2 4 6 8],:);
y_v_4 = y_v(:,[2 4 6 8],:);
mean_sub_h = mean(y_h_4,3); 
mean_sub_v = mean(y_v_4,3);
mean_all = zeros(22,4);
mean_all(:,1)  = mean(mean_sub_h(:,[1 3]),2);
mean_all(:,2)  = mean(mean_sub_h(:,[2 4]),2);
mean_all(:,3)  = mean(mean_sub_v(:,[1 3]),2);
mean_all(:,4)  = mean(mean_sub_v(:,[2 4]),2);

se_en = zeros(22,2,4);
for kk = 1:22;
    se_en(kk,:,1)= std(reshape(y_h_4(kk,[1 3],:),2,9)')/sqrt(sub_n);
    se_en(kk,:,2)= std(reshape(y_v_4(kk,[1 3],:),2,9)')/sqrt(sub_n);
    se_en(kk,:,3)= std(reshape(y_h_4(kk,[2 4],:),2,9)')/sqrt(sub_n);
    se_en(kk,:,4)= std(reshape(y_v_4(kk,[2 4],:),2,9)')/sqrt(sub_n);
end 
se_en = reshape(mean(se_en,2),22,4);

close all
figure('Renderer', 'painters', 'Position', [10 10 1500 900]);

for whichRoi = 1:numel(Coh.roi)
    subplot(4,6,whichRoi)
    hold on
for subid = 1:length(sub)       
    title(Coh.roi{whichRoi})
    scatter(1:4,y_h_4(whichRoi,:,subid),80,'o','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceAlpha',0.3);

    xlim([.5 2.5])
     ylim([1.2 2.1])
    xticks([1:2])
    xticklabels({'L','D'}); 
    xlabel('Presented motion direction')
    ylabel('Entropy (nats)')
end
 errorbar(1:2,mean_all(whichRoi,[1 2]),se_en(whichRoi,[1 3]),'ko','linewidth',1);
 scatter(1:2,mean_all(whichRoi,[1 2]),65,'ko','MarkerFaceColor',[1 1 1],'linewidth',1);

end

figure('Renderer', 'painters', 'Position', [10 10 1500 900]);

for whichRoi = 1:numel(Coh.roi)
    subplot(4,6,whichRoi)
    hold on
for subid = 1:length(sub)       
    title(Coh.roi{whichRoi})
    scatter(1:4,y_v_4(whichRoi,:,subid),80,'o','filled','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0.5 0.5 0.5],'linewidth',1,'MarkerFaceAlpha',0.3);


    xlim([.5 2.5])
     ylim([1.2 2.1])
    xticks([1:2])
    xticklabels({'L','D'}); 
    xlabel('Presented motion direction')
    ylabel('Entropy (nats)')
end
 errorbar(1:2,mean_all(whichRoi,[3 4]),se_en(whichRoi,[2 4]),'ko','linewidth',1);
 scatter(1:2,mean_all(whichRoi,[3 4]),65,'ko','MarkerFaceColor',[1 1 1],'linewidth',1);

end

%% figure 10 no use
close all


f3 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,mean_y_h,0.4,'k','linewidth',1);
bar2 = bar((1:roi_n)+0.4,mean_y_v,0.4,'k','linewidth',1);
errorbar(bar1.XData,bar1.YData,se_en(:,1),'k.','linewidth',1);

errorbar(bar2.XData,bar2.YData,se_en(:,2),'k.','linewidth',1);
bar1.FaceColor = [0 0.1 0];
bar1.FaceAlpha = 0.4;
bar2.FaceColor = [0 0.1 0];
bar2.FaceAlpha = 0.1;
legend('Horizontal','Vertical')
f3.CurrentAxes.XTick = (1:roi_n)+0.25;
f3.CurrentAxes.XTickLabel = Coh.roi;
box on 
ylim([1.2 2])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Entropy (nats)')
title({'Decoding uncertainty averaged across nine subjects'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)


%% pairwise 4 figures bar 
close all
clear all
clc


BASE = '/Users/Puti/Box Sync/plotf (pw1246@nyu.edu)/';
sub = {'01','02','03ny','04','05','06','28','29','48'};
%  sub = {'05'};
roi_n = 22;
ses_h = 'ses-0102';
ses_v = 'ses-0304';

acc = zeros(roi_n,length(sub),4);
acc_chance = zeros(roi_n,length(sub),2);
se = zeros(roi_n,length(sub),2);
y_h = zeros(roi_n,8,length(sub),2); % roi by subject by condition(2D/3D)
y_v = zeros(roi_n,8,length(sub),2);

for subid = 1:length(sub)
    Hor = load([BASE 'sub-02' sub{subid} '-' ses_h '-classify-pair4.mat']);
    Ver = load([BASE 'sub-02' sub{subid} '-' ses_v '-classify-pair4.mat']);
   
    acc(:,subid,1) = Hor.result(:,1); % 2D - hor
    acc(:,subid,2) = Hor.result(:,2); % 3D - hor
    acc(:,subid,3) = Ver.result(:,1); % 2D - ver
    acc(:,subid,4) = Ver.result(:,2); % 3D - ver
    
    acc_chance(:,subid,1) = Hor.result(:,5); % 2D - hor
    acc_chance(:,subid,2) = Hor.result(:,6); % 3D - hor
    acc_chance(:,subid,3) = Ver.result(:,5); % 2D - ver
    acc_chance(:,subid,4) = Ver.result(:,6); % 3D - ver
       
end

mean_pair_acc = zeros(roi_n,4,2);
mean_pair_se = zeros(roi_n,4,2);

for ii = 1:4
    mean_pair_acc(:,ii,1) = mean(acc(:,:,ii),2); % mean accuracy
    mean_pair_acc(:,ii,2) = mean(acc_chance(:,:,ii),2); % mean accuracy chance
    mean_pair_se(:,ii,1) = std(acc(:,:,ii)')'/sqrt(sub_n); % se
    mean_pair_se(:,ii,2) = std(acc_chance(:,:,ii)')'/sqrt(sub_n); %  se chance
end

if sub_n ==1
    clear mean_pair_se
    mean_pair_se(:,1,1) = Hor.result(:,3); %2d h
    mean_pair_se(:,2,1) = Hor.result(:,4); %3d
    mean_pair_se(:,1,2) = Hor.result(:,7); %chance
    mean_pair_se(:,2,2) = Hor.result(:,8); %chance
    mean_pair_se(:,3,1) = Ver.result(:,3); %2d v
    mean_pair_se(:,4,1) = Ver.result(:,4);
    mean_pair_se(:,3,2) = Ver.result(:,7); %chance
    mean_pair_se(:,4,2) = Ver.result(:,8); %chance
else
end

%% 4 figures bar 
for oo =1 
    
close all
f3 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,mean_pair_acc(:,1,1),0.4,'k','linewidth',1); %2D hor
bar2 = bar((1:roi_n)+0.4,mean_pair_acc(:,2,1),0.4,'k','linewidth',1); %3D hor
errorbar(bar1.XData,bar1.YData,mean_pair_se(:,1,1) ,'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,mean_pair_se(:,2,1),'k.','linewidth',1);
plot([0 23],[50 50],'k--','linewidth',1);
% bar1.FaceColor = [0 0.1 0];
% bar1.FaceAlpha = 0.4;
% bar2.FaceColor = [0 0.1 0];
% bar2.FaceAlpha = 0.1;
bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [0.4 0.6 1];
bar1.FaceAlpha = 0.8;bar2.FaceAlpha = 0.8/3;

legend('2D','3D')
f3.CurrentAxes.XTick = (1:roi_n)+0.25;
f3.CurrentAxes.XTickLabel = Hor.roi;
box on 
ylim([45 90])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Horizontal pair-wise decoding accuracy'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)

f11 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,mean_pair_acc(:,3,1),0.4,'k','linewidth',1); %2D ver
bar2 = bar((1:roi_n)+0.4,mean_pair_acc(:,4,1),0.4,'k','linewidth',1); %3D ver
errorbar(bar1.XData,bar1.YData,mean_pair_se(:,3,1) ,'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,mean_pair_se(:,4,1),'k.','linewidth',1);
plot([0 23],[50 50],'k--','linewidth',1);

bar1.FaceColor = [1 0.2 0.2];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.75;bar2.FaceAlpha = 0.75/3;

% bar1.FaceColor = [0 0.1 0];
% bar1.FaceAlpha = 0.4;
% bar2.FaceColor = [0 0.1 0];
% bar2.FaceAlpha = 0.1;
legend('2D','3D')
f11.CurrentAxes.XTick = (1:roi_n)+0.25;
f11.CurrentAxes.XTickLabel = Hor.roi;
box on 
ylim([45 90])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Vertical pair-wise decoding accuracy'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)

f12 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,mean_pair_acc(:,1,1),0.4,'k','linewidth',1); %  2D hor
bar2 = bar((1:roi_n)+0.4,mean_pair_acc(:,3,1),0.4,'k','linewidth',1); % 2D  ver
errorbar(bar1.XData,bar1.YData,mean_pair_se(:,1,1) ,'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,mean_pair_se(:,3,1),'k.','linewidth',1);
plot([0 23],[50 50],'k--','linewidth',1);
% bar1.FaceColor = [0 0.1 0];
% bar1.FaceAlpha = 0.4;
% bar2.FaceColor = [0 0.1 0];
% bar2.FaceAlpha = 0.1;
bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.8;bar2.FaceAlpha = 0.75;
legend('Horizontal','Vertical')
f12.CurrentAxes.XTick = (1:roi_n)+0.25;
f12.CurrentAxes.XTickLabel = Hor.roi;
box on 
ylim([45 90])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'2D pair-wise decoding accuracy'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)

f13 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:roi_n,mean_pair_acc(:,2,1),0.4,'k','linewidth',1); %  3D hor
bar2 = bar((1:roi_n)+0.4,mean_pair_acc(:,4,1),0.4,'k','linewidth',1); % 3D  ver
errorbar(bar1.XData,bar1.YData,mean_pair_se(:,2,1) ,'k.','linewidth',1);
errorbar(bar2.XData,bar2.YData,mean_pair_se(:,4,1),'k.','linewidth',1);
plot([0 23],[50 50],'k--','linewidth',1);

bar1.FaceColor = [0.4 0.6 1];bar2.FaceColor = [1 0.2 0.2];
bar1.FaceAlpha = 0.8/3;bar2.FaceAlpha = 0.75/3;
% 
% bar1.FaceColor = [0 0.1 0];
% bar1.FaceAlpha = 0.4;
% bar2.FaceColor = [0 0.1 0];
% bar2.FaceAlpha = 0.1;
legend('Horizontal','Vertical')
f13.CurrentAxes.XTick = (1:roi_n)+0.25;
f13.CurrentAxes.XTickLabel = Hor.roi;
box on 
ylim([45 90])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'3D pair-wise decoding accuracy'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)
end
%% scatter 
clc
close all
s = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on

for oo = 1;
    x = (1:22)';
y = (mean_pair_acc(:,1,2)+ mean_pair_acc(:,2,2))/2;
dy = (mean_pair_se(:,1,2)+ mean_pair_se(:,2,2))/2;  % made-up error values
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 0.9 0.9],'linestyle','none');


for subid = 1:sub_n
    scatter(1:22,acc(:,subid,1),80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.3);    
end

plot([x, x]', [mean_pair_acc(:,1,1)-mean_pair_se(:,1,1), mean_pair_acc(:,1,1)+mean_pair_se(:,1,1)]', '-k','LineWidth',4)
scatter(1:22,mean_pair_acc(:,1,1),75,'ko','MarkerFaceColor',[1 1 1],'LineWidth',4);
s.CurrentAxes.XTick = (1:roi_n)+0.25;
s.CurrentAxes.XTickLabel = Hor.roi;
box on
ylim([42 100])
xlim([0.5 22.5])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Pair-wise decoding accuracy'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)



%  add toward versus away
s2 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on

for oo = 1;
    x = (1:22)'+0.32;

    for subid = 1:sub_n
        scatter(x,acc(:,subid,2),80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.3);       
    end
    cp = plot([x, x]', [mean_pair_acc(:,2,1)-mean_pair_se(:,2,1), mean_pair_acc(:,2,1)+mean_pair_se(:,2,1)]','Color',[0.5 0.5 0.5],'LineWidth',4)

    plot([(1:21)'+0.67,(1:21)'+0.67]', [zeros(21,1), 100*ones(21,1)]', '-k','LineWidth',0.1) 
    scatter(x,mean_pair_acc(:,2,1),75,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',4);
    s2.CurrentAxes.XTick = (x)+0.25-0.3;
    s2.CurrentAxes.XTickLabel = Hor.roi;
    box on
    ylim([42 100])
    xlim([0.5 22.5])
    xtickangle(45)
    h=gca; h.XAxis.TickLength = [0 0];
    ylabel('Percentage Correct (%)')
    title({'Pair-wise decoding accuracy'})
    xlabel('ROI')
    set(gca,'FontSize',15,'linewidth',1)
    
end
set(gcf, 'color', 'none');   
set(gca, 'color', 'none');
end
%% scatter  vertical pairwise
clc
close all
s = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on

for oo = 1;
    x = (1:22)';
y = (mean_pair_acc(:,3,2)+ mean_pair_acc(:,4,2))/2;
dy = (mean_pair_se(:,3,2)+ mean_pair_se(:,4,2))/2;  
fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 0.9 0.9],'linestyle','none');


for subid = 1:sub_n
    scatter(1:22,acc(:,subid,3),80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.3);    
end

plot([x, x]', [mean_pair_acc(:,3,1)-mean_pair_se(:,3,1), mean_pair_acc(:,3,1)+mean_pair_se(:,3,1)]', '-k','LineWidth',4)
scatter(1:22,mean_pair_acc(:,3,1),75,'ko','MarkerFaceColor',[1 1 1],'LineWidth',4);
s.CurrentAxes.XTick = (1:roi_n)+0.25;
s.CurrentAxes.XTickLabel = Hor.roi;
box on
ylim([42 100])
xlim([0.5 22.5])
xtickangle(45)
h=gca; h.XAxis.TickLength = [0 0];
ylabel('Percentage Correct (%)')
title({'Pair-wise decoding accuracy'})
xlabel('ROI')
set(gca,'FontSize',15,'linewidth',1)



%  add toward versus away
s2 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on

for oo = 1;
    x = (1:22)'+0.32;

    for subid = 1:sub_n
        scatter(x,acc(:,subid,4),80,'o','filled','MarkerEdgeColor',[1 1 1],'linewidth',1,'MarkerFaceAlpha',0.3);       
    end
    cp = plot([x, x]', [mean_pair_acc(:,4,1)-mean_pair_se(:,4,1), mean_pair_acc(:,4,1)+mean_pair_se(:,4,1)]','Color',[0.5 0.5 0.5],'LineWidth',4)

    plot([(1:21)'+0.67,(1:21)'+0.67]', [zeros(21,1), 100*ones(21,1)]', '-k','LineWidth',0.1) 
    scatter(x,mean_pair_acc(:,4,1),75,'o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',4);
    s2.CurrentAxes.XTick = (x)+0.25-0.3;
    s2.CurrentAxes.XTickLabel = Hor.roi;
    box on
    ylim([42 100])
    xlim([0.5 22.5])
    xtickangle(45)
    h=gca; h.XAxis.TickLength = [0 0];
    ylabel('Percentage Correct (%)')
    title({'Pair-wise decoding accuracy'})
    xlabel('ROI')
    set(gca,'FontSize',15,'linewidth',1)
    
end
set(gcf, 'color', 'none');   
set(gca, 'color', 'none');
end

%% bias  lateral/depth 



close all
clear all
clc


BASE = '/Users/pw1246/Box Sync/plotf (pw1246@nyu.edu)/';
sub = {'01','04','05','06','28','29'};
%  sub = {'05'};
roi_n = 22;
ses_h = 'ses-0102';
ses_v = 'ses-0304';

acc = zeros(3,sub_n,22,4);
acc_chance = zeros(3,sub_n,22,4);

se = zeros(roi_n,length(sub),2);

for subid = 1:length(sub)
    Hor = load([BASE 'sub-02' sub{subid} '-' ses_h '-classify-bias4.mat']);
   
    for mm = 1:4
        for  rr = 1:22
    acc(:,subid,rr,mm) = Hor.result(rr,:,mm); % acc bias
    acc_chance(:,subid,rr,mm) = Hor.result(rr,:,mm+4); % chance
        end
    end 
end

mean_acc = zeros(roi_n,3,4);
mean_acc_chance = zeros(roi_n,3,4);
mean_se = zeros(roi_n,3,4);
mean_se_chance = zeros(roi_n,3,4);

for rr = 1:22
for mm = 1:4    
    mean_acc(rr,:,mm) = mean(acc(:,:,rr,mm),2);
    mean_acc_chance(rr,:,mm) = mean(acc_chance(:,:,rr,mm),2);   
    mean_se(rr,:,mm) = std(acc(:,:,rr,mm)')/sqrt(sub_n);
    mean_se_chance(rr,:,mm) = std(acc(:,:,rr,mm)')/sqrt(sub_n);
end  
end

all = mean(mean_acc,3);
all_chance = mean(mean_acc_chance,3);
all_se_chance = mean(mean_se_chance,3);
all_se = mean(mean_se,3);
% lateral bias V1 hMT IPS0 
close all

plotwhichroi = [1 9 11];
f2 = figure('Renderer', 'painters', 'Position', [10 10 300 298]);
hold on

% 
% bar1 = bar([1.2 2 2.8],mean_acc(plotwhichroi,1,1),0.2,'linewidth',1);
% bar2 = bar([1.2 2 2.8]+0.18,mean_acc(plotwhichroi,2,1),0.2,'linewidth',1);
% bar3 = bar([1.2 2 2.8]+0.348,mean_acc(plotwhichroi,3,1),0.2,'linewidth',1);

h(1) = bar([1.2 2 2.8],all(plotwhichroi,1),0.2,'linewidth',1);
h(2) = bar([1.2 2 2.8]+0.1735,all(plotwhichroi,2),0.2,'linewidth',1);
h(3) = bar([1.2 2 2.8]+0.348,all(plotwhichroi,3),0.2,'linewidth',1);
xx = [h(1).XData h(2).XData h(3).XData];
% x = [0 sort(xx) 4]';
% y = [100/3 ;reshape(all_chance([1 9 11],:)',9,1); 100/3];
% dy = [0 ;reshape(all_se_chance([1 9 11],:)',9,1) ;0]; 
% 
% f = fill([x;flipud(x)],[y-dy;flipud(y+dy)],[.9 0.9 0.9],'linestyle','none');
% f.FaceAlpha = 0.6;

 h(1).FaceAlpha = 0.7;h(2).FaceAlpha = 0.2;h(3).FaceAlpha = 0.7;
 plot([0 4],[33.3 33.3],'k--','linewidth',1);

f2.CurrentAxes.XTick = [1.2 2 2.8]+0.13;
f2.CurrentAxes.XTickLabel = Hor.roi(plotwhichroi);
ylabel('Percentage decoded as (%)')
set(gca,'FontSize',15)
title({'Oblique decoding bias'})


% 

% l.Position = l.Position + [0.02 0.025 0 0];

 x = sort(xx)';
yy = [reshape(all([1 9 11],:)',9,1);];
ee = [reshape(all_se([1 9 11],:)',9,1)]; 

plot([x, x]', [yy-ee, yy+ee]', '-k','linewidth',1);


xlabel('ROI')
box on
ax = gca;
ax.LineWidth = 2;
ylim([20 53])
xlim([0.7 3.5])

legend(h([1,2,3]),{'2D','Oblique','3D'})
