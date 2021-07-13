clear all
close all
clc





%% new vs. old fmriprep data

new = load('/Users/pw1246/Desktop/BIDS_6_7/code/result/sub-0201-ses-0304-classify.mat')
old = load('/Users/pw1246/Desktop/BIDS_6_7/code/result/sub-0201-ses-0304-22-mean-classify.mat')

f1 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:length(new.roiname),new.result(:,1),0.4,'k');
bar2 = bar((1:length(new.roiname))+0.4,old.result(:,1),0.4,'k');
bar1.FaceAlpha = 0.05;
bar2.FaceAlpha = 0.3;
legend('new','old')
f1.CurrentAxes.XTick = (1:22)+0.3;
f1.CurrentAxes.XTickLabel = new.roiname;
set(gca,'FontSize',15)
xtickangle(45)


%% matlab classify (sub vs. sub)

sub1 = load('/Users/pw1246/Desktop/motion/code/result/sub-0201-ses-0102-classify.mat')
sub2 = load('/Users/pw1246/Desktop/motion/code/result/sub-0202-ses-0102-classify.mat')

f1 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
bar1 = bar(1:length(sub1.roiname),sub1.result(:,1),0.4,'k');
bar2 = bar((1:length(sub1.roiname))+0.4,sub2.result(:,1),0.4,'k');
bar1.FaceAlpha = 0.05;
bar2.FaceAlpha = 0.3;
legend('sub1','sub2')
f1.CurrentAxes.XTick = (1:22)+0.3;
f1.CurrentAxes.XTickLabel = sub1.roiname;
set(gca,'FontSize',15)
xtickangle(45)

%% Accuracy across ROIs (averaged across subjects)
clear all
close all
clc

sub = {'01','02','29'};
acc = zeros(22,length(sub));
for subid = 1:length(sub)
loadsub = load(['/Users/pw1246/Desktop/motion/code/result/sub-02' sub{subid} '-ses-0304-classify.mat']);
acc(:,subid) = loadsub.result(:,1);
end


meanacc = mean(acc,2);
se = std(acc')/sqrt(length(sub));


f1 = figure('Renderer', 'painters', 'Position', [10 10 950 310]);
hold on
plot(acc,'.','markersize',30,'linewidth',2)
bar1 = bar(1:size(meanacc,1),meanacc,0.4,'k');
errorbar(bar1.XData,bar1.YData,se,'k.','linewidth',2);
%bar2 = bar((1:size(meanacc,1))+0.4,meanacc,0.4,'k');
bar1.FaceAlpha = 0.05;
%bar2.FaceAlpha = 0.3;
%legend('sub1','sub2')
f1.CurrentAxes.XTick = (1:22)+0.3;
f1.CurrentAxes.XTickLabel = loadsub.roiname;
legend(sub)
set(gca,'FontSize',15)
xtickangle(45)



%%


clear all
clc

sub = {'02'};
CNFM = cell(22,length(sub));
for subid = 1:length(sub)
loadsub = load(['/Users/pw1246/Desktop/motion/code/result/sub-02' sub{subid} '-ses-02-classify.mat']);
CNFM(:,subid) = loadsub.cnfm;
end

cnfm = cell(size(loadsub.cnfm,1),1);
    for kk = 1:size(loadsub.cnfm,1);
        sum_fm = zeros(8,8);
        for subid = 1:length(sub)
            temp_fm =  CNFM{kk,subid};
            sum_fm = sum_fm + temp_fm;            
        end
        cnfm{kk,1} = sum_fm./length(sub);
    end
    
    
     figure('Renderer', 'painters', 'Position', [10 10 1150 700]);
     ha1 = tight_subplot(1,4,[-0.13 0.05],[-0.05 -0.068],[0.05 .05]);    
    RoI = [1 2 9 11];
    titleroi = loadsub.roiname(RoI);
    for k1 = 1:length(RoI);
        axes(ha1(k1));
        kk = RoI(k1);

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
%         yticks([1:9]);
%         yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%         xticks([1:9]);
%         xticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
%        
                yticks([1:9]);
        yticklabels(cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)} {char(8594)}]))
        xticks([1:9]);

        xticklabels(cellstr([{char(8600)} {char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)}]))
       
        
        text(0.42:8.42,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);

        set(gca,'YDir','normal')
        ylabel('decoded direction')
        xlabel('presented direction')
        axis square;
        Title = [titleroi{k1}];       
        title(Title);
        set(gca,'FontSize',15)
        
        hold off
    end
