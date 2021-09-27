function plot_confusionM(cnfm,roi,voxelsize,whichroi)
figure('Renderer', 'painters', 'Position', [10 10 250+240*ceil(sqrt(numel(whichroi))) 900]);
ha1 = tight_subplot(4,6,[0.01 0.02],[0.1 0.05],[0.07 .05]);

titleroi = roi(whichroi);
for mm = 1:length(whichroi);
    axes(ha1(mm));
    kk = whichroi(mm);
    
    combine_cfmx = round(cnfm{kk,1}.*100);
    
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
    text(0.55:8.55,1:9,num2str(diag(combine_cfmx),2),'FontSize',15);
    set(gca,'YDir','normal')
%     ylabel('decoded direction')
%     xlabel('presented direction')
    axis square;
    
    if voxelsize == 0
        title(roi{kk});
    else
    Title = [roi{kk} ' (' num2str(voxelsize(kk)) ' voxels)'];
    title(Title);
    end
   box on 
    set(gca,'FontSize',15,'LineWidth',2)
    hold off
end

delete(ha1(length(whichroi)+1:end));
end