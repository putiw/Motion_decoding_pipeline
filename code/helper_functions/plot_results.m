function plot_results(user, parameters, rois, acc)
    % Extract ROI names
%     for mm = 1:length(mask_files)
%        roi_names(mm) = extractBetween(mask_files(mm).name,'sampled_','.');
%     end

    figure; hold on;
    bar(1:length(rois),mean(acc))
    if length(acc) > 1
    er = errorbar(1:length(rois), mean(acc), std(acc)./sqrt(size(acc,1)));
    else
        er = bar(1:length(rois), mean(acc));
    end
    %er.Color = [0 0 0];
    er.LineStyle = 'none';
    line([.5 length(rois)],[.125 .125],'Color','k');
    set(gca,'Xtick',1:length(rois))
    xticklabels(rois)
    xtickangle(45)
    ylabel('Decoding accuracy')
    xlabel('ROI')
    title(['Decoding Accuracy' ' sub-' parameters.subjects{1} ' ses-' parameters.sessions{1} ' ' parameters.detrend_method]);
    % xlim([1 2.5])
    % ylim([0 1])
    set(gca,'FontSize',14)
    saveas(gcf,fullfile(user.output_path, ['Decoding_accuracy' '_sub-' parameters.subjects{1} '_ses-' parameters.sessions{1} '_date-' datestr(now,30) '.pdf']))
end