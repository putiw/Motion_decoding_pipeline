function visualize_fmap(projectDir,whichHemi,fsnativeidx,dotsize)

patch = [projectDir '/derivatives/freesurfer/fsaverage/surf/' lower(whichHemi) 'h.cortex.patch.flat'];
inflate = [projectDir '/derivatives/freesurfer/fsaverage/surf/' lower(whichHemi) 'h.inflated'];
patch = read_patch(patch);
inflate = fs_read_surf(inflate);
vertices = nan(size(inflate.coord))';
vertices(patch.ind+1, :) = [patch.x; patch.y; patch.z]';
scatter(vertices(:,1),vertices(:,2),2,[0.8 0.8 0.8],'filled');
hold on
scatter(vertices(find(fsnativeidx),1),vertices(find(fsnativeidx),2),dotsize,fsnativeidx(find(fsnativeidx),:),'filled');
axis equal
axis off
cc = colormap(flipud(jet));
% caxis([15 30])
cb = colorbar();
cb.TickLabels = cellstr([{char(8594)} {char(8599)} {char(8593)} {char(8598)} {char(8592)} {char(8601)} {char(8595)} {char(8600)}]);
cb.FontSize = 20;
camroll(-90)
title([whichHemi],'FontSize',25)
end

