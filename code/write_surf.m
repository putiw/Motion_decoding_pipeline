% clear all;clc;close all;
% projectDir = setup_user('Puti');
% 
% % subject/session/etc info
% sub = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
% ses = {'01','02'};
% run = [1:10]';
% hemi = {'L','R'};
% scantype = 'fsaverage6';
% detrendcase = 'fft'; % E.g. 'fft', 'linear', 'roi-average', etc
% cutoff = 1/72; % Hz, filter out frequencies below this value
% 
% % for fsaverge 6, full 20 runs of data has the size of 160x40962 per
% % hemisphere per subject
% DATA =load_surf(projectDir,sub,ses,run,scantype,hemi,detrendcase,cutoff);
stim_label = repmat([[4:8 1:3]';[5:-1:1 8:-1:6]'],10,1);
stim_label(:,2) = repmat([[8 1:7]';[1 8:-1:2]'],10,1);
load('DATA.mat')
%%
L_i_all = zeros(9,40962); R_i_all = zeros(9,40962);
L_r_all = zeros(9,40962); R_r_all = zeros(9,40962);
L_x_all = zeros(9,40962); R_x_all = zeros(9,40962);
L_y_all = zeros(9,40962); R_y_all = zeros(9,40962);
for ii = 1:9
    L_samples = DATA(:,:,1,ii);
    R_samples = DATA(:,:,2,ii);
    L_hemi = zeros(8,size(L_samples,2));
    R_hemi = zeros(8,size(R_samples,2));
    for dir = 1:8
        L_hemi(dir,:) = abs(mean(L_samples(stim_label(:,1) == dir,:)));
        R_hemi(dir,:) = abs(mean(R_samples(stim_label(:,2) == dir,:)));
    end
    [L_m L_i] = max(L_hemi);
    L_nanidx = isnan(L_m);
    L_i(L_nanidx) = 0;
    L_m(L_nanidx) = 0;
    [R_m R_i] = max(R_hemi);
    R_nanidx = isnan(R_m);
    R_i(R_nanidx) = 0;
    R_m(R_nanidx) = 0;
    L_i(L_m<=prctile(L_m(~L_nanidx),90)) = 0;
    R_i(R_m<=prctile(R_m(~R_nanidx),90)) = 0;
    L_i_all(ii,:) = L_i;
    R_i_all(ii,:) = R_i;

%     L_r_all(ii,find(L_i_all(ii,:))) = L_hemi(nonzeros(L_i_all(ii,:)),find(L_i_all(ii,:)));
%     R_r_all(ii,find(L_i_all(ii,:))) = R_hemi(nonzeros(R_i_all(ii,:)),find(R_i_all(ii,:)));
%     aaa = nonzeros(L_i_all(ii,:));
% aaa = (aaa-1)*pi/4;
%     L_x_all(ii,:) = cos(aaa)';
%         aaa = nonzeros(R_i_all(ii,:));
% aaa = (aaa-1)*pi/4;
%     R_x_all(ii,:) = cos(aaa)';
%         aaa = nonzeros(L_i_all(ii,:));
% aaa = (aaa-1)*pi/4;
%     L_y_all(ii,:) = sin(aaa)';
%         aaa = nonzeros(R_i_all(ii,:));
% aaa = (aaa-1)*pi/4;
%     R_y_all(ii,:) = sin(aaa)';
end

sub = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID
%%
for SUB = 1;
close all
figure('Renderer', 'painters', 'Position', [10 10 1000 500]);
ha1 = tight_subplot(1,2,[0.01 0.02],[0.1 0.05],[0.07 .05]);
axes(ha1(1));
 fsnativeidx = cvntransfertosubject('fsaverage6','fsaverage',L_i_all(SUB,:)','lh','nearest','orig','orig');
 visualize_fmap(projectDir,'l',fsnativeidx,2);
%  filename = fullfile([pwd '/maps_png/' sub{SUB} '_l.png']);
 titlename = [sub{SUB} ' - L'];
 title(titlename)
%  saveas(gcf,filename)
%  close all
%  f2 = figure;
axes(ha1(2));
 fsnativeidx = cvntransfertosubject('fsaverage6','fsaverage',R_i_all(SUB,:)','lh','nearest','orig','orig');
 visualize_fmap(projectDir,'r',fsnativeidx,2);
  titlename = [sub{SUB} ' - R'];
 title(titlename)
  
 filename = fullfile([pwd '/maps_png/' sub{SUB} '_horizontal.png']);
  saveas(gcf,filename)
end

% %% max direction map
% fs_savemgz(sub, L_i','l.max.mgz', [pwd,'/maps'], 'lh');
% fs_savemgz(sub, R_i', 'r.max.mgz', [pwd,'/maps'], 'rh');
% 
% %% response map to each direction
% fs_savemgz(sub, L_hemi(1,:)','l.1.mgz', [pwd,'/maps'], 'lh');
% fs_savemgz(sub, R_hemi(1,:)', 'r.1.mgz', [pwd,'/maps'], 'rh');
% 
% %%
% fv_surf