function [samples_l samples_r stim_label] = load_surf(BASE,sub,ses,run)

scanID = repmat((1:10)',1,4);
if sub(end) == '4' % this subject uses different scan index
    scanID(:,1) = [(5:13)';15];
    scanID(:,3) = [(2:11)'];
end

% Allocate data
DATA = cell(numel(ses).*numel(run),2);
stim_label=[];
% true trial design matrix [even run, odd run]
stim = [[5:-1:1 8:-1:6]' [4:8 1:3]'];

for whichSession = 1:numel(ses)
    for whichRun = 1:numel(run)

        runidx = scanID(whichRun,str2double(ses(whichSession)));
        RUN = num2str(runidx);

        datapath_l = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{whichSession},'/func/', ...
            sub,'_ses-',ses{whichSession},'_task-3dmotion_run-',RUN,'_space-fsaverage6_hemi-L_bold.func.gii'];
        datapath_r = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{whichSession},'/func/', ...
            sub,'_ses-',ses{whichSession},'_task-3dmotion_run-',RUN,'_space-fsaverage6_hemi-R_bold.func.gii'];
        Func_l = hcp_readfunc(datapath_l);
        Func_r = hcp_readfunc(datapath_r);

        % Drop initial frames to eliminate transients and reach steady state
        framesToDrop = 10;
        Func_l = Func_l(:,framesToDrop+1:end); % Drop n frames
        Func_r = Func_r(:,framesToDrop+1:end); % Drop n frames
        numFrames = size(Func_l,2);


        samples_r = Func_r(:,:)';
        roiSize_r = size(Func_r,1);
        samples_l = Func_l(:,:)';
        roiSize_l = size(Func_l,1);

        %%  detrend + normalize
        %samples = detrend(samples,1); % linear detrend

        % fft-based detrend (works better than linear detrend)
        fmriFFT = fft(samples_l);
        fmriFFT(1:5,:) = zeros(5,roiSize_l);
        fmriFFT(end-4:end,:) = zeros(5,roiSize_l);
        samples_l = real(ifft(fmriFFT));
        fmriFFT = fft(samples_r);
        fmriFFT(1:5,:) = zeros(5,roiSize_r);
        fmriFFT(end-4:end,:) = zeros(5,roiSize_r);
        samples_r = real(ifft(fmriFFT));

        % samples = samples-mean(samples,2); % ROI-average based detrend
        % Used in TAFKAP. Does not do much beyond linear detrend for
        % classify

        samples_l = (samples_l(1:2:end-1,:) + samples_l(2:2:end,:)) ./2; % take average of every 2 TRs
        samples_r = (samples_r(1:2:end-1,:) + samples_r(2:2:end,:)) ./2;
        samples_l = normalize(samples_l); % z-score sample
        samples_r = normalize(samples_r);
        % needed for TAFKAP, does not do much for classify

        % Organize as scan-based or trial-based analysis

        samples_l = squeeze(mean(reshape(samples_l,8,15,[]),2)); % average every 8th datapoint
        samples_r = squeeze(mean(reshape(samples_r,8,15,[]),2));
        % output is downsampled from TRs (240) x voxels to average response per run per direction (8) x voxels
        DATA{numel(run)*(whichSession-1)+whichRun,1} = samples_l;
        DATA{numel(run)*(whichSession-1)+whichRun,2} = samples_r;
   
        stim_label =  [stim_label; stim(:,mod(runidx,2)+1)]; % assign stimulus labels based on odd/even run
    end
end
samples_l = cell2mat(DATA(:,1));
samples_r = cell2mat(DATA(:,2));
end

