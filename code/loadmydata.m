function [dataset stim_label] = loadmydata(BASE,sub,ses,run,roi)
%
% Load and process (detrend, normalize) data for classification

% subject scan # for each session
scanID = repmat((1:10)',1,4);
if sub == 'sub-0204'; % this subject uses different scan index
    scanID(:,1) = [(5:13)';15];
    scanID(:,3) = [(2:11)'];
end

% scan name 
if sub == 'sub-0203' | sub == 'sub-0204'| sub == 'sub-0205' | sub == 'sub-0206'
    taskname = 'TASK'; % New York data
else
    taskname = '3dmotion'; % Abu Dhabi data
end

% Allocate data
DATA = cell(numel(ses).*numel(run),numel(roi));
stim_label=[];
% true trial design matrix [even run, odd run]
stim = [[5:-1:1 8:-1:6]' [4:8 1:3]'];

for session = 1:length(ses)
    for runn = 1:numel(run)
        
        runidx = scanID(runn,str2double(ses(session)));
        RUN = num2str(runidx);
        
        datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
            sub,'_ses-',ses{session},'_task-',taskname,'_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];        
        
        disp(['Loading: ' datapath]);
        Func = niftiread(fullfile(datapath));
        
        % Drop initial frames to eliminate transients and reach steady state
        framesToDrop = 10;
        Func = Func(:,:,:,framesToDrop+1:end); % Drop n frames
        numFrames = size(Func,4);
        
        % Extract timecourses within the ROIs
        for roidx = 1:numel(roi)
            
            roiPath = [BASE,'derivatives/fmriprep/',sub,'/ses-01/anat/rois/', ...
                sub,'_space-T1w_downsampled_',roi{roidx},'.nii.gz'];
            
            ROI = niftiread(fullfile(roiPath));
            roiSize = length(find(ROI));
            [x y z] = ind2sub(size(ROI),find(ROI));
            
            % Extract raw intensity timeseries
            samples = zeros(numFrames,roiSize);
            for voxel = 1:roiSize
                samples(:,voxel) = squeeze(Func(x(voxel),y(voxel),z(voxel),:));
            end
            
            %%  detrend + normalize
            samples = detrend(samples,1); % linear detrend
            
            % tSeries = tSeries-mean(tSeries,2); % ROI-average based detrend
            % Used in TAFKAP. Does not do much beyond linear detrend for
            % classify
            
            % fft-based detrend (works less well than linear detrend)
            % fmriFFT = fft(tSeries);
            % fmriFFT(1:5,:) = zeros(5,roiSize);
            % fmriFFT(end-4:end,:) = zeros(5,roiSize);
            % tSeries = real(ifft(fmriFFT));
            
            samples = (samples(1:2:end-1,:) + samples(2:2:end,:)) ./2; % take average of every 2 TRs
            
            samples = normalize(samples); % z-score samples
            % needed for TAFKAP, does not do much for classify
            
            % Scan-based analysis (as opposed to trial-based)
            samples = squeeze(mean(reshape(samples,8,15,[]),2)); % average every 8th datapoint
                  
            % output is downsampled from TRs (240) x voxels to average response per run per direction (8) x voxels
            DATA{(session-1)*numel(run)+runn,roidx} = samples;
            
        end
        
        % Return stimulus labels
        stim_label =  [stim_label; stim(:,mod(runidx,2)+1)]; % assign stimulus labels based on odd/even run
    end
    
end

dataset = cell(1,numel(roi));
for whichRoi = 1:numel(roi)
    dataset{whichRoi} = cell2mat(DATA(:,whichRoi));
end

end