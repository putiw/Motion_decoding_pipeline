function DATA = loadmydata(BASE,sub,ses,run,roiname)
%
% Load and process (detrend, normalize) data for classification

% Allocate data
DATA = cell(numel(ses).*numel(run),numel(roiname));

runidx = 0;
for session = 1:length(ses)
    for runn = 1:numel(run)
        runidx = runidx+1;
        if sub(end) == '4' && session == 1
            RUN = num2str(runn+1);
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
            %             if runn == 10;
            %                 RUN = num2str(runn+5);
            %                 datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
            %                     sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
            %             end
            
        elseif (sub(end) == '3' || sub(end) == '5' || sub(end) == '6') || (sub(end) == '4' && session == 2)
            RUN = num2str(run(runn));
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
            
        else % all other subjects/sessions
            RUN = num2str(run(runn));
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-3dmotion_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
        end
        
        disp(['Loading: ' datapath]);
        Func = niftiread(fullfile(datapath));
             
        % Drop initial frames to eliminate transients and reach steady state
        framesToDrop = 10;
        Func = Func(:,:,:,framesToDrop+1:end); % Drop n frames
        numFrames = size(Func,4);
        
        % Extract timecourses within the ROIs
        for roidx = 1:numel(roiname)
            
            roiPath = [BASE,'derivatives/fmriprep/',sub,'/ses-01/anat/rois/', ...
                sub,'_space-T1w_downsampled_',roiname{roidx},'.nii.gz'];
            
            roi = niftiread(fullfile(roiPath));
            roiSize = length(find(roi));
            [x y z] = ind2sub(size(roi),find(roi));
            
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
            
            % TODO: Return stimulus labels
            
            % output is downsampled from TRs (240) x voxels to average response per run per direction (8) x voxels
            DATA{runidx,roidx} = samples;
            
        end   
    end
end