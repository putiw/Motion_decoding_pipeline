function DATA = loadmydata_prep(sub,ses,run,BASE,roiname)
%
% Load and process (detrend, normalize) data for classification

% Allocate data
DATA = cell(numel(run),numel(roiname));

runidx = 0;
for session = 1:length(ses)
    for runn = 1:numel(run(:,session))
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
            
            tSeries = zeros(numFrames,roiSize);
            
            % raw intensity
            for voxel = 1:roiSize
                tSeries(:,voxel) = squeeze(Func(x(voxel),y(voxel),z(voxel),:));
            end
            
            %%  detrend + normalize
            tSeries = detrend(tSeries,1); % linear detrend
            
            % tSeries = tSeries-mean(tSeries,2); % ROI-average based detrend
            % needed for TAFKAP, does not do much beyond linear detrend for
            % classify
            
            % fft-based detrend (works less well than linear detrend)
            %             fmriFFT = fft(tSeries);
            %             fmriFFT(1:5,:) = zeros(5,roiSize);
            %             fmriFFT(end-4:end,:) = zeros(5,roiSize);
            %             tSeries = real(ifft(fmriFFT));
            
            tSeries = normalize(tSeries); % z-score time-series
            % needed for TAFKAP, does not do much beyond linear detrend for
            % classify
            
            DATA{runidx,roidx} = tSeries;
            
        end   
    end
end