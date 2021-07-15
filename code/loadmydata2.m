function [DATA] = loadmydata2(sub,ses,run,BASE,roiname)

DATA = cell(numel(run),numel(roiname));
runidx = 0;

for session = 1:length(ses)
    for runn = 1:numel(run)
        runidx = runidx +1;
        if sub(end) == '4' & session == 1;
            RUN = num2str(runn+4);
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
            if runn == 10;
                RUN = num2str(runn+5);
                datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                    sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
            end
            
        elseif (sub(end) == '3' | sub(end) == '5' | sub(end) == '6') | (sub(end) == '4' & session == 2);
            RUN = num2str(run(runn));
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
        else
            RUN = num2str(run(runn));
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-3dmotion_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
        end
        
        disp(['Loading: ' datapath]);
        Func = niftiread(fullfile(datapath));
        
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
            
            temp_tseries = zeros(numFrames,roiSize);
            baseline = zeros(1,roiSize);
            
            % raw intensity
            for voxel = 1:roiSize
                temp_tseries(:,voxel) = squeeze(Func(x(voxel),y(voxel),z(voxel),:));
            end
            
            
            % percentTseries = zeros(numFrames,roiSize);
            percentTseries = temp_tseries;
            % convert raw to percent change (not needed when detrending)
            %             for voxel = 1:roiSize
            %                 baseline = mean(temp_tseries(:,voxel));
            %                 percentTseries(:,voxel) = 100 * (temp_tseries(:,voxel)/baseline - 1);
            %             end
            
            % fft-based detrend (works less well than linear detrend)
            %             fmriFFT = fft(percentTseries);
            %             fmriFFT(1:5,:) = zeros(5,roiSize);
            %             fmriFFT(end-4:end,:) = zeros(5,roiSize);
            %             temp = real(ifft(fmriFFT));
            
            temp = detrend(percentTseries,1); % linear detrend
            
            % temp = percentTseries-mean(percentTseries,2); % ROI-average based detrend
            % needed for TAFKAP, does not do much beyond linear detrend for
            % classify
            
            temp = normalize(temp); % z-score time-series
            % needed for TAFKAP, does not do much beyond linear detrend for
            % classify
            
            DATA{runidx,roidx} = temp;
            
        end
        
    end
end