
function [DATA_detrend,DATA_detrend_norm,DATA_fft, Func2] = loadmydata_prep(sub,ses,run,BASE,roiname)
DATA_detrend = cell(numel(run),numel(roiname));
DATA_detrend_norm = cell(numel(run),numel(roiname));
DATA_fft = cell(numel(run),numel(roiname));
%github
DATA = cell(numel(run),numel(roiname));
runidx = 0;

for session = 1:length(ses)
    for runn = 1:numel(run(:,session))
        runidx = runidx +1;
        if sub(end) == '4' & session == 1;
            RUN = num2str(runn+1);
            datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
                sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
%             if runn == 10;
%                 RUN = num2str(runn+5);
%                 datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{session},'/func/', ...
%                     sub,'_ses-',ses{session},'_task-TASK_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
%             end
            
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

        Func2{runidx} = Func;
        
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

            % raw intensity
            for voxel = 1:roiSize
                temp_tseries(:,voxel) = squeeze(Func(x(voxel),y(voxel),z(voxel),:));
            end
            
            %%  detrend + normalize
            temp_detrend = detrend(temp_tseries,1); % linear detrend
            % temp = percentTseries-mean(percentTseries,2); % ROI-average based detrend
            % needed for TAFKAP, does not do much beyond linear detrend for
            % classify
            temp_norm = normalize(temp_detrend); % z-score time-series
            % needed for TAFKAP, does not do much beyond linear detrend for
            % classify            
            DATA{runidx,roidx} = temp_norm;       
            
            %% High-pass fft
             percentTseries = zeros(numFrames,roiSize);
             baseline = zeros(1,roiSize);
            % convert raw to percent change (not needed when detrending)
                        for voxel = 1:roiSize
                            baseline = mean(temp_tseries(:,voxel));
                            percentTseries(:,voxel) = 100 * (temp_tseries(:,voxel)/baseline - 1);
                        end           
            % fft-based detrend (works less well than linear detrend)
                        fmriFFT = fft(percentTseries);
                        fmriFFT(1:5,:) = zeros(5,roiSize);
                        fmriFFT(end-4:end,:) = zeros(5,roiSize);
                        tempfft = real(ifft(fmriFFT));
          
            DATA{runidx,roidx} = tempfft;
            
            DATA_detrend{runidx,roidx} = temp_detrend;
            DATA_detrend_norm{runidx,roidx} = temp_norm; 
            DATA_fft{runidx,roidx} = tempfft;
            
        end
        
    end
end