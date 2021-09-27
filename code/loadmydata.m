function [dataset stim_label] = loadmydata(BASE,sub,ses,run,roi,p)
%
% Load and process (detrend, normalize) data for classification
%
% Ultimately we want to return dataset as a cell array: sub x ses x run x
% roi

% subject scan # for each session
scanID = repmat((1:10)',1,4);
if sub(end) == '4' % this subject uses different scan index
    scanID(:,1) = [(5:13)';15];
    scanID(:,3) = [(2:11)'];
end

% scan name

switch sub
    case {'0203ny'}
        taskname = '3dmotion';
    case {'0204'}
        taskname = '3dmotion';
    case {'0205'}
        taskname = '3dmotion';
    case {'0206'}
        taskname = '3dmotion';
    otherwise
        taskname = 'TASK';
end

           
% Allocate data
DATA = cell(numel(ses).*numel(run),numel(roi));
stim_label=[];
% true trial design matrix [even run, odd run]
stim = [[5:-1:1 8:-1:6]' [4:8 1:3]'];

for whichSession = 1:numel(ses)
    for whichRun = 1:numel(run)
        
        runidx = scanID(whichRun,str2double(ses(whichSession)));
        RUN = num2str(runidx);
        
        datapath = [BASE,'derivatives/fmriprep/',sub,'/ses-',ses{whichSession},'/func/', ...
            sub,'_ses-',ses{whichSession},'_task-',taskname,'_run-',RUN,'_space-T1w_desc-preproc_bold.nii.gz'];
        
        disp(['Loading: ' datapath]);
        Func = niftiread(fullfile(datapath));
        
        switch p.sample_unit
            case {'raw'}
                
                dataset{(whichSession-1)*numel(run)+whichRun} = Func; % and done
            
            otherwise 
                
                % Drop initial frames to eliminate transients and reach steady state
                framesToDrop = 10;
                Func = Func(:,:,:,framesToDrop+1:end); % Drop n frames
                numFrames = size(Func,4);
                
                % Extract timecourses within the ROIs
                for whichRoi = 1:numel(roi)
                    
                    roiPath = [BASE,'derivatives/fmriprep/',sub,'/ses-01/anat/rois/', ...
                        sub,'_space-T1w_downsampled_',roi{whichRoi},'.nii.gz'];
                    
                    ROI = niftiread(fullfile(roiPath));
                    roiSize = length(find(ROI));
                    [x y z] = ind2sub(size(ROI),find(ROI));
                    
                    % Extract raw intensity timeseries
                    samples = zeros(numFrames,roiSize);
                    for voxel = 1:roiSize
                        samples(:,voxel) = squeeze(Func(x(voxel),y(voxel),z(voxel),:));
                    end
                    
                    %%  detrend + normalize
                    %samples = detrend(samples,1); % linear detrend
                    
                    % fft-based detrend (works better than linear detrend)
                    fmriFFT = fft(samples);
                    fmriFFT(1:5,:) = zeros(5,roiSize);
                    fmriFFT(end-4:end,:) = zeros(5,roiSize);
                    samples = real(ifft(fmriFFT));
                    
                    % samples = samples-mean(samples,2); % ROI-average based detrend
                    % Used in TAFKAP. Does not do much beyond linear detrend for
                    % classify
                    
                    samples = (samples(1:2:end-1,:) + samples(2:2:end,:)) ./2; % take average of every 2 TRs
                    
                    samples = normalize(samples); % z-score samples
                    % needed for TAFKAP, does not do much for classify
                    
                    % Organize as scan-based or trial-based analysis
                    switch p.sample_unit
                        case {'scan'}
                            samples = squeeze(mean(reshape(samples,8,15,[]),2)); % average every 8th datapoint
                        case {'trial'}
                            % do nothing, keep samples
                        otherwise
                            error('Unknown sample unit')
                    end
                    
                    % output is downsampled from TRs (240) x voxels to average response per run per direction (8) x voxels
                    DATA{numel(run)*(whichSession-1)+whichRun,whichRoi} = samples;
                end % end of roi 
        end % end of run
        
        % Return stimulus labels
        switch p.sample_unit
            case {'scan'}
                stim_label =  [stim_label; stim(:,mod(runidx,2)+1)]; % assign stimulus labels based on odd/even run
            case {'trial'}
                nRepeats = 15; % repeats of nDirections per scan
                stim_label =  [stim_label; repmat(stim(:,mod(runidx,2)+1),nRepeats,1)]; % assign stimulus labels based on odd/even run
            case {'raw'}
                % do nothing
            otherwise
                error('Unknown sample unit')
        end
        
    end
    
end

dataset = cell(1,numel(roi));
for whichRoi = 1:numel(roi)
    dataset{whichRoi} = cell2mat(DATA(:,whichRoi));
end

end