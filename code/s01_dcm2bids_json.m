% Convert NYUAD dcm images to .nii in BIDS format
%
% Adapted from Aaron Cochrane (@wisc.edu)
clear all

% Depending on the install method, fmriprep-docker is in different locations
PATH = getenv('PATH'); setenv('PATH', ['/opt/anaconda3/bin:/usr/local/bin:' PATH]); % 

projectDir  = '/Volumes/Vision/MRI/Decoding';
% projectDir  = '~/Desktop/motion';
sub         = '0248'; %; '0903'; %%'0248'; % '203' has an unexpected epi naming convention
ses         = {'01','02','03','04'}; % {'01','02'}; %%'201019a'; %'201020a'; %'160725a'; %'140821a'; % '151106a'; %

%% Run dcm2bids in the shell wrapped in matlab

for ii = 1:length(ses)
    temp = dir(fullfile(projectDir,['sourcedata/sub-', sub '_ses-' ses{ii} '_Br_*']));
    dcmDir = fullfile(temp.folder, temp.name);
    % !dcm2bids_scaffold
    % !dcm2bids_helper -d dcmDir
    % config = [projectDir '/code/bids_convert.json'];
    config = fullfile(pwd, 'bids_convert.json');
    
    system(['dcm2bids -d ' dcmDir ...
        ' -o ' fullfile(projectDir, 'rawdata') ...
        ' -p ' sub ' -s ' ses{ii} ...
        ' -c ' config ' --forceDcm2niix --clobber']);
end

%% Fix fmap json files by adding run information
% dcm2bids will not set the correct intendedFor field in fmap json files
% See https://docs.google.com/document/d/19oWXHbcYH55vZZ6wickjPrLElfaEl4SelRcAHxJwgZU/edit\

% TODO: Use json matlab toolbox - https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab-a-toolbox-to-encode-decode-json-files
% jsonlab replaces all / with \/ so let's hold off for now
% addpath(genpath('~/Documents/GitHub/jsonlab'));

for nses = 1:length(ses)
    jsons = dir(fullfile(projectDir, 'rawdata', ['sub-' sub], ['ses-' ses{nses}], 'fmap', '*.json'));
    for ff = 1:length(jsons)
        fname = fullfile(jsons(ff).folder, jsons(ff).name);
        %     fid = fopen(fname);
        %     str = fread(fid,inf);
        %     str = char(str');
        %     fclose(fid);
        str = fileread(fname);
        val = jsondecode(str);
        
        % jsonlab
        % dat = loadjson(fname);
        
        % TODO: Update to more general case
        num_runs = 10;
        kk = 1;
        for ii = 1:length(val.IntendedFor)
            for jj = 1:num_runs
                % TODO: Update to more general case
                to_add{kk} = insertAfter(val.IntendedFor{ii}, '_task-3dmotion', ['_run-' num2str(jj,'%02.f')]);
                kk = kk+1;
            end
        end
        val.IntendedFor = to_add;
        str = jsonencode(val);
        % Make the json output file more human readable
        str = strrep(str, ',"', sprintf(',\n"'));
        str = strrep(str, '[{', sprintf('[\n{\n'));
        str = strrep(str, '}]', sprintf('\n}\n]'));
        
        fid = fopen(fname,'w');
        fwrite(fid,str);
        fclose(fid);
        
        % savejson('', dat, [fname '.test']);
    end
end

%% Run fmri_prep

% optional: fmriprep version upgrade 
% system('sudo -H pip3 install fmriprep-docker --upgrade');

system(['fmriprep-docker' ...
    ' ' projectDir '/rawdata' ...
    ' ' projectDir '/derivatives' ...
    ' participant --participant-label ' sub ...
    ' --fs-license-file /Applications/freesurfer/license.txt' ...
    ' --output-space T1w MNI152NLin2009cAsym fsnative:den-41k fsaverage:den-41k' ...
    ' --cifti-output' ]); %...
    %' --skip_bids_validation']);
    


% res-native: the original BOLD resolution
% does not seem to work for fsnative:res-native, produces ~230k vertices
% fsaverage:den-41k is equivalent to fsaverage6
% ' --output-space T1w:res-native MNI152NLin2009cAsym:res-native fsnative:res-native fsaverage:den-41k' ...

% surface area of each hemisphere is about 1200 cm^2
% each of our voxels cover an area of ~ .04 cm^2
% so we need a vertex resolution of ~30k/hemisphere
% fsaverage6 (40,960 vertices per hemisphere)
