% Convert NYUAD dcm images to .nii in BIDS format
%
% Adapted from Aaron Cochrane (@wisc.edu)

clear all
PATH = getenv('PATH'); setenv('PATH', [PATH ':/usr/local/bin']); % Set path

%projectDir  = '/Volumes/Vision/MRI/Decoding';
projectDir  = '/Users/pw1246/Desktop/motion';
sub         = '0205'; %'br'; %'hm'; %'ah'; %'rl'; % 'ds'; %
ses         = {'01','02'}; %{'01','02','03','04'}; %'201019a'; %'201020a'; %'160725a'; %'140821a'; % '151106a'; %

%% Run dcm2bids in the shell wrapped in matlab

for ii = 1:length(ses)
    temp = dir(fullfile(projectDir,['sourcedata/sub-', sub '_ses-' ses{ii} '_Br_*']));
    dcmDir = fullfile(temp.folder, temp.name);
    % !dcm2bids_scaffold
    % !dcm2bids_helper -d dcmDir
    config = [projectDir '/code/bids_convert.json'];
    
    system(['dcm2bids -d ' dcmDir ...
        ' -o ' projectDir ...
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
    jsons = dir(fullfile(projectDir, ['sub-' sub], ['ses-' ses{nses}], 'fmap', '*.json'));
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
clc
system(['fmriprep-docker' ...
    ' ' projectDir ...
    ' ' projectDir '/derivatives' ...
    ' participant --participant-label ' sub ...
    ' --fs-license-file /Applications/freesurfer/license.txt' ...
    ' --output-space T1w fsnative MNI152NLin2009cAsym']);

% TODO: Force mriprep LTS version 20.2.x
