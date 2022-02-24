% dcm2bids for all projects
clear all
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/nyuad/opt/anaconda3/bin:/usr/local/bin/']);

%% User defined variables

MRIDir = '/Users/pw1246/Desktop/MRI/';                       % path to MRI directory
currentDir = pwd;
ProjectName = 'Decoding'; % Decoding % Sample_dMRI   % project folder name
TaskName = 'Br_3D'; % 'Br_3D' % 'Prf'                          % task name in the sourcedata dicom folder

subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID



for sub_idx = 7

sub = subject{sub_idx};

ses = {'01','02','03','04'};                                                % session ID
ses = {'01'};

%% Step 1. Run dcm2bids in the shell wrapped in matlab and modify the .json in the fmap folder
projectDir = [MRIDir ProjectName];
tic
for ii = 1:length(ses)

    dcmDir = [projectDir '/sourcedata/' sub '_ses-' ses{ii} '_' TaskName '*/']; % load the directory of the dicoms
    config = [currentDir '/bids_convert.json'];


    % from dicoms to nifti
    system(['dcm2bids -d ' dcmDir ...
        ' -o ' projectDir '/rawdata' ...
        ' -p ' sub ' -s ' ses{ii} ...
        ' -c ' config ' --forceDcm2niix --clobber']);
end

%%
for ii = 1:length(ses)

    % Fix fmap json files by adding run information
    jsons = dir(fullfile(projectDir, '/rawdata/', sub, ['ses-' ses{ii}], 'fmap', '*.json'));
    for ff = 1:length(jsons) % for each json file
        fname = fullfile(jsons(ff).folder, jsons(ff).name);
        str = fileread(fname);
        val = jsondecode(str);

        funcfiles = dir(fullfile(projectDir, '/rawdata/', sub, ['ses-' ses{ii}], 'func', '*.nii.gz'));
        funcfiles = struct2cell(funcfiles);
        funcfiles = strcat(['ses-' ses{ii}],'/func/', funcfiles(1,:)');
        val.IntendedFor = funcfiles;

        str = jsonencode(val);
        % Make the json output fi       le more human readable
        str = strrep(str, ',"', sprintf(',\n"'));
        str = strrep(str, '[{', sprintf('[\n{\n'));
        str = strrep(str, '}]', sprintf('\n}\n]'));

        fid = fopen(fname,'w');
        fwrite(fid,str);
        fclose(fid);

    end

end


end
toc



%% add task name 

subject = {'sub-0201','sub-0202','sub-0204','sub-0205','sub-0206','sub-0228','sub-0229','sub-0248','sub-0903'};                                          % subject ID


for sub_idx = 7

sub = subject{sub_idx};
ses = {'01','02','03','04'};                                                % session ID

for ii = 1:length(ses)

    % Fix func json files by adding task name information
    func_jsons = dir(fullfile(projectDir, '/rawdata/', sub, ['ses-' ses{ii}], 'func', '*.json'));
    for ff = 1:length(func_jsons) % for each json file
        fname = fullfile(func_jsons(ff).folder, func_jsons(ff).name);
        str = fileread(fname);
        val = jsondecode(str);

        val.TaskName = "3dmotion";

        str = jsonencode(val);
        % Make the json output fi       le more human readable
        str = strrep(str, ',"', sprintf(',\n"'));
        str = strrep(str, '[{', sprintf('[\n{\n'));
        str = strrep(str, '}]', sprintf('\n}\n]'));

        fid = fopen(fname,'w');
        fwrite(fid,str);
        fclose(fid);

    end

end


end