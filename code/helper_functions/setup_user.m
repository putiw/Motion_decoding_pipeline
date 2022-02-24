function projectDir = setup_user(username)

% user specific locations
switch(username)
    case 'm1'
        projectDir=  '';
    case 'm2'
        projectDir=  '/Users/mariamelsahhar/Dropbox (RVL)/ToShare/Decoding';
        githubDir = '/Users/mariamelsahhar/Documents/GitHub';
        freesurferDir = '/Applications/freesurfer/7.2.0';
    case 'Server'
        projectDir=  '/Volumes/Vision/MRI/Decoding';
    case {'Bas', 'bas'}
        projectDir=  '~/Dropbox/MRI/Decoding';
        freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '~/Documents/GitHub';
    case 'Puti'
        projectDir=  '/Users/pw1246/Desktop/MRI/Decoding';
        
       freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '~/Documents/GitHub';
    case 'Omnia'
        projectDir='/Users/omh7815/Documents/Decoding';
        githubDir = '~/Documents/GitHub';
    case {'class'}
        projectDir=  '/Users/nyuad/Dropbox/fMRI/Decoding';
        freesurferDir = '/Applications/freesurfer/7.2.0';
        githubDir = '/Users/nyuad/Dropbox/GitHub';
end
projectDir = char(py.os.path.realpath(py.os.path.expanduser(projectDir))); % convert relative to absolute path

% setup toolboxes
addpath(genpath(fullfile(pwd,'Toolbox')));
addpath(genpath(fullfile(githubDir, 'gifti'))); % https://github.com/gllmflndn/gifti
addpath(genpath(fullfile(githubDir, 'cvncode'))); % https://github.com/cvnlab/cvncode

% freesurfer settings
setenv('FREESURFER_HOME', freesurferDir);
addpath(genpath(fullfile(freesurferDir, 'matlab')));
setenv('SUBJECTS_DIR', [projectDir '/derivatives/freesurfer']); 
