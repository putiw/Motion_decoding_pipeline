clear all
close  all
clc

setenv('FREESURFER_HOME', '/Applications/freesurfer');


BASE = ['/Users/pw1246/Desktop/motion/derivatives/fmriprep/'];

sub = 'sub-0203';
ses = 'ses-01';
run = 'run-1'; % reference  run 

task = '3dmotion'; % sub 0201 0202 0229
task = 'TASK'; % sub 0203 0204 0205 0206

acq = '';   % sub 0201 0202 0229
acq = 'acq-highres_'; % sub 0203 0204 0205 0206

thr = '0.5';

func_MNI = [BASE sub '/' ses '/func/' sub '_' ses '_task-' task '_' run '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii'];
REF_func = [BASE sub '/' ses '/func/' sub '_' ses '_task-' task '_' run '_space-MNI152NLin2009cAsym_desc-preproc_bold_TEMP_REF.nii.gz'];

system(['mcflirt -in ' func_MNI ' -out ' REF_func ' -refvol 125 -plots']);

REF_downsample = [BASE sub '/ses-01/func/' sub '_' ses '_task-' task '_' run '_space-MNI152NLin2009cAsym_desc-preproc_bold_TEMP_REF_125.nii.gz'];

system(['fslroi ' REF_func ' ' REF_downsample ' 125 1']);

GM_downsample = [BASE sub '/' ses '/anat/' sub '_' ses '_MNI_GM_downsample.nii.gz'];

system(['/Applications/freesurfer/7.1.1/bin/mri_convert -rl ' REF_downsample ' -rt nearest ' BASE sub '/ses-01/anat/' sub '_' ses '_' acq 'space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz ' GM_downsample]);
 
GM_binarize = [BASE sub '/' ses '/anat/' sub '_' ses '_MNI_GM_downsample_bi_' thr '.nii.gz'];
system(['/Applications/freesurfer/7.1.1/bin/mri_binarize --i ' GM_downsample ' --o ' GM_binarize ' --min ' thr]);
