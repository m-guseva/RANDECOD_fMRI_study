clear 
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FIRST-LEVEL GLM                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs through a series of steps for each participant's 
% preprocessed fMRI data to conduct the first level univariate GLM analysis
% using block onsets as predictors.
% 
% 1. Data Preparation and Organization
%     Load preprocessed data, exclude 5 datasets (see manuscript)
%     Load motion parameters and behavioral data
% 2. Model Specification
%     Specify model for the GLM analysis 
% 3. Model Estimation
%     Estimate the parameters of the model
% 4. Contrast Creation
%    Contrasts of interest are defined
% 5. Normalization
%    Normalize contrast maps to the standard Montreal Neurological Institute (MNI) template


%% Data preparation

preprocessed_dataFolder = '/.../preprocessed_data'
cd(preprocessed_dataFolder)

preprocessed_data = dir
preprocessed_data = preprocessed_data(~ismember({preprocessed_data.name},{'.','..','._.DS_Store', '.DS_Store'}));
preprocessed_data = string({preprocessed_data.name}');

% remove excluded datasets
exclude_list = [10267, 10415, 10460, 10464, 10534];
preprocessed_data = double(setdiff(preprocessed_data, string(exclude_list)));

for subject_number = 1:length(preprocessed_data)
%% Housekeeping, fetch files, create new directories
% Create output directory
outputFolder = sprintf('/.../firstLeveLGLM/%d', preprocessed_data(subject_number));
mkdir(outputFolder);

% Get preprocessed functional scans
funcFolder = sprintf('%s/%d/1_niftiFiles/func', preprocessed_dataFolder, preprocessed_data(subject_number));
anatFolder = sprintf('%s/%d/1_niftiFiles/anat', preprocessed_dataFolder, preprocessed_data(subject_number));
scans= cellstr(spm_select('FPList', funcFolder, 's_u_run*'));

% Get motion parameters:
motionParameters = cellstr(spm_select('FPList', funcFolder, '.txt'));

% Load behavioral data
load('/.../behavioralData.mat');

%% Read behavioral data to get trial onsets

clearvars trial_onsets block_onsets
responseData = behavioralData.response{behavioralData.CCNB_ID == preprocessed_data(subject_number)};
for i = 1:6
trial_onsets{i} = responseData.fixationStart(responseData.runNumber == i);
block_onsets{i} = responseData.fixationStart(responseData.trialNumber == 1 & responseData.runNumber == i);
end


%% Model specification
matlabbatch = [];

matlabbatch{1}.spm.stats.fmri_spec.dir = {outputFolder};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

for i =1:6
matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = spm_select('expand', scans(i));
matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).name = 'block_onset';
matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).onset = block_onsets{i};
matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).duration = 90;


matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).orth = 1;

matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = motionParameters(i);
matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
end


matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

spm_jobman('run', matlabbatch)

%% Model estimation 
matlabbatch = [];
spmFile= cellstr(spm_select('FPList', outputFolder, 'SPM.mat'))

matlabbatch{1}.spm.stats.fmri_est.spmmat = spmFile;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch)
%% Contrast
matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {sprintf('%s/SPM.mat', outputFolder)};

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'contrast';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch)

%% Normalization of relevant contrast maps
matlabbatch = [];

matlabbatch{1}.spm.spatial.normalise.write.subj.def = {sprintf('%s/y_anat.nii', anatFolder)};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {sprintf('%s/con_0001.nii,1', outputFolder)};
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run', matlabbatch)
end