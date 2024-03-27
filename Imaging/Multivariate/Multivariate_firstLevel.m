clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FIRST-LEVEL GLM                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs through a series of steps for each participant's 
% preprocessed fMRI data to conduct a first level univariate GLM analysis 
% using choices as predictors.
% 
% 1. Data Preparation and Organization
%     Load preprocessed data, exclude 5 datasets (see manuscript)
%     Load motion parameters and behavioral data
% 2. Model Specification
%     Specify model for the GLM analysis 
% 3. Model Estimation
%     Estimate the parameters of the model

%% Data preparation
rootFolder = "/.../"

preprocessed_dataFolder = rootFolder + 'preprocessed_data'
cd(preprocessed_dataFolder)

preprocessed_data = dir
preprocessed_data = preprocessed_data(~ismember({preprocessed_data.name},{'.','..','._.DS_Store', '.DS_Store'}));
preprocessed_data = string({preprocessed_data.name}');

% remove exclusions
exclude_list = [10267, 10415, 10460, 10464, 10534];
preprocessed_data = double(setdiff(preprocessed_data, string(exclude_list)));

for subject_number = 1:length(preprocessed_data)

%% Housekeeping, fetch files, create new directories
% Create output directory
outputFolder = sprintf('%sDecodingAnalysis/DecodingTrialwise/firstLevelGLM_trialWise/%d', rootFolder, preprocessed_data(subject_number));
mkdir(outputFolder);

% Get preprocessed functional scans
funcFolder = sprintf('%s/%d/1_niftiFiles/func', preprocessed_dataFolder, preprocessed_data(subject_number));
anatFolder = sprintf('%s/%d/1_niftiFiles/anat', preprocessed_dataFolder, preprocessed_data(subject_number));
scans= cellstr(spm_select('FPList', funcFolder, '^u_run*'));

% Get motion parameters:
motionParameters = cellstr(spm_select('FPList', funcFolder, '.txt'));

% Get behavioral data
load(rootFolder +'behavioralData.mat');

%% Read behavioral data to get trial onsets

clearvars choice_onsetsRIGHT choice_onsetsLEFT
responseData = behavioralData.response{behavioralData.CCNB_ID == preprocessed_data(subject_number)};

responseData.keys_recode = nan(1,length(responseData.keys))';
responseData.keys_recode(responseData.keys == 1) = -1;
responseData.keys_recode(responseData.keys == 2) = -1;
responseData.keys_recode(responseData.keys == 3) = 1;
responseData.keys_recode(responseData.keys == 4) = 1;


for i = 1:6
    choice_onsetsLEFT{i} = responseData.choiceStart(responseData.runNumber == i & responseData.keys_recode == -1);
    choice_onsetsRIGHT{i} = responseData.choiceStart(responseData.runNumber == i & responseData.keys_recode == 1);
end

%% Model specification
matlabbatch = [];

matlabbatch{1}.spm.stats.fmri_spec.dir = {outputFolder};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 1.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

for i =1:6
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans = scans(i);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).name = 'choice_onsetsLEFT';
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).onset = choice_onsetsLEFT{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(1).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(2).name = 'choice_onsetsRIGHT';
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(2).onset = choice_onsetsRIGHT{i};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(2).duration = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond(2).orth = 1;
    
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

%%%%% Model estimation %%%%%
matlabbatch = [];
spmFile= cellstr(spm_select('FPList', outputFolder, 'SPM.mat'))

matlabbatch{1}.spm.stats.fmri_est.spmmat = spmFile;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run', matlabbatch)
end
toc