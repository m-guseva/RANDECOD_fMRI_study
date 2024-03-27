clear
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Searchlight Decoding Procedure               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script loads the beta maps from the first level GLM (on unsmoothed
% data) and performs a searchlight decoding procedure that trains a model
% to predict button presses on left-out test data, creating prediction
% accuracy maps.

%% Data preparation
rootFolder = "/.../"

load(rootFolder + 'behavioralData.mat')
firstLevelGLM_folder = sprintf('%sDecodingAnalysis/DecodingTrialwise/firstLevelGLM_trialWise', rootFolder)
cd(firstLevelGLM_folder)

firstLevelGLM_folder = dir
firstLevelGLM_folder = firstLevelGLM_folder(~ismember({firstLevelGLM_folder.name},{'.','..','._.DS_Store', '.DS_Store'}));
firstLevelGLM_folder = string({firstLevelGLM_folder.name}');


%% Decoding 
for subject_number = 1:length(firstLevelGLM_folder)
% Set defaults
cfg = decoding_defaults;

% Set the analysis that should be performed 
cfg.analysis = 'searchlight';
cfg.searchlight.radius = 4; 

% Set the output directory where data will be saved
cfg.results.dir = sprintf('%sDecodingAnalysis/DecodingTrialwise/decoding_trialwise_results/%s', rootFolder,firstLevelGLM_folder(subject_number))

beta_loc = sprintf('%sDecodingAnalysis/DecodingTrialwise/firstLevelGLM_trialWise/%s', rootFolder, firstLevelGLM_folder(subject_number))

labelname1 = 'choice_onsetsLEFT'
labelname2 = 'choice_onsetsRIGHT'

cfg.results.output = 'accuracy_minus_chance'

regressor_names = design_from_spm(beta_loc)
cfg = decoding_describe_data(cfg,{labelname1 labelname2},[-1 1],regressor_names,beta_loc);

cfg.design = make_design_cv(cfg);

results = decoding(cfg);
end