clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    SECOND-LEVEL DECODING                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script performs a second-level analysis on the prediction accuracy maps to test for differences between conditions.
% 
% 1. Data Setup
%     Load data, create output folders
% 2. Data Preprocessing
%     Normalize and smooth the accuracy maps
% 3. ANOVA Analysis
%     Specify and estimate an ANOVA model based on task groups and covariates
% 4. Contrast Analysis
%     Define contrasts for comparing different conditions
    
%% Data setup

rootFolder = "/.../"

cd(rootFolder + "DecodingAnalysis/DecodingTrialwise")

decoding_data = dir(rootFolder + "DecodingAnalysis/DecodingTrialwise/decoding_trialwise_results")
decoding_data = decoding_data(~ismember({decoding_data.name},{'.','..','._.DS_Store', '.DS_Store'}));
decoding_data = string({decoding_data.name}');

% Get behavioral data
load(rootFolder + "behavioralData.mat");

% Create output folder structure
anovaFolder = rootFolder+"DecodingAnalysis/DecodingTrialwise/decoding_trialwise_anova";

if exist(anovaFolder) == 0
    mkdir(anovaFolder)
end


%% Normalize and smooth all images:

for subject_number = 1:length(decoding_data)
    
    subjectFolder = sprintf('%sDecodingAnalysis/DecodingTrialwise/decoding_trialwise_results/%s/', rootFolder, decoding_data(subject_number));
    accuracy_map= spm_select('FPList', subjectFolder, '^res_accuracy_minus_chance.nii$');
    
    % Normalize:
    matlabbatch = [];
    
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = {sprintf('%spreprocessed_data/%s/1_niftiFiles/anat/y_anat.nii', rootFolder, decoding_data(subject_number))};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = cellstr(accuracy_map);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
        78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'n_';
    spm_jobman('run', matlabbatch)
    
    
    % Smooth:
    n_accuracy_map = spm_select('FPList', subjectFolder, '^n_res_accuracy_minus_chance.nii$')
    matlabbatch = [];
    
    matlabbatch{1}.spm.spatial.smooth.data =  cellstr(n_accuracy_map)
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's_';
    
    spm_jobman('run', matlabbatch)
    
end


%% Run ANOVA

%% Get all normalized and smoothed accuracy maps
for subject_number = 1:length(decoding_data)
    subject_number    
    subjectFolder = sprintf('%sDecodingAnalysis/DecodingTrialwise/decoding_trialwise_results/%s/', rootFolder, decoding_data(subject_number));
    accuracy_maps{subject_number}= spm_select('FPList', subjectFolder, '^s_n_res_accuracy_minus_chance.nii$');
      
end

accuracy_maps = string(vertcat(accuracy_maps{:}));

%% Get covariates, put everything in one big table to sort by condition

% Extract likert responses
for sbjct = 1:height(behavioralData)
focusResponse(sbjct) = mode(behavioralData.response{sbjct}.focusResponse);
instructionsResponse(sbjct) = mode(behavioralData.response{sbjct}.instructionsResponse);
end

% Create anova_data table with relevant columns and sort it
anova_data = behavioralData(:, {'condition','age', 'sex'});
anova_data = [anova_data, table(focusResponse', instructionsResponse', 'VariableNames', {'focusResponse', 'instructionsResponse'}),table(accuracy_maps)]
anova_data = sortrows(anova_data, 'condition')


%% ANOVA with three groups

% Model estimation ANOVA
matlabbatch = []
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(anovaFolder);
 
for cond = 1:3
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(cond).scans = cellstr(anova_data.accuracy_maps(anova_data.condition == cond));
end
 
matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
 
matlabbatch{1}.spm.stats.factorial_design.cov(1).c = anova_data.age;
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(2).c = double(anova_data.sex); % coding: 1=female, 2=male
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(3).c = anova_data.focusResponse;
matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'focusResponse';
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 1;

matlabbatch{1}.spm.stats.factorial_design.cov(4).c = anova_data.instructionsResponse;
matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'instructionsResponse';
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 1;
 
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run', matlabbatch)

% Model estimation
matlabbatch = [];
spmFile= cellstr(spm_select('FPList', anovaFolder, 'SPM.mat'))
 
matlabbatch{1}.spm.stats.fmri_est.spmmat = spmFile;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
 
spm_jobman('run', matlabbatch)

%% Contrasts 

matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {sprintf('%s/SPM.mat', anovaFolder)};

% TESTS WITH COMPARISONS AGAINST OTHER CONDITIONS
% f-test:
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'ftest_comparison';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 -0.5 -0.5; -0.5 1 -0.5; -0.5 -0.5 1];
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';

% t-tests:
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = '1vsall';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1 -0.5 -0.5];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = '2vsall';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [-0.5 1 -0.5];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = '3vsall';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [-0.5 -0.5 1];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

% SIMPLE WITHOUT COMPARING AGAINST OTHER CONDITIONS:

matlabbatch{1}.spm.stats.con.consess{5}.fcon.name = 'ftest_simple';
matlabbatch{1}.spm.stats.con.consess{5}.fcon.weights = [1 0 0; 0 1 0; 0 0 1];
matlabbatch{1}.spm.stats.con.consess{5}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = '1_simple';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [1 0 0];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = '2_simple';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 1 0];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = '3_simple';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 1];
matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

% Pairwise comparisons
matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = '1v3';
matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [1 0 -1];
matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = '3v1';
matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [-1 0 1];
matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = '1v2';
matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [1 -1 ];
matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = '2v1';
matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = '2v3';
matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = [0 1 -1];
matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = '3v2';
matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = [0 -1 1];
matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch)


