clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  SECOND-LEVEL ANALYSIS                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script performs a second-level analysis, steps include:
% 1. Loading Data and Setting Up Folders
%         Load behavioral data and set up the folder structure
%         Retrieve first-level GLM results from the specified folder.
% 2. Preparing Data for ANOVA
%         Extract likert responses and create a table with relevant columns from behavioral data and contrasts maps from the first-level analysis.
% 3. ANOVA Model Specification & Estimation
% 4. Contrasts
%         Define contrasts of interest to test specific hypotheses about the effects of task condition
    
%% Load data, set up folders
 rootFolder = "/.../"             

% Load behavioral data
load(rootFolder + 'behavioralData.mat');

% Get first level results
firstLevelGLM_dataFolder = rootFolder + 'firstLevelGLM';
cd(firstLevelGLM_dataFolder)
firstLevelGLM_data = dir;
firstLevelGLM_data = firstLevelGLM_data(~ismember({firstLevelGLM_data.name},{'.','..','._.DS_Store', '.DS_Store'}));
firstLevelGLM_data = string({firstLevelGLM_data.name}');

% Create output folder structure
outputFolder = rootFolder+'secondLevelGLM/';
anovaFolder= outputFolder + 'anova';

cd(outputFolder)

if exist(anovaFolder) == 0
    mkdir(anovaFolder)
end


%% Get all contrast maps
for subject_number = 1:length(firstLevelGLM_data)
    subject_number
    subjectFolder = sprintf('%s/%s/', firstLevelGLM_dataFolder, firstLevelGLM_data(subject_number));
     con_map{subject_number}= spm_select('FPList', subjectFolder, 'wcon*');
end

con_map = string(vertcat(con_map{:}));

%% Get covariates, put everything in one big table to sort by condition

% Extract likert responses
for sbjct = 1:height(behavioralData)
focusResponse(sbjct) = mode(behavioralData.response{sbjct}.focusResponse);
instructionsResponse(sbjct) = mode(behavioralData.response{sbjct}.instructionsResponse);
end

% Create anova_data table with relevant columns and sort it
anova_data = behavioralData(:, {'condition','age', 'sex'});
anova_data = [anova_data, table(focusResponse', instructionsResponse', 'VariableNames', {'focusResponse', 'instructionsResponse'}),table(con_map)]
anova_data = sortrows(anova_data, 'condition')



%% ANOVA with three groups

% Model specification
matlabbatch = []
matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(anovaFolder);
 
for cond = 1:3
matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(cond).scans = cellstr(anova_data.con_map(anova_data.condition == cond));
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

matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = '1vsall_reverse';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [-1 0.5 0.5];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = '2vsall_reverse';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0.5 -1 0.5];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = '3vsall_reverse';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0.5 0.5 -1];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

% SIMPLE WITHOUT COMPARING AGAINST OTHER CONDITIONS:

matlabbatch{1}.spm.stats.con.consess{8}.fcon.name = 'ftest_simple';
matlabbatch{1}.spm.stats.con.consess{8}.fcon.weights = [1 0 0; 0 1 0; 0 0 1];
matlabbatch{1}.spm.stats.con.consess{8}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = '1_simple';
matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [1 0 0];
matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = '2_simple';
matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [0 1 0];
matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{11}.tcon.name = '3_simple';
matlabbatch{1}.spm.stats.con.consess{11}.tcon.weights = [0 0 1];
matlabbatch{1}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

% Pairwise comparisons
matlabbatch{1}.spm.stats.con.consess{12}.tcon.name = '1v3';
matlabbatch{1}.spm.stats.con.consess{12}.tcon.weights = [1 0 -1];
matlabbatch{1}.spm.stats.con.consess{12}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{13}.tcon.name = '3v1';
matlabbatch{1}.spm.stats.con.consess{13}.tcon.weights = [-1 0 1];
matlabbatch{1}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{14}.tcon.name = '1v2';
matlabbatch{1}.spm.stats.con.consess{14}.tcon.weights = [1 -1 ];
matlabbatch{1}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{15}.tcon.name = '2v1';
matlabbatch{1}.spm.stats.con.consess{15}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{15}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{16}.tcon.name = '2v3';
matlabbatch{1}.spm.stats.con.consess{16}.tcon.weights = [0 1 -1];
matlabbatch{1}.spm.stats.con.consess{16}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{17}.tcon.name = '3v2';
matlabbatch{1}.spm.stats.con.consess{17}.tcon.weights = [0 -1 1];
matlabbatch{1}.spm.stats.con.consess{17}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 0;

spm_jobman('run', matlabbatch)


