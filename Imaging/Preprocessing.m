clear
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        PREPROCESSING                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs a set of steps over each particpant's fMRI scans to use
% for subsequent univariate analyses. Steps include:

% 1. DICOM to NIFTI Conversion
% 2. Data Organization
%         Renaming and moving files to appropriate folders for better organization.
%         Checking and ensuring the correct number of runs.
% 3. 3D to 4D Conversion
% 4. Field Map Correction
%         Calculating VDM (voxel displacement map) files for field map correction 
% 5. Realignment and Unwarping
%         Realignment and unwarping of functional runs using the calculated VDM images 
% 6. Coregistration
%         Coregistration of functional image and anatomical image 
% 7. Segmentation
%         Anatomical image segmentation into different tissue classes 
% 8. Smoothing
%         Spatial smoothing of functional images using a Gaussian kernel



%% Determine location of files

project_folder    = '.../ProjectFolder'
cd(project_folder)

%Determine number of files:
allFiles = dir(sprintf('%s/RANDECOD_DICOM', project_folder));
visibleFiles = allFiles(arrayfun(@(x) ~strcmp(x.name(1),'.'),allFiles)) %remove hidden files

% Constants
total_subject     = length(visibleFiles); %number of participants.
total_run         = 6; %number of runs
%% Preprocessing steps loop over all participants:


for sbjct = 1:length(visibleFiles) %Loop preprocessing steps over each participant
    tic
    subject_number = visibleFiles(sbjct).name
    sprintf('Starting preprocessing of subject: %d', sbjct)

    %% Housekeeping, accessing raw data folder and creating new output directory

    %select raw dcm files
    dcmFolder = sprintf('%s/RANDECOD_DICOM/%s', project_folder, subject_number)
    cd(dcmFolder)
    dicomFiles = spm_select('FPListRec', dcmFolder, '.dcm');

    %Create new output directory:
    outputFolder = sprintf('%s/preprocessed_data/%s', project_folder, subject_number)
    mkdir(outputFolder)

    %Create directory where converted nifti files go:
    niftiFolder = sprintf('%s/1_niftiFiles', outputFolder)
    mkdir(niftiFolder)
    niftiFolder = convertStringsToChars(niftiFolder)

    %% DICOM TO NIFTI CONVERSION

    %create the batch
    matlabbatch = [];
    matlabbatch{1}.spm.util.import.dicom.data = cellstr(dicomFiles)
    matlabbatch{1}.spm.util.import.dicom.root = 'series';
    matlabbatch{1}.spm.util.import.dicom.outdir = {niftiFolder};
    matlabbatch{1}.spm.util.import.dicom.protfilter = '.*';
    matlabbatch{1}.spm.util.import.dicom.convopts.format = 'nii';
    matlabbatch{1}.spm.util.import.dicom.convopts.meta = 0;
    matlabbatch{1}.spm.util.import.dicom.convopts.icedims = 0;

    %run batch
    spm_jobman('run', matlabbatch)

    %% move files and create appropriate data structure

    fprintf('Move files and create appropriate data structure \n')
    cd(niftiFolder)
    %rename anatomy folder to anat
    movefile(dir('anat*').name, 'anat')
    cd('anat')
    movefile(dir('MF*').name, 'anat.nii')

    cd(niftiFolder)
    %rename localizer folder to localizer
    movefile(dir('loc*').name, 'localizer')

    %CHECK: Are there 6 runs?
    if numel(dir('ep2d*')) ~= 6
        error("NOT CORRECT NUMBER OF RUNS")
    end

    %rename EPIs:
    mkdir('func')
    epis = dir('ep2d*')
    epis = {epis.name}
    epiNames = ["run_01", "run_02", "run_03", "run_04", "run_05", "run_06"]
    for run =1:length(epis)
        movefile(epis{run}, epiNames(run))
        movefile(epiNames(run), 'func')
    end

    %rename fieldnamps
    cd(niftiFolder)
    mkdir('fieldmaps')
    gre = dir('gre_field*')
    gre= {gre.name}
    cd(niftiFolder+"/"+string(gre{1}))
    magMapShortTE = dir('*1-0.nii*').name
    magMapLongTE = dir('*2-0.nii*').name
    %rename and move maps
    movefile(magMapShortTE, niftiFolder+"/fieldmaps/magMapShortTE.nii")
    movefile(magMapLongTE, niftiFolder+"/fieldmaps/magMapLongTE.nii")
    cd(niftiFolder+"/"+string(gre{2}))
    phaseMap = dir('*1-0.nii*').name
    movefile(phaseMap, niftiFolder+"/fieldmaps/phaseMap.nii")
    %remove old folders
    cd(niftiFolder)
    rmdir('gre_*')


    %% Convert 3D scans to 4D
    fprintf('Merging 3D files into 4D\n')

    %select folder with functional runs:
    funcFolder = sprintf('%s/func', niftiFolder)

    % Check if there are correct number of volumes in each run:
    for run = 1:length(epis)
        vols = dir(sprintf('%s/%s',funcFolder, epiNames(run)));
        vols = vols(arrayfun(@(x) ~strcmp(x.name(1),'.'),vols));
        numVols = numel(vols);

        if numVols ~= 230
            error(sprintf("NUMBER OF VOLUMES IN RUN %d IS %d!!",run, numVols))
        end
    end


    %select specific folder:
    for run =1:length(epis)

        funcRun = epiNames(run)
        funcRunFolder =sprintf('%s/%s', funcFolder, funcRun)
        files = spm_select('FPListRec', funcRunFolder, '.nii');
        matlabbatch                       = [];
        matlabbatch{1}.spm.util.cat.vols  = cellstr(files)
        matlabbatch{1}.spm.util.cat.name  = convertStringsToChars(epiNames(run));
        matlabbatch{1}.spm.util.cat.dtype = 0;
        spm_jobman('run', matlabbatch);


        %delete 3d scans
        fprintf('Deleting 3d files\n')
        files        = cellstr(files);
        delete(files{:});

        %move func run out of its folder and delete folders
        cd(funcRunFolder)
        movefile(dir('*.nii').name, funcFolder)
    end
    cd(funcFolder)
    for x = 1:length(epiNames)
        rmdir(epiNames(x), 's')
    end


    %% Extract epi_runs and put in cell for later use

    for run = 1:total_run
        path_epi = sprintf('%s/func/%s.nii', niftiFolder, epiNames(run));
        files = spm_select('expand', path_epi);
        epi_run{run} = cellstr(files);
    end

    %% Field map correction 
    fprintf('Field map correction - Create VDM files')
   
    %Extract phase image:
    magMapShortTE = sprintf('%s/fieldmaps/magMapShortTE.nii', niftiFolder)
    phaseMap = sprintf('%s/fieldmaps/phaseMap.nii', niftiFolder)


    matlabbatch = [];
    %Create vdm image:
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = {phaseMap};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {magMapShortTE};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [4.92 7.38];
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 51.2;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'/home/maja/Downloads/spm12/toolbox/FieldMap/T1.nii'};
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;

    for run = 1:total_run
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session(run).epi = cellstr(epi_run{run}(1));
    end
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 0;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 1;

    spm_jobman('run', matlabbatch);


    %% Realign and unwarp with VDM image:

    matlabbatch = [];

    for run = 1:total_run
        matlabbatch{1}.spm.spatial.realignunwarp.data(run).scans = epi_run{run}
        matlabbatch{1}.spm.spatial.realignunwarp.data(run).pmscan = {sprintf('%s/fieldmaps/vdm5_scphaseMap_session%d.nii', niftiFolder, run)};
    end

    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 7;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';

    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 7;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u_';

    spm_jobman('run', matlabbatch);


    %% Coregistration

    %path to mean picture:
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {sprintf('%s/func/meanu_run_01.nii', niftiFolder)};

    %path to anatomy:
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {sprintf('%s/anat/anat.nii', niftiFolder)};

    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    spm_jobman('run', matlabbatch)

    %% Segmentation
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {convertStringsToChars(niftiFolder+"/anat/anat.nii")};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/home/maja/Downloads/spm12/tpm/TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/home/maja/Downloads/spm12/tpm/TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/home/maja/Downloads/spm12/tpm/TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/home/maja/Downloads/spm12/tpm/TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/home/maja/Downloads/spm12/tpm/TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/home/maja/Downloads/spm12/tpm/TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
        NaN NaN NaN];


    spm_jobman('run', matlabbatch)



    %% Smoothing
    files = spm_select('FPList', funcFolder, '^u_run*')
    matlabbatch = [];

    matlabbatch{1}.spm.spatial.smooth.data =  cellstr(files)
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's_';

    spm_jobman('run', matlabbatch)

end


