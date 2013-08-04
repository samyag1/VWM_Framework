function wrapperEstimateModel(wrapperOptions, mrifOptions, subjects, clearInput, clearOutput)

% determine which sessions will be analyzed
if wrapperOptions.flattenSessions
    sessionsToAnalyze = 1;
else
    sessionsToAnalyze = 1:wrapperOptions.sessCount;
end

% TODO - This desperately should be done in setter methods inside an
% options class
% calculate the options that are generated automaticaly
wrapperOptions.runCount = numel(wrapperOptions.estRuns) + numel(wrapperOptions.valRuns);
%wrapperOptions.stimCount = wrapperOptions.sessCount*(numel(wrapperOptions.estRuns)*wrapperOptions.estStim + numel(wrapperOptions.valRuns)*wrapperOptions.valStim);

% delete the input files created by this script so they will be remade
if clearInput
    clearInputFiles(wrapperOptions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Paradigms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create the paradigm files, if they're not already created. These contain
% the order of the stimuli as presented
stimIDMapFilename = wrapperCreateParadigms(wrapperOptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run the model for all subjects
for i = 1:numel(subjects)
    
    % get the current subject folder name
    subj = subjects{i};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Experiment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create the experiment
    experiment = wrapperCreateExperiment(subj, wrapperOptions);
    
    if clearOutput
        for curSession = sessionsToAnalyze
            estFolder = fullfile(experiment.(subj){curSession}.estDir, '*');
            delete(estFolder);
            valFolder = fullfile(experiment.(subj){curSession}.valDir, '*');
            delete(valFolder);
            snrFolder = fullfile(experiment.(subj){curSession}.snrDir, '*');
            delete(snrFolder);
            niftisFolder = fullfile(experiment.(subj){curSession}.niftisDir, '*');
            delete(niftisFolder);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if any of the regressor filenames have an sprintf flag (%s) that
    % means that the subject should be added to them
    for curRegressorFilenameIdx = 1:numel(wrapperOptions.regressorFilenames)
        
        curRegressorFilename = wrapperOptions.regressorFilenames{curRegressorFilenameIdx};        
        if not(isempty(strfind(curRegressorFilename, '%s')))
            wrapperOptions.regressorFilenames{curRegressorFilenameIdx} = sprintf(curRegressorFilename, subj);
        end
    end
    
    if wrapperOptions.useSimpleRegressors
        % create the features to use. This returns a matrix with the regressor
        % values for all the stimuli, in the order
        [features, featureNames, interactionMasks, valCompFeatures] = wrapperCreateFeaturesSimple(subj, wrapperOptions, stimIDMapFilename);
        % simple regressors doesn't support nuissance events...
        excludeValFeatures = [];
    else
        % create the features to use. This returns a matrix with the regressor
        % values for all the stimuli, in the order
        [features, featureNames, interactionMasks, valCompFeatures, excludeValFeatures] = wrapperCreateFeatures(subj, wrapperOptions, stimIDMapFilename);
    end
    
    % this vector indicates the indices of the features that will be used
    % to split apart the prediction accuracies. The wrapperOption
    % regressorValComp is read in the CreateFeatures method and the
    % features created using that regressor are this vector
    mrifOptions.valCompFeatures = valCompFeatures;
    
    % store the feature matrix columns that should be exlcluded from
    % validation prediction because they're nuissance regressors
    mrifOptions.excludeValFeatures = excludeValFeatures;
    
    % the feature names will be used in naming the output files, but since
    % the weights estimated will not include the excludeValFeatures, then
    % we must pass in the features minus those excluded
    useFeatureNames = featureNames;
    useFeatureNames(excludeValFeatures) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % only do voxel selection if the options say to
    voxIdxs = [];
    if wrapperOptions.doVoxelSelection
        % choose which voxels will be estimated. This loads a brain mask file
        % and selects the voxels identfied as brain in tht mask
        voxIdxs = wrapperVoxelSelection(wrapperOptions, subj);
    end
    
    % fit the FIR model which uses elastic net and leave one out
    % cross-validation to estimate the weights
    mriFit(experiment, subj, sessionsToAnalyze, features, mrifOptions, voxIdxs);
    
    % create niftis files of the beta weights and correlation coefficients
    % by concatonating all the model fits
    wrapperCreateOutputFiles(experiment, subj, sessionsToAnalyze, useFeatureNames, interactionMasks, mrifOptions);
end