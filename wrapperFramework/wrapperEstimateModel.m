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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run the model for all subjects
for i = 1:numel(subjects)
    
    % get the current subject folder name
    subj = subjects{i};
    
    % create the paradigm and features folders
    modelFolder = fullfile(wrapperOptions.outputFolder, wrapperOptions.modelName);
    wrapperOptions.paradigmFolder = fullfile(modelFolder, subj, 'paradigms');
    wrapperOptions.featuresFolder = fullfile(modelFolder, subj, 'features');
    
    % if the model folder doesn't exist, create it
    if ~isdir(modelFolder)
        mkdir(modelFolder)
    end
    if ~isdir(wrapperOptions.paradigmFolder)
        mkdir(wrapperOptions.paradigmFolder)
    end
    if ~isdir(wrapperOptions.featuresFolder)
        mkdir(wrapperOptions.featuresFolder)
    end
    
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
        [features, featureNames, interactionMasks, valCompFeatures] = wrapperCreateFeaturesSimple(subj, ...
            stimIDMapFilename, ...
            'features', ...
            wrapperOptions.featuresFolder, ...
            wrapperOptions.regressorsFolder, ...
            wrapperOptions.regressorFilenames, ...
            wrapperOptions.estReps, ...
            wrapperOptions.valReps, ...
            wrapperOptions.stimTotalThreshold);

        %%%%%%%%%%WARNING HACK!!!%%%%%%%%%%%%%%%%%%%%%%%%%
        % Just assume that the nuissance event is in the first files passed in, and
        % that it indexes the column in that first file. FIX THIS
        excludeValFeatures = wrapperOptions.regressorNuissanceEvents;
    else
        % create the features to use. This returns a matrix with the regressor
        % values for all the stimuli, in the order
        [features, featureNames, interactionMasks, valCompFeatures, excludeValFeatures] = wrapperCreateFeatures(subj, wrapperOptions, stimIDMapFilename);
    end
    
    featuresNuissanceEst = [];
    featureNuissanceNamesEst = [];
    if numel(wrapperOptions.regressorRemoveNuissanceFilenamesFIREst) > 0
        % create a features matrix for nuissance regressors that are to be
        % removed before the model, but have FIR bins
        [featuresNuissanceEst, featureNuissanceNamesEst, ~, ~] = wrapperCreateFeaturesSimple(subj, ...
            stimIDMapFilename, ...
            'featuresNuissanceEst', ...
            wrapperOptions.featuresFolder, ...
            wrapperOptions.regressorsFolder, ...
            wrapperOptions.regressorRemoveNuissanceFilenamesFIREst, ...
            wrapperOptions.estReps, ...
            wrapperOptions.valReps, ...
            wrapperOptions.stimTotalThreshold);
    end
    
    featuresNuissanceVal = [];
    featureNuissanceNamesVal = [];
    if numel(wrapperOptions.regressorRemoveNuissanceFilenamesFIRVal) > 0
        % create a features matrix for nuissance regressors that are to be
        % removed before the model, but have FIR bins
        [featuresNuissanceVal, featureNuissanceNamesVal, ~, ~] = wrapperCreateFeaturesSimple(subj, ...
            stimIDMapFilename, ...
            'featuresNuissanceVal', ...
            wrapperOptions.featuresFolder, ...
            wrapperOptions.regressorsFolder, ...
            wrapperOptions.regressorRemoveNuissanceFilenamesFIRVal, ...
            wrapperOptions.estReps, ...
            wrapperOptions.valReps, ...
            wrapperOptions.stimTotalThreshold);
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

    % transfer the nuissance regressor filenames from the wrapper to the
    % mrifOptions, as that's where they're used
    mrifOptions.removeNuissanceFilenamesEst = wrapperOptions.removeNuissanceFilenamesEst;
    mrifOptions.removeNuissanceFilenamesVal = wrapperOptions.removeNuissanceFilenamesVal;
    mrifOptions.modelNuissanceFilenames = wrapperOptions.modelNuissanceFilenames;

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
    
    % save the options before running the model
    optionsFilename = fullfile(modelFolder, 'options.mat');
    save(optionsFilename, 'wrapperOptions', 'mrifOptions');
    
    % fit the FIR model which uses elastic net and leave one out
    % cross-validation to estimate the weights
    mriFit(experiment, subj, sessionsToAnalyze, features, featuresNuissanceEst, featuresNuissanceVal, mrifOptions, voxIdxs);
    
    % create niftis files of the beta weights and correlation coefficients
    % by concatonating all the model fits
    wrapperCreateOutputFiles(experiment, subj, sessionsToAnalyze, useFeatureNames, interactionMasks, wrapperOptions.writeTopVoxels, mrifOptions);
end