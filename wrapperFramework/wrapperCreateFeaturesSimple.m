function [features, featureNames, interactionMasks, valCompFeatures] = wrapperCreateFeaturesSimple(subject, stimIDMapFilename, featureMatName, featuresFolder, regressorsFolder, regressorFilenames, estReps, valReps, stimTotalThreshold)

% make the filename for the features structure
featuresFilename = sprintf('%s_%s.mat', featureMatName, subject);
featuresFilename = fullfile(featuresFolder, featuresFilename);

% if the features file has already been created, then load it and return
if exist(featuresFilename, 'file')
    load(featuresFilename, 'features', 'featureNames', 'interactionMasks', 'valCompFeatures');
    return;
end

% iterate through the regressor filenames loading each one and 
% calculate the level count and start indices for each
featureNames = {};
featureMaps = {};
regressorFilenameCount = numel(regressorFilenames);
for curRegressorIdx = 1:regressorFilenameCount

    % get the current regressor name and filename
    regressorFilename = regressorFilenames{curRegressorIdx};
    regressorFilename = fullfile(regressorsFolder, regressorFilename);
    
    % determine the file type (CSV, excel, mat) and load
    % appropriately
    [folder, title, ext] = fileparts(regressorFilename);
    switch(ext)
        case '.mat'
            % define the variable names in the mat file
            STIM_FILENAMES = 'fNms';
            FEATURES = 'X';
            FEATURE_NAMES = 'featureNames';
            
            % determine the variables in the mat file
            varNames = who('-FILE', regressorFilename);
            assert(ismember(varNames, FEATURES));
            assert(ismember(varNames, STIM_FILENAMES));
            
            % try and load the data, stim filenames and regressor names
            curRegressorFileData = load(regressorFilename, FEATURES);
            curRegressorFileData = curRegressorFileData.(FEATURES);
            curStimFilenames = load(regressorFilename, STIM_FILENAMES);
            curStimFilenames = curStimFilenames.(STIM_FILENAMES);
            
            % convert the data into a cell array, as required below
            curRegressorFileData = mat2cell(curRegressorFileData, ones(size(curRegressorFileData,1),1), ones(size(curRegressorFileData,2),1));

            if ismember(varNames, FEATURE_NAMES)
                curHeaderVals = load(regressorFilename, FEATURE_NAMES);
                curHeaderVals = curHeaderVals.(FEATURE_NAMES);
            else
                curHeaderVals = cell(size(curRegressorFileData, 2),1);
                for curFeature = 1:numel(curHeaderVals)
                    curHeaderVals{curFeature} = sprintf('feature%04i', curFeature);
                end
            end
        case '.csv'
            % load the regressors mapping between stimulus filename and value
            curRegressorFileData = csvimport(regressorFilename, 'noheader', true);
            
            % assume the first row's regressor value is a string that means it's the
            % column header, so store and remove the first row
            assert(ischar(curRegressorFileData{1,2}) && isempty(str2num(curRegressorFileData{1,2})))
            curHeaderVals = curRegressorFileData(1,:);
            curHeaderVals(1) = [];
            curRegressorFileData(1,:) = [];
            curStimFilenames = curRegressorFileData(:,1);
            curRegressorFileData(:,1) = [];
            
            % It is assumed that the regressors take on numeric values, so convert
            % if they're not already
            if ~isnumeric(curRegressorFileData{1,1})
                curRegressorFileData = cellfun( @str2double, curRegressorFileData, 'UniformOutput', false );
            end
        otherwise
            error('Invalid regressor filetype: %s', regressorFilename);
    end
    
    % for each regressor stored in this file, create a map
    for curColumn = 1:numel(curHeaderVals)

        % get the name of the current regressor
        curRegressorName = curHeaderVals{curColumn};
        
        % Populate a mapping between filename and regressor value. This
        % considers the case that different values are assigned to different
        % repetitions of the same stim. When this occurs, the regressor file
        % will contain multiple rows with the same stim filename, and the order
        % they are within the file is the order they'll be used in the design
        % matrix.
        curRegressorMap = containers.Map();
        for curImageIdx = 1:size(curStimFilenames,1)
            
            % get the filename and value of the current stim
            curFilename = curStimFilenames{curImageIdx,1};
            curVal = curRegressorFileData{curImageIdx,curColumn};
            
            % if this filename is already in the map, then append the current
            % value to the list of values already in the map
            if curRegressorMap.isKey(curFilename)
                curVal = [curRegressorMap(curFilename) curVal];
            end
            
            % store the current value in the map
            curRegressorMap(curFilename) = curVal;
        end
        
        % store the regressors map and name
        featureMaps{end+1} = curRegressorMap;
        featureNames{end+1} = curRegressorName;
    end
end

% the features count is the next start index, whether that comes from the
% first regressor loop above (when there are no interactions), or from the
% interaction loop above
featureCount = numel(featureNames);

% read in all the stimID map that maps stimIDs as they are used in the
% paradigms to filenames. This is needed because the features array must be
% ordered by the stimIDs in the paradigm files
load(stimIDMapFilename, 'stimIDMap');
stimIDFilenames = stimIDMap.Keys;

% define the features matrix that stores the regressor values for all the
% stimuli, where the order is determined by the values in the stimOrder
% files
maxReps = max(estReps, valReps);
features = zeros(stimIDMap.Count, featureCount, maxReps);

% iterate through the unique IDs and add a row to the features
% vector for each one.
for curKeyIdx = 1:stimIDMap.Count
    
    % get the filename of the current stim
    curStimFilename = stimIDFilenames{curKeyIdx};
    curStimIdx = stimIDMap(curStimFilename);
    
    % populate all the rep values for the current
    for curRep = 1:maxReps
        
        % now iterate through all the regressors and add them to the features
        % matrix
        for curStimRegressorIdx = 1:featureCount
            
            % get the current filename to regressor value map
            curRegressorMap = featureMaps{curStimRegressorIdx};
            
            % get the current regressor's list of values for this stim
            curRegressorValueList = curRegressorMap(curStimFilename);
            
            % if the list is the max number of reps for this experiment
            % (probably the validation reps), then just store those
            if numel(curRegressorValueList) == maxReps
                curRegressorValue = curRegressorValueList(curRep);
                % otherwise create a list the length of maxReps populated with NaNs
            else
                % and if this regressor is the same for all features (because
                % there's only 1 copy of it in the regressor file), then use it
                curRegressorValue = 0;
                if numel(curRegressorValueList) == 1
                    curRegressorValue = curRegressorValueList;
                    % otherwise it's likely a per trial estimation regressor, so
                    % use a value if there is one for this rep, otherwise let it
                    % defualt to 0
                elseif curRep <= numel(curRegressorValueList)
                    curRegressorValue = curRegressorValueList(curRep);
                end
            end
            
            % otherwise it's continuous, so store the current regressor
            % value in the correct cell of the features matrix. Add 1 to
            % the start index since all continuous regressors in essence
            % have 1 continous level
            features(curStimIdx, curStimRegressorIdx, curRep) = curRegressorValue;
        end
    end
end

% TODO - verify that all features have at least 1 stimulus associated with
% it. If some don't, then remove them from the features matrix and
% featureNames vector, and write out a warning message

% SET THESE TO EMPTY SINCE THESE AREN't HANDLED IN THIS SIMPLE VERSION
% during the calculation of prediction accuracy the pipeline can use
% various levels of a categorical regressor to split the prediction
% accuracies. See if such a regressor has been specified, and determine the
% indices into the features matrix corresponding with the levels it takes
valCompFeatures = [];
interactionMasks = [];

    
% determine if any features have less than the thresholded number of
% stimuli, and remove those
stimTots = sum((features(:,:,1)~=0),1);
belowThresholdIdxs = find(stimTots < stimTotalThreshold);
belowThresholdFeatures = {};
if numel(belowThresholdIdxs) > 0
    fprintf('The following features have less stimuli than the given threshold of %i', stimTotalThreshold);
    
    % store the featurenames that are below thredholds for
    % documentation
    belowThresholdFeatures = featureNames(belowThresholdIdxs);
    
    % write out the features below the threshold
    for curBelowFeature = 1:numel(belowThresholdFeatures)
        fprinf(belowThresholdFeatures{curBelowFeature});
    end
    
    % remove the features from the features matrix
    features(:,belowThresholdIdxs,:) = [];
    featureNames(belowThresholdIdxs) = [];
end

% write out the features matrix
save(featuresFilename, 'features', 'featureNames', 'interactionMasks', 'valCompFeatures', 'belowThresholdFeatures', 'stimIDFilenames');

end