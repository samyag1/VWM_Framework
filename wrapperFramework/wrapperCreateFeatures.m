function [features, featureNames, interactionMasks, valCompFeatures, excludeValFeatures] = wrapperCreateFeatures(subject, options, stimIDMapFilename)

% make the filename for the features structure
featuresFilename = sprintf('features_%s.mat', subject);
featuresFilename = fullfile(options.featuresFolder, featuresFilename);

% if the features file has already been created, then load it and return
if exist(featuresFilename, 'file')
    load(featuresFilename, 'features', 'featureNames', 'interactionMasks', 'valCompFeatures', 'excludeValFeatures');
    return;
end

% make sure all the regresssor fields are the same length
assert(numel(options.regressorFilenames) == numel(options.regressorCategorical));

% deteremine the number of regressors, events and interactions
standardRegressorCount = numel(options.regressorFilenames);
usedRegressors = [options.regressorEvents options.regressorNuissanceEvents];
usedRegressors = unique(usedRegressors);
usedRegressorCount = numel(usedRegressors);
interactionRegressorCount = numel(options.regressorInteractions);
totalRegressorCount = usedRegressorCount + interactionRegressorCount;

% iterate through the regressor filenames loading each one and 
% calculate the level count and start indices for each
regressorMaps = cell(1,standardRegressorCount);
regressorLevelCount = zeros(1,standardRegressorCount);
for curRegressorIdx = 1:standardRegressorCount

    % get the current regressor name and filename
    regressorFilename = options.regressorFilenames{curRegressorIdx};
    regressorFilename = fullfile(options.regressorsFolder, regressorFilename);
    
    % TODO - determine the file type (CSV, excel, mat) and load
    % appropriately
    
    % load the regressors mapping between stimulus filename and value
    curRegressorCSV = csvimport(regressorFilename, 'noheader', true);
    % if the first row's regressor value is a string that means it's the
    % column header, so remove the first row
    if ischar(curRegressorCSV{1,2}) && isempty(str2num(curRegressorCSV{1,2}))
        curRegressorCSV(1,:) = [];
    end
    % It is assumed that the regressors take on numeric values, so convert
    % if they're not already
    if ~isnumeric(curRegressorCSV{1,2})
      curRegressorCSV(:,2) = cellfun( @str2double, curRegressorCSV(:,2), 'UniformOutput', false );        
    end
    
    % Populate a mapping between filename and regressor value. This
    % considers the case that different values are assigned to different
    % repetitions of the same stim. When this occurs, the regressor file
    % will contain multiple rows with the same stim filename, and the order
    % they are within the file is the order they'll be used in the design
    % matrix.
    curRegressorMap = containers.Map();
    for curImageIdx = 1:size(curRegressorCSV,1)
        
        % get the filename and value of the current stim
        curFilename = curRegressorCSV{curImageIdx,1};
        curVal = curRegressorCSV{curImageIdx,2};
        
        % TODO - This dictionary is SLOW, it just does a linear search.
        % Find a dictionary that uses hashing.
        
        % if this filename is already in the map, then append the current
        % value to the list of values already in the map
        if curRegressorMap.isKey(curFilename)
            curVal = [curRegressorMap(curFilename) curVal];
        end
        
        % store the current value in the map
        curRegressorMap(curFilename) = curVal;
    end
    
    % store the regressors map
    regressorMaps{curRegressorIdx} = curRegressorMap;
    
    % figure out how many levels this regressor file has
    if options.regressorCategorical(curRegressorIdx)
        % figure out how many unique values this regressor takes. If the
        % value is zero interpret that to mean no value should be assigned
        % for the stim
        uniqueRegressorVals = unique(cell2mat(curRegressorMap.values));
        uniqueRegressorVals(uniqueRegressorVals==0) = [];
        regressorValsCount = numel(uniqueRegressorVals);
        regressorNameCount = numel(options.regressorCategoricalNames{curRegressorIdx});
        
        % if this fails it means that there are more levels of the
        % regressor in the map file than defined in the names
        assert(regressorValsCount <= regressorNameCount);
        
        % store that count in an array
        regressorLevelCount(curRegressorIdx) = regressorNameCount;
    else
        % otherwise it's not categorical, so it only takes one level
        regressorLevelCount(curRegressorIdx) = 1;        
    end
end

% this is a list of the feature columns to exclude from validation
% prediction because they are nuissance regressors
excludeValFeatures = [];

% now iterate through all the events, which are just the regressors that
% are to be added as features, not in interactions, and add each one's
% name, featureNames and start indices to their repsective lists
regressorStartIndices = zeros(1,totalRegressorCount);
regressorNames = cell(1,totalRegressorCount);
featureNames = {};
for curEventIdx = 1:usedRegressorCount
    
    % get the index of the regressor for the current event
    eventRegressorIdx = usedRegressors(curEventIdx);
    
    % use a helper function to recursively create the feature names for the
    % current regressor which will add all category names if it's a
    % categorical regressor, then add them to the main featureNames list
    curEventNames = createFeatureNames([eventRegressorIdx], options);
    featureNames = [featureNames curEventNames];
    
    % store the name of the current regressor for documentation
    regressorNames{curEventIdx} = options.regressorNames(eventRegressorIdx);
    
    % if there are more regressors than the currentone, then update the start indice for
    % the next regressor, which is just the current start index pluse the
    % level count
    nextStartIdx = regressorStartIndices(curEventIdx) + regressorLevelCount(eventRegressorIdx);
    if curEventIdx < totalRegressorCount
        regressorStartIndices(curEventIdx+1) = nextStartIdx;
    end
    
    % if the current used regressor is a nuissance regressor, then add it
    % to the valExlcudeFeatures list so it isn't used in predictions
    if ismember(options.regressorNuissanceEvents,eventRegressorIdx)
        curRegressorFeatures = [regressorStartIndices(curEventIdx)+1:regressorStartIndices(curEventIdx)+regressorLevelCount(eventRegressorIdx)];
        excludeValFeatures = [excludeValFeatures curRegressorFeatures];
    end
end

% during the calculation of prediction accuracy the pipeline can use
% various levels of a categorical regressor to split the prediction
% accuracies. See if such a regressor has been specified, and determine the
% indices into the features matrix corresponding with the levels it takes
valCompFeatures = [];
if not(isempty(options.regressorValComp))
    
    % find the index of the regressor to use for the validation comparison
    valCompIdx = find(strcmp(options.regressorNames, options.regressorValComp));
    if isempty(valCompIdx)
        error('Non-existant regressorValComp specified: %s. Should be in the regressorName cell array.', options.regressorValComp);
    end
    
    if not(ismember(options.regressorEvents,valCompIdx))
        error('regressorValComp: %s is not specified in regressorEvents when it needs to be. Either change the regressorValComp or add it to the regressorEvents vector.', valCompidx);
    end
    
    % only categorical regressors count
    if options.regressorCategorical(valCompIdx) == true
    
        valCompStart = regressorStartIndices(valCompIdx) + 1;
        valCompEnd = valCompStart + numel(options.regressorCategoricalNames{valCompIdx})-1;
        valCompFeatures = valCompStart:valCompEnd;
    else
        error('A regressorValComp: %s was set that isnt a categegorical regressor. Invalid.', options.regressorValComp);
    end
end

% iterate through the interactions and determine the start indices for each
% of them
featureCount = 0;
interactionLeveCount = zeros(1,interactionRegressorCount);
for curInteractionIdx = 1:interactionRegressorCount
    
    % get the current interaction list
    curInteractionList = options.regressorInteractions{curInteractionIdx};
    
    % iterate through the regressors in this interaction and multiple each
    % one's level count to determine how many columns the current
    % interaction will need
    curInteractionLevelCount = 1;
    curInteractionName = '';
    for curRegressorIdx = 1:numel(curInteractionList)
        
        % get the current regressor ID 
        curID = curInteractionList(curRegressorIdx);
        curName = options.regressorNames{curID};
        
        if curRegressorIdx == 1
            curInteractionName = curName;
        else
            curInteractionName = [curInteractionName ' x ' curName];
        end
        
        % get the level count
        curLevelCount = regressorLevelCount(curID);
        
        % multiple it by the count for this interaction term
        curInteractionLevelCount = curInteractionLevelCount * curLevelCount;
    end
    
    % store the name of the current regressor for documentation
    regressorNames{usedRegressorCount+curInteractionIdx} = curInteractionName;
    
    % store the number of features in the current interaction
    interactionLevelCount(curInteractionIdx) = curInteractionLevelCount;
    
    % use a helper function to recursively create the feature names for all
    % the interaction terms of this interaction, then add them to the main
    % featureNames list
    curInteractionNames = createFeatureNames(curInteractionList, options);
    featureNames = [featureNames curInteractionNames];
    
    % if this isn't the last interaction term, then update the start indice for
    % the next interaction term, which is just the current start index pluse the
    % level count
    totalRegressorIdx = usedRegressorCount + curInteractionIdx;
    nextStartIdx = regressorStartIndices(totalRegressorIdx) + curInteractionLevelCount;
    if curInteractionIdx < interactionRegressorCount
        regressorStartIndices(totalRegressorIdx+1) = nextStartIdx;
    end
end

% the features count is the next start index, whether that comes from the
% first regressor loop above (when there are no interactions), or from the
% interaction loop above
featureCount = nextStartIdx;

% read in all the stimID map that maps stimIDs as they are used in the
% paradigms to filenames. This is needed because the features array must be
% ordered by the stimIDs in the paradigm files
load(stimIDMapFilename, 'stimIDMap');
stimIDFilenames = stimIDMap.Keys;
stimIDFilenames = sort(stimIDFilenames);

% define the features matrix that stores the regressor values for all the
% stimuli, where the order is determined by the values in the stimOrder
% files
maxReps = max(options.estReps, options.valReps);
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
        for curStimRegressorIdx = 1:numel(usedRegressors)
            
            % some of these lists contain values in the order of used
            % regressors (the start indices) and others in the order of main
            % events like the categorical list.
            curStimEventIdx = usedRegressors(curStimRegressorIdx);
            
            % get the current filename to regressor value map
            curRegressorMap = regressorMaps{curStimRegressorIdx};
            
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
            
            % get the start index in the features matrix for the current
            % regressor
            curStartIndex = regressorStartIndices(curStimRegressorIdx);
            
            % if the current regressor is categorical then we'll put a 1 in the
            % level of that regressor that is repressented in the map
            if options.regressorCategorical(curStimEventIdx)
                % if the regressor is categorical then a zero value means
                % there is no regressor for the current stim, so leave it
                % alone
                if curRegressorValue ~= 0
                    features(curStimIdx, curStartIndex + curRegressorValue, curRep) = 1;
                end
            else
                % otherwise it's continuous, so store the current regressor
                % value in the correct cell of the features matrix. Add 1 to
                % the start index since all continuous regressors in essence
                % have 1 continous level
                features(curStimIdx, curStartIndex + 1, curRep) = curRegressorValue;
            end
        end
        
        % now add all the interaction term values for the current stim
        for curStimInteractionIdx = 1:interactionRegressorCount
            
            % get the current interaction list which contains the standard
            % regressors IDs to use in the interaction
            curInteractionList = options.regressorInteractions{curStimInteractionIdx};
            
            % iterate through all the terms IN REVERSE to determine which column to update
            featureInteractionValue = 1;
            featureInteractionIdx = 0;
            prevLevelCount = 1;
            for curStandardIDIdx = numel(curInteractionList):-1:1
                
                % get the standard regressor ID for the current regressor in
                % the interaction term
                curStandardID = curInteractionList(curStandardIDIdx);
                
                % get the current filename to regressor value map
                curStandardMap = regressorMaps{curStandardID};
                
                % get the current regressor's value for this stim
                curStandardValueList = curStandardMap(curStimFilename);
                
                % if the list is the max number of reps for this experiment
                % (probably the validation reps), then just store those
                if numel(curStandardValueList) == maxReps
                    curStandardValue = curStandardValueList(curRep);
                    % otherwise create a list the length of maxReps populated with NaNs
                else
                    % and if this regressor is the same for all features (because
                    % there's only 1 copy of it in the regressor file), then use it
                    curStandardValue = 0;
                    if numel(curStandardValueList) == 1
                        curStandardValue = curStandardValueList;
                        % otherwise it's likely a per trial estimation regressor, so
                        % use a value if there is one for this rep, otherwise let it
                        % defualt to nan
                    elseif curRep <= numel(curStandardValueList)
                        curStandardValue = curStandardValueList(curRep);
                    end
                end
                
                % if the regressor is categorical, then the index into the
                % features matrix must be updated
                if options.regressorCategorical(curStandardID)
                
                    % if any categorical value in the interaction has a 0
                    % value that means to exclude this feature of this
                    % stim, so set the interaction value to -1 so it is
                    % ignored later and break out
                    if curStandardValue == 0
                        featureInteractionValue = nan;
                        break;
                    end
                    % the index into the features matrix for this interaction
                    % term is calculated by multiplying the value of the
                    % current categorical value by the previous (i.e. higher
                    % since we're moving in reverse) categorical
                    % levels for this interaction, and then summing all these
                    % values. For example, if there are 2 categorical
                    % regressors (A qand B), A with 20 levels and B with 3, then
                    % their the order will be: A1-B1, A1-B2,
                    % A1-B3, ...,A20-B1,A20-B2,A20-B3.
                    featureInteractionIdx = featureInteractionIdx + (prevLevelCount * (curStandardValue-1));
                    curLevelCount = regressorLevelCount(curStandardID);
                    prevLevelCount = prevLevelCount * curLevelCount;
                    % otherwise it's continuous, so the value to be put into the
                    % matrix must be updated
                else
                    featureInteractionValue = featureInteractionValue * curStandardValue;
                end
            end
            
            % add the starting index for this interaction, plus one since matlab is one based indexing
            curInteractionStartIndex = regressorStartIndices(curStimInteractionIdx+usedRegressorCount);
            featureInteractionIdx = featureInteractionIdx + curInteractionStartIndex + 1;
            
            % update the features matrix if there's a valid interaction
            % value
            if not(isnan(featureInteractionValue))
                features(curStimIdx, featureInteractionIdx, curRep) = featureInteractionValue;
            end
        end
    end
end

% now calculate the interaction masks. if regressors are not made into events, but only interaction terms then
% the original event can be pulled out of the interaction term by adding
% the features of all levels of the interaction. To do this, specify a 2
% element vector, with the first value indicating the interaction to use,
% and the second value indicating the term within that interactiono to use.
% If the term specified is not a categorical variable then all the values
% of that interaction term will be added.
maskCount = numel(options.regressorInteractionMasks);
interactionMasks = {};
for curMask = 1:maskCount
    
    % get the current mask specifications
    curMaskSpecs = options.regressorInteractionMasks{curMask};
    
    % get the index of the event to pull out of the interaction
    curInteractionIdx = curMaskSpecs(1);
    curInteractionList = options.regressorInteractions{curInteractionIdx};
    curInteractionTermIdx = curMaskSpecs(2);
    curInteractionEventIdx = curInteractionList(curInteractionTermIdx);
    
    % determine the location of the current interaction's features in the
    % feature matrix. Add one since it's 0 indexed
    curInteractionStartIdx = regressorStartIndices(usedRegressorCount+curInteractionIdx) + 1;
    curInteractionFeatureCount = interactionLevelCount(curInteractionIdx);
    
    % get the name of the event for the mask
    curEventName = options.regressorNames{curInteractionEventIdx};
    
    % if the mask event is categorical then choose the features which
    % comprise that event
    if options.regressorCategorical(curInteractionEventIdx)

        % get the names of the different levels of the event to be masked
        curEventLevelNames = options.regressorCategoricalNames{curInteractionEventIdx};

        % determine how many levels the event has
        curEventLevelCount = regressorLevelCount(curInteractionEventIdx);

        % can't really explain this other than it works. These are used in
        % the algorithm to determine which features include the given level
        % of the event
        chunkCount = prod(regressorLevelCount(curInteractionList(1:(curInteractionTermIdx-1))));
        chunkSize = prod(regressorLevelCount(curInteractionList((curInteractionTermIdx+1):end)));
        
        % iterate through all the levels of the categorical event and make
        % a mask for each one
        for curEventLevel = 1:curEventLevelCount
        
            % initialize the current mask
            curInteractionMask = [];
            curInteractionMask.mask = zeros(featureCount,1);
            
            % create and store the name of the current level's mask
            curFeatureName = curEventLevelNames{curEventLevel};
            curInteractionMask.name = [curEventName '-' curFeatureName];

            % algorithm to determine which features include the current
            % level of the event. It will set them all to 1 in the mask
            for curChunk = 1:chunkCount
                curChunkStart = ((curChunk-1)*chunkSize*curEventLevelCount + 1) + ((curEventLevel-1)*chunkSize);
                curChunkEnd = curChunkStart + chunkSize -1;
                curInteractionMask.mask(curChunkStart:curChunkEnd) = 1;
            end
            
            % add the current interaction mask to the list
            interactionMasks{end+1} = curInteractionMask;
        end
    % otherwise it's scalar, so include all the items of that interaction
    else
        % initialize the current mask
        curInteractionMask = [];
        
        maskStart = curInteractionStartIdx;
        maskEnd = curInteractionStartIdx+curInteractionFeatureCount-1;
    
        curInteractionMask.mask = zeros(featureCount,1);
        curInteractionMask.mask(maskStart:maskEnd) = 1;        
        curInteractionMask.name = curEventName;
    
        % add the current interaction mask to the list
        interactionMasks{end+1} = curInteractionMask;
    end
end


% TODO - verify that all features have at least 1 stimulus associated with
% it. If some don't, then remove them from the features matrix and
% featureNames vector, and write out a warning message

% write out the features matrix
save(featuresFilename, 'features', 'featureNames', 'regressorNames', 'regressorStartIndices', 'interactionMasks', 'valCompFeatures', 'excludeValFeatures');

end

function featureNames = createFeatureNames(regressorIdxs, options, parentName)

    % one the first call to this function no parentName is passed in, so
    % set it to the empty string
    if ~exist('parentName', 'var')
        parentName = '';
    end

    % get the current regressors name and index
    curIdx = regressorIdxs(1);
    curName = options.regressorNames{curIdx};
    
    % the next level down in the recurssion will not contain the current
    % regressor index
    nextRegressorIdxs = regressorIdxs(2:end);
    
    % define the cell array for this level of feature names
    featureNames = {};
    
    % if the current regressor is categorical then a name must be created
    % for each of it's levels
    if options.regressorCategorical(curIdx)
        
        % get the names of the regressor's levels
        categoryNames = options.regressorCategoricalNames{curIdx};
        categoryCount = numel(categoryNames);
        
        % iterate through all the levels of the categorical regressor and
        % add it to the list of feautre names, then 
        for curLevel = 1:categoryCount
            newFeatureName = [parentName curName '-' categoryNames{curLevel}];
            
            % if there are still more regressors, then keep recursing down
            if numel(nextRegressorIdxs) > 0
                
                % there are children, so this is an interaction term, thus put
                % an x indicating times in the name
                newParentName = [newFeatureName '_x_'];
                
                % call the function recursively
                childFeatures = createFeatureNames(nextRegressorIdxs, options, newParentName);
                featureNames = [featureNames childFeatures];
            % otherwise this is the bottom level, so add the new feature name
            % created here to the list to return
            else
                featureNames{end+1} = newFeatureName;
            end
        end
    else
        % create the new feature name by just concatonating the parent's
        % name with the current name
        newFeatureName = [parentName curName];
        
        % if there are still more regressors, then keep recursing down
        if numel(nextRegressorIdxs) > 0
            
            % there are children, so this is an interaction term, thus put
            % an x indicating times in the name
            newParentName = [newFeatureName '_x_'];
            
            % call the function recursively
            featureNames = createFeatureNames(nextRegressorIdxs, options, newParentName);
        % otherwise this is the bottom level, so add the new feature name
        % created here to the list to return
        else
            featureNames{end+1} = newFeatureName;
        end
    end
end