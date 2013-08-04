function opt= mrifRun(xpmt,chunkIdx,features,opt)
%opt= mrifRun(xpmt,chunkIdx,features,opt)
%------------------------------------------------------------------

if notDefined('xpmt')
     error(sprintf('Need to provide an experimental data \nstructure for a single session!!\n'));
end

if notDefined('opt')
	opt = mrifOpts;
end

condIdx = opt.condIdx;

if ~isfield(opt,'dataPre')
     dataPrefix = 'RCMvolume';
else
     dataPrefix = opt.dataPre;
end

addBias = opt.addBias;


estDir = xpmt.estDir;
valDir = xpmt.valDir;
snrDir = xpmt.snrDir;

runFolds = xpmt.runFold{opt.condIdx};
nRuns = numel(runFolds);
opt.nRuns = nRuns;

% how to do the model estimation
estType = opt.estType;
crossValFolds = opt.crossValFolds;
averageCrossValFolds = opt.averageCrossValFolds;
crossValSignificant = opt.crossValSignificant;

% this is list of the bin offsets to use for the FIR model
binsFIR = opt.binsFIR;

% ***NOTE*** THIS ASSUMES THE BINS ARE REPEATED IN BLOCKS. IF THAT CHANGES
% (TO REPEAT EACH REGRESSOR TOGETHER FOR EXAMPLE) THIS MUST CHANGE TOO!!
% These are the features to be excluded from the design matrix when doing
% validation because they're nuissance regressors to be put in FIR bins. 
% update the excludeValFeatures to reflect the repeated structure of the
% FIR bins
excludeValFeatures = [];
for curExcludeFeature = opt.excludeValFeatures
    for binIdx = 1:numel(binsFIR)
        newFeature = (binIdx-1)*size(features,2) + curExcludeFeature;
        excludeValFeatures = [excludeValFeatures newFeature];
    end
end

% DISPLAY INITIAL INFO
fprintf('\nNOTE: using the following options:\n');
disp(opt)

% SAM - changed this to load the first paradigm in the "stimFiles" field,
% as that allows a paradigm per run, whereas this 'paradigmFile' field is
% an old one that was used when all runs had the same stim order
% LOAD PARADIGM INFO
%tmp = load(xpmt.paradigmFile{opt.condIdx});
load(xpmt.stimFile{opt.condIdx}{1}, 'paradigm');

%-------------------------------------------------------------------
% LOAD DATA ACROSS RUNS
Y = [];
YZ = [];
X = [];
XZ = [];
XNu = [];
XNuZ = [];
stimOrder = [];
driftParams = [];
fprintf('\nLoading Data (%d runs):\n',nRuns);

% LOOP OVER RUNS
stimAsFeatures = opt.mode == 'snr';
XRunsList = {};
XRunsListZ = {};
XRunsNextSame = zeros(nRuns,1);
prevSortedStimOrder = [];
totalStimIDFeatures = 0;
totalTimeSize = 0;
for run = 1:nRuns
    runFold =runFolds{run}; % ASSUMES runFolds IS DEFINED!
    fprintf('\r---Run %d--->:%s',run,runFold);
    
    %-----------------------------------------------------------------
    % DATA LOADING STUFF ...
    Y0 = loadVoxDat(fullfile(runFold,[dataPrefix,'*.nii']),chunkIdx);	
	
    % remove variance do to nuissance regressors from the BOLD signal
    switch lower(opt.mode)
        case 'est'
            if numel(opt.removeNuissanceFilenamesEst) > 0
                Y0N = removeNuisssanceVariance(Y0, runFold, opt.removeNuissanceFilenamesEst);
            else
                Y0N = Y0;
            end
        case {'val', 'snr'}
            if numel(opt.removeNuissanceFilenamesVal) > 0
                Y0N = removeNuisssanceVariance(Y0, runFold, opt.removeNuissanceFilenamesVal);
            else
                Y0N = Y0;
            end
    end
    
    % only do polynomial detrending if the options say to
    switch opt.detrendType
        case 'poly'
            % POLYNOMIAL DETRENDING
            [Y0D,dp,N] = detrendPoly(Y0N,opt.detrendPolyDeg);
            driftParams = [driftParams,dp];
        case {'SG', 'SavitzkyGolay'}
            % SAVITZKY-GOLAY DETRENDING
            nMinutes = opt.detrendWindow;
            windowSize = round(nMinutes*60/xpmt.TR{opt.condIdx});
            windowSize = windowSize + mod(mod(windowSize,2)+1,2); % ENSURE ODD WINDOW
            windowSize = min(size(Y0,1), windowSize); % ensure window is no bigger than time series
            dp = sgolayfilt(Y0,opt.detrendPolyDeg,windowSize);
            Y0D = Y0 - dp;
            driftParams = [driftParams dp];
            N = [];
        otherwise
            Y0D = Y0N;
            N = [];
    end
    
    % only zscore if the options say to
    if opt.zscore
        % z-score the data
        Y0Z = zScore(Y0D);    
    else
        Y0Z = Y0D;    
    end
    
    Y = [Y; Y0];
    YZ = [YZ; Y0Z];
    clear Y0;    
    
	%--------------------------------------------------------------------
	% CONCATENATE STIMULUS FEATURE MATRICES
	[X0,curStimOrder] = featuresToDesignMatrix(xpmt.stimFile{opt.condIdx}{run}, ...
								features, paradigm.presHz,paradigm.TR,paradigm.nDummy,0,[],binsFIR, stimAsFeatures);

    if numel(opt.modelNuissanceFilenames) > 0
        XNu0 = readNuissanceFiles(runFold, opt.modelNuissanceFilenames);
        if opt.zscore
            XNu0Z = zscore(XNu0);
        else
            XNu0Z = XNu0;
        end
    else
        XNu0 = [];
        XNu0Z = [];
    end
                            
    % only zscore if the options say to
    if opt.zscore
        X0Z = zscore(X0);
    else
        X0Z = X0;
    end

    % Sam - TODO - if the X is stimAsFeatures, then concatonating doesn't
    % work, I need to build an empty matrix the size of all stim in all
    % runs and fill in the "diagonals" with the matrices created above
    if stimAsFeatures
        XRunsListZ{end+1} = X0Z;
        XRunsList{end+1} = X0;
        curSortedStimOrder = sort(unique(curStimOrder));
        curSortedStimOrder(curSortedStimOrder == 0) = [];
        
        if not(isempty(prevSortedStimOrder)) && sum(curSortedStimOrder == prevSortedStimOrder) == numel(curSortedStimOrder)
            XRunsNextSame(run-1) = true;
        else
            totalStimIDFeatures = totalStimIDFeatures + size(X0Z,2);
        end            
        prevSortedStimOrder = curSortedStimOrder;
    else
        X = [X;X0];
        XZ = [XZ;X0Z];
        XNu = [XNu;XNu0];
        XNuZ = [XNuZ;XNu0Z];
    end

    stimOrder = [stimOrder;curStimOrder];

end % (END RUN LOOP)s

fprintf('\nDone.\n');
fprintf('\nProcessing:');
switch lower(opt.mode)
    % - HRF/RESPONSE AMPLITUDE ESTIMATION
    case {'est'}
		estFile = fullfile(estDir,sprintf('%04d.mat',opt.chunkNum));
        
        if exist(estFile, 'file')
            fprintf('Estimation Chunk: %i already estimated. Skipping estimation', opt.chunkNum);
        else
            fprintf('Calculating feature weights.\n');
            model=mrifEstimate(YZ,XZ,addBias, estType, crossValFolds, nRuns, 'Runs', averageCrossValFolds, XNuZ, excludeValFeatures, crossValSignificant);
            model.voxFit = chunkIdx;
            model.driftBasis = N;
            model.driftParams = driftParams;
            model.info = paradigm;
            
            fprintf('\nSaving estimated parameters to:\n%s\n',estFile);
            save(estFile,'model');
        end
    % VALIDATION RESPONSES ESTIMATION
    case {'val'}       
		valFile = fullfile(valDir,sprintf('%04d.mat',opt.chunkNum));

        if exist(valFile, 'file')
            fprintf('Validation Chunk: %i already validated. Skipping validation', opt.chunkNum);
        else
            estFile = fullfile(estDir,sprintf('%04d.mat',opt.chunkNum));
            weights =loadSubField(estFile,'model','weights');
            
            % take the mean of the each validation run and calculate
            % prediction accuracy on those values. This assumes that the
            % validation images had random images previous to them so the
            % effect of the previous stimulus averages out
            model = [];
            if opt.averageValRuns
                
                % get the number of features being modeling, disregarding
                % FIR binning
                featureCount = size(features, 2);
                
                % iterate over the FIR bins used, averaging each one before
                % estimating predictions
                for curFIRBinOffsetIdx = 1:numel(opt.binsFIR)
                    
                    % calculate the mean value per stimuli for the current
                    % Bin offset. These will be returned in the stimOrder
                    % order.
                    curFIRBinOffset = opt.binsFIR(curFIRBinOffsetIdx);
                    [meanVar, uniqueStimOrder] = calculateMeanVariance(YZ, stimOrder, curFIRBinOffset);
                    YMeanZ = meanVar.meanAll;
    
                    % TODO - BROKEN - after I aded the capacity to have
                    % different feature values per repetition of a stim
                    % this no longer works, or at least, it assumes that
                    % the first value is correct.
                    % get the features for the val stim
                    XMean = features(uniqueStimOrder,:,1);
                    XMeanZ = zscore(XMean);
                    
                    % get the weights from the current bin offset
                    curBinStart = (curFIRBinOffsetIdx-1) * featureCount + 1;
                    curBinEnd = curFIRBinOffsetIdx*featureCount;
                    curBinWeights = weights(curBinStart:curBinEnd,:);
  
                    if addBias
                        try
                            bias = loadSubField(estFile,'model','biasTerms');
                            curBinWeights = [mean(bias,1);curBinWeights];
                        catch
                            addBias = 0;
                        end
                    end
                    
                    % calculate prediction accuracy for the current FIR Bin
                    curModel = mrifForward(YMeanZ, XMeanZ, curBinWeights, addBias, excludeValFeatures);

                    % store the current FIR bin's validation values into the model
                    % to save to file
                    ccField = sprintf('cc_bin%i',curFIRBinOffset);
                    cIField = sprintf('cI_bin%i',curFIRBinOffset);
                    predField = sprintf('pred_bin%i',curFIRBinOffset);
                    residField = sprintf('resid_bin%i',curFIRBinOffset);
                    model.(ccField) = curModel.cc;
                    model.(cIField) = curModel.cI;
                    model.(predField) = curModel.pred;
                    model.(residField) = curModel.resid;                    
                    
                    % add the correlation coefficients by feature level for 
                    % the features specified in the options
                    for curCompFeature = opt.valCompFeatures
                        compFeatureTimeIdxs = find(XMean(:,curCompFeature)==1);
                        [cc, cI] = ccMatrix(curModel.pred(compFeatureTimeIdxs,:),YMeanZ(compFeatureTimeIdxs,:),1);
                        ccFieldFeat = sprintf('cc_bin%i_valCompFeature%04i', curFIRBinOffset, curCompFeature);
                        model.(ccFieldFeat) = cc;
                    end
                end                
            else
  
                if addBias
                    try
                        bias = loadSubField(estFile,'model','biasTerms');
                        weights = [mean(bias,1);weights];
                    catch
                        addBias = 0;
                    end
                end
                
                fprintf('Calculating validation responses.\n')
                curModel = mrifForward(YZ,XZ,weights,addBias, excludeValFeatures);
                
                model.cc = curModel.cc;
                model.cI = curModel.cI;
                model.pred = curModel.pred;
                model.resid = curModel.resid;
            end
            
            % create the feature x FIR bin prediction matrix
            model.averageValRuns = opt.averageValRuns
            model.binsFIR = opt.binsFIR;
            model.voxFit = chunkIdx;
            model.driftBasis = N;
            model.driftParams = driftParams;
            model.info = paradigm;
            
            fprintf('Saving validation models to:\n   %s\n',valFile);
            save(valFile,'model');
        end
    case {'snr'}
		snrFile = fullfile(snrDir,sprintf('%04d.mat',opt.chunkNum));

        if exist(snrFile, 'file')
            fprintf('Signal to Noise Chunk: %i already calculated. Skipping SNR', opt.chunkNum);
        else
            fprintf('Calculating signal to noise.\n');

            % determine the unique stim IDs used here
            uniqueStim = unique(stimOrder);
            uniqueStim = uniqueStim(uniqueStim ~= 0);
            uniqueStim = sort(uniqueStim);
            uniqueStimCount = numel(uniqueStim);
                
            % determine the number of voxels
            voxCount = size(YZ,2);
            dataCount = size(YZ,1);
            
            % determine the repetition count OF each stimuli
            repetitionCount = numel(find(stimOrder==uniqueStim(1)));
            
            % instead calculate use Mark's method for calculating R2 on a
            % linear regression fit of stimuli as regressors
            if opt.averageValRuns
                % TODO - Sam - debug the below code, add a calculation of
                % the F-stat, then update the writeOutput function to make
                % niftis of all these

                % iterate over the FIR bins used, averaging each one before
                % estimating predictions
                model = [];
                for curFIRBinOffsetIdx = 1:numel(binsFIR)
                    
                    % calculate the mean value per stimuli for the current
                    % Bin offset. These will be returned in the stimOrder
                    % order.
                    curFIRBinOffset = binsFIR(curFIRBinOffsetIdx);
                    
                    % allocate the matrix to store the repetion values for each stim
                    outputMatrix = NaN(uniqueStimCount, repetitionCount, voxCount);
                    
                    % iterate over the unique stim values and add each one to the final
                    % matrix
                    for curRunStimIdx = 1:uniqueStimCount
                        
                        % get the current stimID
                        curStim = uniqueStim(curRunStimIdx);
                        
                        % store the values of all repetitions where the stim was presented
                        % plus the current FIR Bin offset
                        stimTimeIndices = find(stimOrder == curStim);
                        stimTimeIndices = stimTimeIndices + curFIRBinOffset;
                        stimTimeIndices = stimTimeIndices(stimTimeIndices <= dataCount);
                        curStimValues = YZ(stimTimeIndices,:);
% SAM - since I'm using an offset I don't have data points for all stim
% (why not, I should?) so figure out how to handle
                        outputMatrix(curRunStimIdx, 1:size(curStimValues,1),1:size(curStimValues,2)) = curStimValues; 
                    end
                
                    % calculate explainable variance assuming the mean of the
                    % repeats for each image is the "true" signal
                    fieldR2 = sprintf('r2_bin%i', curFIRBinOffset);
                    model.(fieldR2) = markR2(outputMatrix);
                end
            % calculate predictions based on average correlation 
            % of a model fit
            else
                % since we're making the stim actual features, make a design
                % matrix that is as wide as the stim count for all the runs and
                % set the design matrices above along the diagonals
                XZ = zeros(dataCount, totalStimIDFeatures);
                timeCounter = 0;
                featureCounter = 0;
                for curRunXIdx = 1:numel(XRunsList)
                    
                    curRunX = XRunsList{curRunXIdx};

                    curRunStartTimeIdx = timeCounter + 1;
                    timeCounter = timeCounter + size(curRunX,1);
                    curRunEndTimeIdx = timeCounter;

                    curRunStartFeatureIdx = featureCounter + 1;
                    curRunEndFeatureIdx = featureCounter + size(curRunX,2);
                    if not(XRunsNextSame(curRunXIdx))
                        featureCounter = curRunEndFeatureIdx;
                    end
                                            
                    XZ(curRunStartTimeIdx:curRunEndTimeIdx,curRunStartFeatureIdx:curRunEndFeatureIdx) = curRunX;
                end
                
              %  XZ = zscore(XZ);
                
                volCount = size(YZ,1);
            
                % create a list indicating which cross-validation fold each
                % volume will be included in. Use 2 & 3 TRs past the
                % stimuli presentation of each repetition
                repetitionCounts = ones(1,uniqueStimCount);
                crossValAssignments = zeros(1,volCount);
                prevAssignment = 0;
                offset = 2;
                for curVolume = 1:(volCount-offset)

                    % if the current stim plus it's offset is 0, that means
                    % no stim was presented at this time point, so use the
                    % previouosly assigned fold
                    if stimOrder(curVolume) == 0
                        crossValAssignments(curVolume+offset) = prevAssignment;
                    % otherwise this is a stim
                    else
                        % figure out the index of this unique stim
                        stimIdx = find(uniqueStim == stimOrder(curVolume));
                        
                        % only include repetitions that are not excluded
                        curCrossVal = repetitionCounts(stimIdx);
                        if not(ismember(opt.noiseCeilingIgnoreReps, curCrossVal))
                            % set the cross val assignment as the repetition of
                            % this stimuli
                            crossValAssignments(curVolume+offset) = curCrossVal;
                            
                            % store this assignment for use in the next null TR
                            prevAssignment = curCrossVal;
                        else
                            % otherwise just make this trial (and the
                            % following) 0 so it is ignored
                            crossValAssignments(curVolume+offset) = 0;
                            prevAssignment = 0;
                        end
                        
                        % update the repetition count for this unique stim
                        repetitionCounts(stimIdx) = repetitionCounts(stimIdx) + 1;
                    end
                end
                
                % This is estimating weights for stimuli as features, which
                % should give an estimate of explainable variance
                model=mrifEstimate(YZ,XZ,addBias, estType, repetitionCount, nRuns, crossValAssignments);
            end
            
            if opt.writeMeanVar

                % iterate over the FIR bins used, averaging each one before
                % estimating predictions
                for curFIRBinOffsetIdx = 1:numel(opt.binsFIR)
                    
                    % calculate the mean value per stimuli for the current
                    % Bin offset. These will be returned in the stimOrder
                    % order.
                    curFIRBinOffset = opt.binsFIR(curFIRBinOffsetIdx);
                    [curMeanVar, uniqueStimOrder] = calculateMeanVariance(Y, stimOrder, curFIRBinOffset);
                    
                    % store the current FIR bin's validation values into the model
                    % to save to file
                    meanField = sprintf('mean_bin%i',curFIRBinOffset);
                    varField = sprintf('var_bin%i',curFIRBinOffset);
                    model.(meanField) = curMeanVar.meanAll;
                    model.(varField) = curMeanVar.varAll;
                end
                
                % only one ordering of the stim occurs, so save that here
                model.stimOrder = uniqueStimOrder;
            end
            
            % store the general values here
            model.voxFit = chunkIdx;
            model.info = paradigm;
            if opt.polyDeg > 0
                model.driftBasis = N;
                model.driftParams = driftParams;
            end
                            
            fprintf('Saving snr data to:\n%s\n',snrFile);
            save(snrFile, 'model');
        end
end
     