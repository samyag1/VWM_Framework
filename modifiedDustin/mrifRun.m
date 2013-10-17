function opt= mrifRun(xpmt,chunkIdx,features,featuresNuissanceEst,featuresNuissanceVal,opt)
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
sessionDir = xpmt.sessionDir;

runFolds = xpmt.runFold{opt.condIdx};
nRuns = numel(runFolds);
opt.nRuns = nRuns;

% how to do the model estimation
estType = opt.estType;
crossValFolds = opt.crossValFolds;
modelSelectionIndices = opt.modelSelectionIndices;
averageCrossValFolds = opt.averageCrossValFolds;
crossValSignificant = opt.crossValSignificant;
useSingleLambda = opt.useSingleLambda;
doPermutationTesting = opt.doPermutationTesting;
permutationIters = opt.permutationIters;
hrfType = opt.hrfType;

% determine the FIR Bins variable based on the HRF estimation type
switch(hrfType)
    % if a deconvolution model is being used then we want to create a design
    % matrix with the features at the TR where the stimuli appeared, so set the
    % FIR Bin variable to 0, which will do just that
    case {'Deconvolve'}
        binsFIR = 0;
    
    % Otherwise use the fir bins passed in by the user
    case {'FIR'}
        
        % this is list of the bin offsets to use for the FIR model
        binsFIR = opt.binsFIR;
    otherwise
        error(sprintf('Invalid hrfType parameter provided: %s. Choices are \"Deconvolve\" and \"FIR"', hrfType));
end

onlyUseStimVols = opt.onlyUseStimVols;
onlyUseStimVolsOffset = opt.onlyUseStimVolsOffset;

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
            
            % create the design matrix for nuissance regressors with FIR Bins
            XNuissanceFIREst = [];
            if numel(featuresNuissanceEst) > 0
                [XNuissanceFIREst,~] = featuresToDesignMatrix(xpmt.stimFile{opt.condIdx}{run}, ...
                    featuresNuissanceEst, paradigm.presHz,paradigm.TR,paradigm.nDummy,0,[],binsFIR, stimAsFeatures);
            end
            
            if numel(opt.removeNuissanceFilenamesEst) > 0 || numel(XNuissanceFIREst) > 0
                Y0N = removeNuisssanceVariance(Y0, runFold, opt.removeNuissanceFilenamesEst, XNuissanceFIREst);
            else
                Y0N = Y0;
            end
        case {'val', 'snr'}
            
            % create the design matrix for nuissance regressors with FIR Bins
            XNuissanceFIRVal = [];
            if numel(featuresNuissanceVal) > 0
                [XNuissanceFIRVal,~] = featuresToDesignMatrix(xpmt.stimFile{opt.condIdx}{run}, ...
                    featuresNuissanceVal, paradigm.presHz,paradigm.TR,paradigm.nDummy,0,[],binsFIR, stimAsFeatures);
            end
            
            if numel(opt.removeNuissanceFilenamesVal) > 0 || numel(XNuissanceFIRVal) > 0
                Y0N = removeNuisssanceVariance(Y0, runFold, opt.removeNuissanceFilenamesVal, XNuissanceFIRVal);
            else
                Y0N = Y0;
            end
        otherwise
            Y0N = Y0;
    end
    clear Y0;
    
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
            windowSize = min(size(Y0N,1), windowSize); % ensure window is no bigger than time series
            dp = sgolayfilt(Y0N,opt.detrendPolyDeg,windowSize);
            Y0D = Y0N - dp;
            driftParams = [driftParams dp];
            N = [];
        otherwise
            Y0D = Y0N;
            N = [];
    end
    clear Y0N;
    
    % only zscore if the options say to
    if opt.zscore
        % z-score the data
        Y0Z = zScore(Y0D);    
    else
        Y0Z = Y0D;    
    end
    
	%--------------------------------------------------------------------
	% CONCATENATE STIMULUS FEATURE MATRICES
	[X0,curStimOrder] = featuresToDesignMatrix(xpmt.stimFile{opt.condIdx}{run}, ...
								features, paradigm.presHz,paradigm.TR,paradigm.nDummy,0,[],binsFIR, stimAsFeatures);
                            
    if numel(opt.modelNuissanceFilenames) > 0
        XNu0 = readNuissanceFiles(runFold, opt.modelNuissanceFilenames);
    else
        XNu0 = [];
    end

    % remove all trials except those witha given offset from stim
    % presentation to account for variance due to image on/off
    if onlyUseStimVols
        useVols = find(curStimOrder~=0)+onlyUseStimVolsOffset;
        Y0D = Y0D(useVols,:);
        Y0Z =Y0Z(useVols,:);
        X0 = X0(useVols,:);
        curStimOrder = curStimOrder(useVols,:);
        if numel(XNu0) > 0
            XNu0 = XNu0(useVols,:);
        end
    end
                            
    % only zscore if the options say to
    if opt.zscore
        X0Z = zscore(X0);

        if numel(opt.modelNuissanceFilenames) > 0
            XNu0Z = zscore(XNu0);
        else
            XNu0Z = XNu0;        
        end
    else
        X0Z = X0;
        XNu0Z = XNu0;
    end

    % Sam - TODO - if the X is stimAsFeatures, then concatonating doesn't
    % work, I need to build an empty matrix the size of all stim in all
    % runs and fill in the diagonal blocks with the matrices created above
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
    
    Y = [Y; Y0D];
    YZ = [YZ; Y0Z];

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

            % write out the design matrices
            if opt.writeDesignMatrix
                designMatFile = fullfile(sessionDir,'designMatrixEst.mat');
                if not(exist(designMatFile, 'file'))
                    save(designMatFile, 'X', 'XZ', 'XNu', 'XNuZ');
                end
            end
            
            fprintf('Calculating feature weights.\n');

            switch(hrfType)            
                case {'FIR'}
    
                    % hold out 10% of the data for use in model/voxel selection
                    % if modelSelectionIndices are given, which represent
                    % the runs to hold out
                    if ~isempty(modelSelectionIndices)
                        nObs = size(YZ,1);
                        obsPerRun = nObs / runCount;
                        selectionIndices = zeros(nObs,1);
                        for i = modelSelectionIndices
                            startInd = (i-1)*obsPerRun+1;
                            endInd = i*obsPerRun;
                            selectionIndices(startInd:endInd) = 1;
                        end
                        modelSelectionY = YZ(selectionIndices,:);
                        modelSelectionX = XZ(selectionIndices,:);
                        modelEstimationY = YZ(~selectionIndices,:);
                        modelEstimationX = XZ(~selectionIndices,:);
                    else
                        modelSelectionY = [];
                        modelSelectionX = [];
                        modelEstimationY = YZ;
                        modelEstimationX = XZ;
                    end
                    
                    % estimate the model with the FIR binned design matrix
                    model=mrifEstimate(modelEstimationY,modelEstimationX,addBias, estType, crossValFolds, nRuns, 'Runs', averageCrossValFolds, useSingleLambda, XNuZ, excludeValFeatures, crossValSignificant);
                    
                    % if there is model selection data, then predict it
                    if ~isempty(modelSelectionY)
                        % calculate the prediction accuracy
                        modelSelection = mrifForward(modelSelectionY,modelSelectionX,model.weights,addBias,[]);
                        model.heldOutCC = modelSelection.cc;
                        model.heldOutCI = modelSelection.cI;
                        model.heldOutPred = modelSelection.pred;
                        model.heldOutResid = modelSelection.resid;
                        
                        % if permutation testing is requested, do it here
                        % and store the results on the model to be written
                        % to disk
                        if doPermutationTesting
                            [nullDist, pVals] = permutationTest(modelSelectionY,respAmpSelectionX,model.weights,addBias,[],permutationIters, model.heldOutCC); 
                            model.heldOutPVals = pVals;
                            model.heldOutNullDist = nullDist;
                        end
                    end
                    
                case {'Deconvolve'}
                    
                    % if the options says to load the hrf values from file,
                    % then do so
                    if opt.loadHRFFromFile
                        hrfFile = fullfile(opt.hrfFolder,sprintf('%04d.mat',opt.chunkNum));
                        hrf =loadSubField(hrfFile,'model','hrf');
                        hrfBasisWeights =loadSubField(hrfFile,'model','hrfBasisWeights');
                        hrfBasis =loadSubField(hrfFile,'model','hrfBasis');
                    % otherwise estimate the hrfs now
                    else
                        
                        % determine the hrf function to use in the
                        % deconvolution
                        hrfFit = estimateHRF(YZ, X, paradigm.TR);
                        hrf = hrfFit.hrf;
                        hrfBasisWeights = hrfFit.basisParams;
                        hrfBasis = hrfFit.basis;
                    end
                    
                    % deconvolve the hrf from the raw data creating
                    % response amplitudes, the response of the voxels to
                    % each stimuli
                    [respAmp, respStimOrder] = deconvolveHRF(YZ, stimOrder, hrf);

                    % hold out 10% of the data for use in model/voxel selection
                    % For now, require model selection indices. If none are
                    % provided that means no data is held out
                    if ~isempty(modelSelectionIndices)
                        respSelectionStim = arrayfun(@(x)find(respStimOrder==x,1),modelSelectionIndices);
                        respEstimationStim = respStimOrder;
                        respEstimationStim(respSelectionStim) = [];
                        modelSelectionData = respAmp(respSelectionStim);
                        modelEstimationData = respAmp;
                        modelEstimationData(respSelectionStim) = [];
                    else
                        respSelectionStim = [];
                        respEstimationStim = respStimOrder;
                        modelSelectionData = [];
                        modelEstimationData = respAmp;
                    end
                                                            
                    % use only the features of the stimuli in the response
                    % amplitudes to create the design matrix onto which
                    % we'll regress the response amplitudes
                    % ASSUME that the first repitition's set of features is
                    % correct. Can't do per repitition features in this
                    % model
                    respAmpSelectionFeatures = features(respSelectionStim,:,1);
                    respAmpEstimationFeatures = features(respEstimationStim,:,1);
                    
                    % estimate the features model on the response
                    % amplitudes, and not the raw BOLD time series
                    model=mrifEstimate(modelEstimationData,respAmpEstimationFeatures,0, estType, crossValFolds, nRuns, 'Runs', averageCrossValFolds, useSingleLambda, [], [], crossValSignificant);
                    
                    % if there is held out data, then predict it
                    if ~isempty(modelSelectionData)
                        % calculate the prediction accuracy
                        modelSelection = mrifForward(modelSelectionData,respAmpSelectionFeatures,model.weights,addBias,[]);
                        model.heldOutCC = modelSelection.cc;
                        model.heldOutCI = modelSelection.cI;
                        model.heldOutPred = modelSelection.pred;
                        model.heldOutResid = modelSelection.resid;
                        model.heldOutIndices = respSelectionStim;
                        
                        % if permutation testing is requested, do it here
                        % and store the results on the model to be written
                        % to disk
                        if doPermutationTesting
                            [nullDist, pVals] = permutationTest(modelSelectionData,respAmpSelectionFeatures,model.weights,addBias,[],permutationIters, model.heldOutCC); 
                            model.heldOutPVals = pVals;
                            model.heldOutNullDist = nullDist;
                        end
                    end
                    
                    % store the HRF info
                    model.hrf = hrf;
                    model.hrfBasisWeights = hrfBasisWeights;
                    model.hrfBasis = hrfBasis;
                    model.respAmp = respAmp;
                    model.respStimOrder = respStimOrder;
            end
            
            % now predict the held out chunk and store the accuracy
            
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
            
            % write out the design matrices
            if opt.writeDesignMatrix
                designMatFile = fullfile(sessionDir,'designMatrixVal.mat');
                if not(exist(designMatFile, 'file'))
                    save(designMatFile, 'X', 'XZ', 'XNu', 'XNuZ');
                end
            end
            
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
                switch(hrfType)
                    case {'FIR'}
                        
                        % calculate the prediction accuracy
%                        curModel = mrifForward(YZ,XZ,weights,addBias, excludeValFeatures);
 
                        % store the variables to user in prediction
                        predY = YZ;
                        predX = XZ;
                        predExcludeValFeatures = excludeValFeatures;
                        
                    case {'Deconvolve'}
                        
                        % load the hrf's estimated during the estimation
                        % step, one for each voxel
                        voxHRFs = loadSubField(estFile,'model','hrf');
                        
                        % deconvolve the hrf from the raw data creating
                        % response amplitudes, the response of the voxels to
                        % each stimuli
                        [respAmp, respStimOrder] = deconvolveHRF(YZ, stimOrder, voxHRFs);
                        
                        % use only the features of the stimuli in the response
                        % amplitudes to create the design matrix onto which
                        % we'll regress the response amplitudes
                        respAmpFeatures = features(respStimOrder,:,1);
                        
                        % store the variables to use in prediction
                        predY = resAmp;
                        predX = respAmpFeatures;
                        predExcludeValFeatures = [];
                        
                        % calculate the prediction accuracy
                        %curModel = mrifForward(respAmp,respAmpFeatures,weights,addBias,[]);
                end
                                
                % calculate the prediction accuracy
                model = mrifForward(predY,predX,weights,addBias,predExcludeValFeatures);
                
                % if permutation testing is requested, do it here
                % and store the results on the model to be written
                % to disk
                if doPermutationTesting
                    [nullDist, pVals] = permutationTest(predY,predX,weights,addBias,predExcludeValFeatures,permutationIters,model.cc);
                    model.pVals = pVals;
                    model.nullDist = nullDist;
                end
            end
            
            % create the feature x FIR bin prediction matrix
            model.averageValRuns = opt.averageValRuns;
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
                XStim = zeros(dataCount, totalStimIDFeatures);
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
                                            
                    XStim(curRunStartTimeIdx:curRunEndTimeIdx,curRunStartFeatureIdx:curRunEndFeatureIdx) = curRunX;
                end
                
                if opt.zscore
                    XStimZ = zscore(XStim);
                else
                    XStimZ = XStim;
                end
                
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
                model=mrifEstimate(YZ,XStimZ,addBias, estType, repetitionCount, [], nRuns, crossValAssignments, 'None', false, useSingleLambda);
            
                % write out the design matrices
                if opt.writeDesignMatrix
                    designMatFile = fullfile(sessionDir,'designMatrixSNR.mat');
                    if not(exist(designMatFile, 'file'))
                        save(designMatFile, 'XStim', 'XStimZ');
                    end
                end
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
            if opt.detrendPolyDeg > 0
                model.driftBasis = N;
                model.driftParams = driftParams;
            end
                            
            fprintf('Saving snr data to:\n%s\n',snrFile);
            save(snrFile, 'model');
        end
end
     
