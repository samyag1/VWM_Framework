function model = mrifEstimate(Y,X,addBias, estType, crossvalFolds, runCount, crossvalType, averageCrossValFolds, useSingleLambda, XNu, excludeValFeatures,crossValSignificant)
%  model = mrifEstimate(Y,X,addBias);
%------------------------------------------------------------------------
% CROSSVALIDATED OLS
[nObs,nFeat] = size(X);

if notDefined('addBias')
	addBias = 0;
end
if notDefined('estType')
    estType = 'ridge';
end
if notDefined('crossvalFolds')
    crossvalFolds = 10;
end
if notDefined('runCount')
    runCount = 1;
end
if notDefined('crossvalType')
    crossvalType= 'Runs';
end
if notDefined('averageCrossValFolds')
    averageCrossValFolds = false;
end
if notDefined('useSingleLambda')
    useSingleLambda = true;
end
if notDefined('XNu')
    XNu = [];
end
if notDefined('excludeValFeatures')
    excludeValFeatures = [];
end
if notDefined('crossValSignificant')
    crossValSignificant = .05;
end

% if the run count isn't divisible by the CV count then don't
% partition the cross-validation folds by run
if mod(runCount, crossvalFolds) ~= 0
    error(sprintf('Trying to partition the cross-validation folds by runs, but the number of runs: %i, is not divisble by the crossvalFolds: %i provided',runCount, crossvalFolds));
end

% chunk the data up by runs
part = zeros(1,nObs);
obsPerRun = nObs / runCount;
if strcmp(crossvalType,'Runs')
    % can't do cross-validation over runs if there are more folds than runs
    if runCount < crossvalFolds
        error('Trying to do cross-validation over runs when there are are more runs: %i than cross-validation folds: %i.', runCount, crossvalFolds);
    end
    
    runParts = mod(1:runCount,crossvalFolds);
    runParts(runParts==0) = crossvalFolds;
    for i = 0:(runCount-1)
        part(i*obsPerRun+1:(i+1)*obsPerRun) = runParts(i+1);
    end
% or by volumes    
elseif strcmp(crossvalType,'Volumes')
    part = mod(1:nObs,crossvalFolds) + 1;
elseif isa(crossvalType, 'double') && (numel(crossvalType) == nObs)
    part = crossvalType;
else
    error('Invalid cross-validation type specified in mrifEstimate.');
end
    
% only do cross validation if told to
if strcmp(estType, 'ols')
    
    % TODO - add support for the nuissance regressors (bot XNu and
    % excludeValFeatures)
    
    if nObs ~= size(Y,2)
        Y = Y';
    end
    nVox = size(Y,1);
    
    weights = zeros(nFeat,nVox,crossvalFolds);
    pred  = zeros(size(Y));
    bias = zeros(crossvalFolds,nVox);
    if addBias
        weightsSE = zeros(nFeat+1,nVox,crossvalFolds);
    else
        weightsSE = zeros(nFeat,nVox,crossvalFolds);
    end
    
    for curFold = 1:crossvalFolds
        fprintf('\rResample %d out of %d',curFold,crossvalFolds);
        Xest = X(part~=curFold,:);
        Ytmp = Y(:,part~=curFold);
        if addBias
            [weights(:,:,curFold),bias(curFold,:),weightsSE(:,:,curFold)] = ols(Xest,Ytmp,[],1);
            pred(:,part == curFold) = ([ones(numel(find(part==curFold)),1),X(part ==curFold,:)]* ...
                [bias(curFold,:);weights(:,:,curFold)])';
        else
            [weights(:,:,curFold),~,weightsSE(:,:,curFold)] = ols(Xest,Ytmp,[],0);
            pred(:,part == curFold) = (X(part ==curFold,:)*weights(:,:,curFold))';
        end
    end
    
    % CALCULATE RESIDUALS
%    e = pred - Y;
    
    % GET THE CORRELATION COEFFICIENTS
    fprintf('\nCalculating correlations...');
    [cc,CI] = ccMatrix(pred,Y,2);
    fprintf('done\n');
    
    sigCC = pval2r(crossValSignificant,numel(find(part==1)));
    nGood = numel(find(cc > sigCC));
    
    model.weights = single(mean(weights,3));
    if addBias
        model.bias = single(mean(bias,1));
    end
    model.weightsSE = single(mean(weightsSE,3));
    model.cc = single(cc);
    model.ccCI = CI;
    model.sigCC = sigCC;
    model.nGood = nGood;
elseif strcmp(estType, 'ridge')

    % we want the Y matrix to be time in rows and voxels in columns, which
    % is the oposite for the above ols model
    if nObs ~= size(Y,1)
        Y = Y';
    end
    nVox = size(Y,2);
    
    % if there are nuissance regressors to add to the model, then
    % concatonate them with the design matrix
    if not(isempty(XNu))
        XFullModel = [X XNu];
        
        % now add the indices of these newly added nuissance regressors to
        % the list of features to exclude from cross-val and validation
        nuissanceRegressorIndices = [size(X,2)+1:size(X,2)+size(XNu,2)];
        excludeValFeatures = [excludeValFeatures nuissanceRegressorIndices];
    else
        XFullModel = X;
    end

    % add the beta0 terms to the design matrix if told to. 
    % MUST BE THE LAST COLUMNS IN DESIGN MATRIX, assumed so below
    if addBias
        % create additional columns (features), one for each run, to act as a
        % bias (beta0) term and concatonate it to the end of the current mat
        biasMat = zeros(nObs,runCount);
        for curRun = 1:runCount
            curRunStart = (curRun-1)*obsPerRun+1;
            curRunEnd = curRun*obsPerRun;
            biasMat(curRunStart:curRunEnd,curRun) = 1;
        end
        XFullModel = [XFullModel biasMat];
    end

    % create the indices that are to be used for prediction. THese will be used
    % for cross-validation of the estimation weights, as well as determine
    % which weights are stored.
    useValFeatures = [1:size(XFullModel,2)]';
    useValFeatures(excludeValFeatures) = [];
        
    % update the feature count to account for the bias terms and nuissance
    % regressors
    nFeat = numel(useValFeatures);
    nNuissanceFeat = numel(excludeValFeatures);
    
    % TODO - make this a parameter
    % use the range of lambdas that Dustin specified
	lambdas = logspace(-10,4,10);
    lambdaCount = numel(lambdas);
    
    % allocate the matrices for the weights and predication values
    weights = zeros(nFeat,nVox,lambdaCount,crossvalFolds);
    if nNuissanceFeat > 0
        nuissanceWeights = zeros(nNuissanceFeat,nVox,lambdaCount,crossvalFolds);
    end
    mseMaster = zeros(lambdaCount,crossvalFolds,nVox);
    ccMaster = zeros(lambdaCount,crossvalFolds,nVox);
%    cIMaster = zeros(lambdaCount,crossvalFolds,nVox,2);
    r2Master = zeros(lambdaCount,crossvalFolds,nVox);
    predMaster = zeros(nObs, nVox, lambdaCount); % the predicted values for each voxel at each time point for each lambda value

    % Do the cross-validation grid-search using a k-Fold
    % cross-validation over the values of alpha and lambda calculated above
    for curFold = 1:crossvalFolds
        fprintf('\rResample %d out of %d',curFold,crossvalFolds);
        
        % figure out the indices (voxels) to use for this CV iteration.
        % Exclude any volumes that have a 0 in the crossValAssignments
        estUseVector = part~=curFold;
        estUseVector(find(part==0)) = 0;
        valUseVector = part==curFold;
        
        % create the x and Y subsets to use in estimation and predication
        Xest = XFullModel(estUseVector,:);
        Yest = Y(estUseVector,:);
        Xval = XFullModel(valUseVector,useValFeatures);
        Yval = Y(valUseVector,:);
        
        % estimate using ridge regression with the lambdas specified. The
        % weights returned include those for the nuissance regressors and
        % bias terms
        curWeights = ridgemulti(Xest, Yest, lambdas);
        
        % store only the regressors of interest, throwing away the weights
        % estimated for the nuissance regressors, but including the bias
        % term weights
        weights(:,:,:,curFold) = curWeights(useValFeatures,:,:);

        if nNuissanceFeat > 0
            nuissanceWeights(:,:,:,curFold) = curWeights(excludeValFeatures,:,:);
        end
        
        % iterate through the lambda values estimated and calculate the
        % prediction the estimated models make for each, then calculate the
        % mean squared error and correlation coregigcations for the model
        for curLambda = 1:lambdaCount
            % calculate the predicted values for the current value of lambda
            pred = Xval*weights(:,:,curLambda,curFold);
            predMaster(valUseVector,:,curLambda) = pred;
            
            % if we're calculating prediction by averaging across folds,
            % then store those values now
            if averageCrossValFolds
                mseMaster(curLambda,curFold,:) = mean((pred - Yval).^2,1);
                [ccTemp, cITemp] = ccMatrix(pred, Yval,1);
                ccMaster(curLambda,curFold,:) = ccTemp;
                %            cIMaster(curLambda,:,:) = reshape(cITemp',nVox,2);
                r2Temp = calcR2(Yval,pred);
                r2Master(curLambda,curFold,:) = r2Temp;
            end
        end
    end
    fprintf('\nDetermining best cross-validation fit...');
    
    % determine the prediction accuracies across all folds
    ccFinal = zeros(nVox,lambdaCount);
%    cIFinal = zeros(nVox,lambdaCount,2);
    mseFinal = zeros(nVox,lambdaCount);
    r2Final = zeros(nVox,lambdaCount);
    
    % calculate the prediction accuracy
    for curVox = 1:nVox
        % if prediction is to be calculated across folds, then take the
        % average across folds
        if averageCrossValFolds
            % take the average mse and cc acorss cross-val folds for each value of
            % lambda
            mseFinal(curVox,:) = mean(mseMaster(:,:,curVox),2);
            ccFinal(curVox,:) = mean(ccMaster(:,:,curVox),2);
            r2Final(curVox,:) = mean(r2Master(:,:,curVox),2);
            
            % otherwise, we'll calculate a single correlation value per
            % lambda value
        else
            for curLambda = 1:lambdaCount
                mseFinal(curVox,curLambda) = mean((predMaster(:,curVox,curLambda) - Y(:,curVox)).^2,1);
                [ccTemp, cITemp] = ccMatrix(predMaster(:,curVox,curLambda), Y(:,curVox),1);
                ccFinal(curVox, curLambda) = ccTemp;
                %               cIMaster(curVox,curLambda,:) = reshape(cITemp',1,2);
                r2Final(curVox,curLambda) = calcR2(Y(:,curVox),predMaster(:,curVox,curLambda));
            end
        end
    end

    % calculate the r value that is significant at the value passed in
    if averageCrossValFolds
        sigCC = pval2r(crossValSignificant,numel(find(part==1)));
    else
        sigCC = pval2r(crossValSignificant,nObs);
    end        
    nGood = numel(find(ccFinal > sigCC));

    % store the prediction accuracy and other model output
    model.mse = single(mseFinal);
    model.cc = single(ccFinal);
%    model.ccCI = single(cIFinal);
    model.r2 = single(r2Final);
    model.sigCC = sigCC;
    model.nGood = nGood;
    model.lambdas = lambdas;
    
    % is we're to use a single lambda for all voxels, then simply write out
    % the full weights to file, so it can be determined in a later step
    % with all the chunks together
    if useSingleLambda
        
        model.weightsFull = weights;
        if nNuissanceFeat > 0
            model.nuissanceWeightsFull = nuissanceWeights;
        else
            model.nuissanceWeightsFull = [];
        end
    % otherwise determine which lambda works best for each voxel
    else
        weightsFinal = zeros(nFeat,nVox);
        if nNuissanceFeat > 0
            nuissanceWeightsFinal = zeros(nNuissanceFeat,nVox);
        end
       lambdasUsed = zeros(nVox,1);
       for curVox = 1:nVox
                        
            % get the lambda index that was best for the current voxel
            [~,bestLambdaIdx] = max(ccFinal(curVox,:));
            lambdasUsed(curVox) = bestLambdaIdx;
        
            % now get the set of sets of weights calculated for the max lambda
            % where there are K sets, one for each cross-validation iteration
            weightsFinal(:,curVox) = mean(squeeze(weights(:,curVox,lambdasUsed(curVox),:)),2);
            
            % if there are nuissance regressors, then store the weights from
            % the same lambdas that predicted best
            if nNuissanceFeat > 0
                nuissanceWeightsFinal(:,curVox) = mean(reshape(nuissanceWeights(:,curVox,lambdasUsed(curVox),:),nNuissanceFeat,crossvalFolds),2);
            end
        end
        
        % if a bias term was added then only save the actual feature weights
        % into the weights field, and store those bias terms seperately
        if addBias
            model.weights = single(weightsFinal(1:end-runCount,:));
            model.biasTerms = single(weightsFinal(end-runCount+1:end,:));
        else
            model.weights = single(weightsFinal);
            model.biasTerms = zeros(runCount,nVox);
        end
        if nNuissanceFeat > 0
            model.nuissanceWeights = single(nuissanceWeightsFinal);
        else
            model.nuissanceWeights = [];
        end
        
        model.lambdasUsed = lambdasUsed;
    end
    
    fprintf('done\n');
    
elseif strcmp(estType, 'elasticNet')
    % Question? don't we need to control for the end of one run into the
    % beginning of another so that the HRF is not convolved over it all? No
    % because we're not fitting an HRF and there are dummy scans at the end
    % of each run, so there will be stimuli being predicted from those last
    % volumes into the first volumes of the next run
     
    % TODO - add support for the nuissance regressors
   
    % we want the Y matrix to be time in rows and voxels in columns, which
    % is the oposite for the above ols model
    if nObs ~= size(Y,1)
        Y = Y';
    end
    nVox = size(Y,2);
    
    % alpha can take values from 0 to 1. use 100 over an equal grid of that
    % range. Start with something slightly higher than 0, as Dustin said
    % this version of Lasso produces strange results at values close to 0.
    alphas = 0.1:.1:1.0;
    alphaCount = numel(alphas);

    tic;
    % use all the data to determine the values of lambda to use
    lambdaCount = 10;
    lambdas = zeros(nVox,lambdaCount);
    for curVox = 1:nVox
        [~, dummyStats] = lasso(X,Y(:,curVox), 'Standardize', false, 'NumLambda', lambdaCount);
        lambdas(curVox,:) = dummyStats.Lambda;
    end
    toc;
    
    % allocate the matrices for the weights and predication values
    weights = zeros(nFeat,lambdaCount,alphaCount,nVox,crossvalFolds);
    %weightsSE = zeros(d,lambdaCount,alphaCount,nVox,crossvalFolds);
    pred  = zeros(alphaCount,lambdaCount,nObs,nVox);
    
    % Do the cross-validation grid-search using a k-Fold
    % cross-validation over the values of alpha and lambda calculated above
    for curFold = 1:crossvalFolds
        fprintf('\rResample %d out of %d',curFold,crossvalFolds);
        
        % figure out the indices (voxels) to use for this CV iteration
        estUseVector = part~=curFold;
        valUseVector = part==curFold;
        Xest = X(estUseVector,:);
        
        % iterate through all the voxels one at a time since elastic
        % net is an iterative algorithm, and so cannot estimate weights
        % for all voxels at one time like OLS
        for curVox = 1:nVox
            
            % get the features and BOLD value for the current voxel
            Yvox = Y(estUseVector,curVox);
            
            % iterate over the possible values of alpha
            for curAlphaIdx = 1:alphaCount
            
                % get the current alpha and lambdas to try
                curAlpha = alphas(curAlphaIdx);
                curLambdas = lambdas(curVox,:);
            
                % estimate using elasticNet with the give alpha
                [curWeights, ~] = lasso(Xest, Yvox, ...
                                               'Alpha', curAlpha, ...
                                               'Lambda', curLambdas, ...
                                               'Standardize', false);
                weights(:,:,curAlphaIdx,curVox,curFold) = curWeights;
                
                % TODO - it looks like SE for each predictor is not
                % returned by lasso. Maybe it can't be calcualted? Ask
                % Dustin about this.
                %weigthsSE 
                
                % calculate the predicted values for all lambda
                % values for the current alpha and voxel.
                pred(curAlphaIdx,:,valUseVector,curVox) = (X(valUseVector,:)*squeeze(weights(:,:,curAlphaIdx,curVox,curFold)))';
            end
        end
    end
    
    fprintf('\nCalculating correlations...');
    
    % calculate the residuals per voxel for all the values of alpha and
    % lambda, calculate the correlation coefficients for all, then choose
    % the values of alpha and lambda that best fits
    lambdasUsed = zeros(nVox,1);
    alphasUsed = zeros(nVox,1);
    ccMaster = zeros(nVox,1);
    cIMaster = zeros(nVox,2);
    cvWeights = zeros(nFeat,crossvalFolds,nVox);
    for curVox = 1:nVox
        ccHyperParams = zeros(alphaCount,lambdaCount);
        cIHyperParams = zeros(alphaCount,lambdaCount,2);
        
        for curAlpha = 1:alphaCount
            for curLambda = 1:lambdaCount
                
                % get the current prediction matrix for lambda and alpha
                curPred = squeeze(pred(curAlpha, curLambda,:,curVox));
                
                % CALCULATE RESIDUALS
                %            e = curPred - Y;
                
                % GET THE CORRELATION COEFFICIENT OF CURRENT CVOXEL
                [ccHyperParams(curAlpha,curLambda),cIHyperParams(curAlpha,curLambda,:)] = ccMatrix(curPred,Y(:,curVox),1);
            end
        end
        
        % find the max correlation correficient, and use those values of alpha
        % and lambda
        [maxIndicesAlpha, maxIndicesLambda] = find(ccHyperParams == max(ccHyperParams(:)));
        maxAlphaIdx  = maxIndicesAlpha(1);
        maxLambdaIdx = maxIndicesLambda(1);
        alphasUsed(curVox) = alphas(maxAlphaIdx);
        lambdasUsed(curVox) = lambdas(maxLambdaIdx);
        
        % now get the set of sets of weights calculated for the max lambda of alpha
        % where there are K sets, one for each cross-validation iteration
        cvWeights(:,:,curVox) = squeeze(weights(:,maxLambdaIdx,maxAlphaIdx,curVox,:));        

        % updated the correlation coeficients and confidence interval
        % vectors for the current voxel
        ccMaster(curVox) = ccHyperParams(maxAlphaIdx,maxLambdaIdx);
        cIMaster(curVox,:) = cIHyperParams(maxAlphaIdx,maxLambdaIdx,:);
    end
    
    sigCC = pval2r(crossValSignificant,numel(find(part==1)));
    nGood = numel(find(ccMaster > sigCC));
    
    model.weights = single(squeeze(mean(cvWeights,2)));
%    model.weightsSE = single(mean(weightsSE,3));
    model.cc = single(ccMaster);
    model.ccCI = cIMaster;
    model.sigCC = sigCC;
    model.nGood = nGood;
    model.alphas = alphasUsed;
    model.lambdas = lambdasUsed;
    
    fprintf('done\n');
else
    error('Invalid model estimation type provided through the mriFitOptions estType field: %s', estType);
end