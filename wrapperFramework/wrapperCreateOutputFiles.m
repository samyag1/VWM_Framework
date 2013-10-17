function wrapperCreateOutputFiles(experiment, subj, sessionIdxs, featureNames, interactionMasks, writeTopVoxels, mrifOptions)

TOP_VOXEL_COUNT = 20;

if notDefined('sessionIdxs')
     sessionIdxs = 1:numel(exper.(subj));
end

% TODO
% - Possibly change the design matrix so that the FIR bins of each feature
%   are next to each other, instead of the entire feature set being repeated

% get the FIR bins that were used for this study. These are 0 based
% determine the FIR Bins variable based on the HRF estimation type
switch(mrifOptions.hrfType)
    % if a deconvolution model is being used then we want to create a design
    % matrix with the features at the TR where the stimuli appeared, so set the
    % FIR Bin variable to 0, which will do just that
    case {'Deconvolve'}
        binsFIR = 0;
    
    % Otherwise use the fir bins passed in by the user
    case {'FIR'}
        
        % this is list of the bin offsets to use for the FIR model
        binsFIR = mrifOptions.binsFIR;
    otherwise
        error(sprintf('Invalid hrfType parameter provided: %s. Choices are \"Deconvolve\" and \"FIR"', hrfType));
end

binCount = numel(binsFIR);
featureCount = numel(featureNames);

% read in the header of the first EPI from the first run to use it's header information
% for the niftis to be written here
firstRunFolder = experiment.(subj){1}.runFold{1}{1};
files = dir(fullfile(firstRunFolder,[mrifOptions.dataPre,'*']));
firstFilename = fullfile(firstRunFolder, files(3).name);
header = spm_vol(firstFilename);

% IMPORTANT set the header's datatype value to be a 32bit float, instead of
% a 16 bit int. BOLD files are recorded in 16-bit ints, but the files being 
% saved here are either beta maps or correlation maps, both of which are floats.
header.dt(1,1) = 16;

% get the voxel dimensions for the niftis to write
volumeVoxCount = prod(header.dim);

% iterate through all the sessions to write
for session = sessionIdxs
    
    xpmt = experiment.(subj){session};
    
    % retrieve the folders where the estimation and validation models are
    % written
    estDir = xpmt.estDir;
    valDir = xpmt.valDir;
    snrDir = xpmt.snrDir;
    niftisDir = xpmt.niftisDir;
    
    % get all the estimated beta weights
    betas = concatSubFields('model','weights',2,estDir);

    % load the indices into the original nifti that were modeled
    voxIdxs = concatSubFields('model','voxFit',1,estDir);
    if numel(voxIdxs) == 0
        voxIdxs = 1:volumeVoxCount;
    end
    
    % allocate a matrix to store the betas calculated using the interaction
    % masks
    interactionMaskCount = numel(interactionMasks);
    voxCount = size(betas,2);
    eventBetas = zeros(interactionMaskCount*binCount,voxCount);
    
    % iterate through the FIR bins and make a map for each bin value of
    % each feature
    for curBin = 1:binCount
        
        % get the betas for the current bin since they're stored in blocks
        % that are featureCount long
        betasStart = (curBin-1)*featureCount + 1;
        betasEnd = curBin*featureCount;
        curBinBetas = betas(betasStart:betasEnd,:);
        
        curEventBetas = zeros(interactionMaskCount, voxCount);
        
        % iterate through the features in the beta weights and make a map for
        % each one
        for curFeatureIdx = 1:numel(featureNames)
            
            % ge the name of the current feature
            curFeatureName = featureNames{curFeatureIdx};
            
            % allocate a matrix to store the betas for the entire volume
            curBetaMap = zeros(1,volumeVoxCount);
            curBetaMap(voxIdxs) = curBinBetas(curFeatureIdx,:);
            
            % write out the current feature's betas
            curFeatureFilename = sprintf('betas-%s-Bin%i.nii', curFeatureName, curBin);
            curFeatureFilename = fullfile(niftisDir, curFeatureFilename);
            writeVol(curFeatureFilename, header, curBetaMap);
        end
        for curInteractionMaskIdx = 1:interactionMaskCount
            
            % get the current mask
            curInteractionMask = interactionMasks{curInteractionMaskIdx};
            
            % create the filename for the current nifti
            curFeatureFilename = sprintf('betas-%s-Bin%i.nii', curInteractionMask.name, curBin);
            curFeatureFilename = fullfile(niftisDir, curFeatureFilename);
            
            % the interaction masks provide a 0&1 mask of the features to
            % sum together to determine the beta for a single event from an
            % interaction term.
            maskBetas = curBinBetas(curInteractionMask.mask==1,:);
            summedBetas = sum(maskBetas,1);

            % store the summed beats for the current event 
            curEventBetas(curInteractionMaskIdx,voxIdxs) = summedBetas;
            
            % allocate a matrix to store the betas for the entire volume
            curInteractionMap = zeros(1,volumeVoxCount);
            curInteractionMap(voxIdxs) = summedBetas;
            
            % write out the current feature's betas
            writeVol(curFeatureFilename, header, curInteractionMap);
        end
        
        % copt the event beats for this bin into the master matrix
        eventBetasStart = (curBin-1)*interactionMaskCount + 1;
        eventBetasEnd = curBin*interactionMaskCount;
        eventBetas(eventBetasStart:eventBetasEnd,:) = curEventBetas;
    end    
    
    % write out the betas
    betasFilename = fullfile(niftisDir, 'betas.mat');
    save(betasFilename, 'betas');
    
    % get the correlation of the corss-validation sets
    ccEst = concatSubFields('model','cc',1,estDir);
    
    % get the indices of the lambdas used for each voxel and then get the
    % correlations used for each voxel
    lambdasUsed = concatSubFields('model','lambdasUsed', 1, estDir);
    ccEstUsed = ccEst(sub2ind(size(ccEst),[1:length(lambdasUsed)]',lambdasUsed));
            
    % allocate a matrix to store the betas for the entire volume
    ccEstUsedVol = zeros(1,volumeVoxCount);
    ccEstUsedVol(voxIdxs) = ccEstUsed;
               
    % write out the estimation correlation coefficients
    ccEstFilename = 'correlation-Est.nii';
    ccEstFilename = fullfile(niftisDir, ccEstFilename);
    writeVol(ccEstFilename, header, ccEstUsedVol);
    
    if mrifOptions.writeMeanVar          
       for curFIRBin = mrifOptions.binsFIR          
            % get the variance accounted for in the validation set
            fieldNameMean = sprintf('mean_bin%i', curFIRBin);
            meanBin = concatSubFields('model', fieldNameMean, 2, snrDir);
            fieldNameVar = sprintf('var_bin%i', curFIRBin);
            varBin = concatSubFields('model', fieldNameVar, 2, snrDir);

            % create the means and vars folders if they don't exist
            meansFolder = fullfile(niftisDir, 'means');
            varsFolder = fullfile(niftisDir, 'vars');
            if ~exist(meansFolder, 'dir')
                mkdir(meansFolder)
            end
            if ~exist(varsFolder, 'dir')
                mkdir(varsFolder)
            end
            
            % iterate through all the stimuli that were shown and write out
            % a 4D data file for each
            stimCount = size(meanBin,1);
            for curStim = 1:stimCount
            
                % create the filenames for mean and variance images
                meanBinFilename = sprintf('mean-Val_bin%i_stim%03i.nii', curFIRBin, curStim);
                meanBinFilename = fullfile(niftisDir, 'means', meanBinFilename);
                varBinFilename = sprintf('var-Val_bin%i_stim%03i.nii', curFIRBin, curStim);
                varBinFilename = fullfile(niftisDir, 'vars', varBinFilename);
            
                % allocate a matrix to store the betas for the entire volume
                meanBinVol = zeros(1,volumeVoxCount);
                meanBinVol(voxIdxs) = meanBin(curStim,:);
                varBinVol = zeros(1,volumeVoxCount);
                varBinVol(voxIdxs) = varBin(curStim,:);

                % write out the files
                writeVol(meanBinFilename, header, meanBinVol);
                writeVol(varBinFilename, header, varBinVol);
            end
       end
    end
    
    % if validation averaging was used then there will be one cc for each
    % FIR bin, and a one cc for each bin x compValFeature, so pull all of
    % those out and write them
    if mrifOptions.averageValRuns
        for curFIRBin = mrifOptions.binsFIR          
            
            % get the expainable variance
            expVarR2FieldName = sprintf('r2_bin%i', curFIRBin);
            expVarR2 = concatSubFields('model', expVarR2FieldName, 1, snrDir);
            
            % allocate a matrix to store the betas for the entire volume
            expVarVolR2 = zeros(1,volumeVoxCount);
            expVarVolR2(voxIdxs) = expVarR2;
            expVarVolR = sqrt(expVarVolR2);
            
            % write out the estimation correlation coefficients
            expVarFilenameR2 = 'explainableVariance_R2.nii';
            expVarFilenameR2 = fullfile(niftisDir, expVarFilenameR2);
            writeVol(expVarFilenameR2, header, expVarVolR2);
            expVarFilenameR = 'explainableVariance_R.nii';
            expVarFilenameR = fullfile(niftisDir, expVarFilenameR);
            writeVol(expVarFilenameR, header, expVarVolR);
            
            % get the variance accounted for in the validation set
            fieldNameBin = sprintf('cc_bin%i', curFIRBin);
            ccValBin = concatSubFields('model', fieldNameBin, 2, valDir);
            
            % allocate a matrix to store the betas for the entire volume
            ccValBinVol = zeros(1,volumeVoxCount);
            ccValBinVol(voxIdxs) = ccValBin;
            
            % write out the estimation correlation coefficients
            ccValBinFilename = sprintf('correlation-Val_bin%i.nii', curFIRBin);
            ccValBinFilename = fullfile(niftisDir, ccValBinFilename);
            writeVol(ccValBinFilename, header, ccValBinVol);

            % now scale the correlation coefficients by the variance map. Variance
            % here is explainable variance, so this is a way to normalize the
            % correlation coefficients to show how well they did compared to what
            % the data allows (i.e. how good our signal is)
            ccValBinNorm = ccValBinVol ./ expVarVolR;
            ccValBinNormFilename = sprintf('correlationNormed-Val_bin%i.nii', curFIRBin);
            ccValBinNormFilename = fullfile(niftisDir, ccValBinNormFilename);
            writeVol(ccValBinNormFilename, header, ccValBinNorm);
            
%             
%             % iterate through all the levels of the feature that the
%             % predictions were split upon
%             for curCompValFeature = mrifOptions.valCompFeatures
%           
%                 % get the variance accounted for in the validation set
%                 fieldNameValBinComp = sprintf('cc_bin%i_valCompFeature%04i', curFIRBin, curCompValFeature);
%                 ccValBinComp = concatSubFields('model', fieldNameValBinComp, 2, valDir);
%                 
%                 % allocate a matrix to store the betas for the entire volume
%                 ccValBinCompVol = zeros(1,volumeVoxCount);
%                 ccValBinCompVol(voxIdxs) = ccValBinComp;
%                 
%                 % write out the estimation correlation coefficients
%                 ccValBinCompFilename = sprintf('correlation-Val_bin%i_%s.nii', curFIRBin, featureNames{curCompValFeature});
%                 ccValBinCompFilename = fullfile(niftisDir, ccValBinCompFilename);
%                 writeVol(ccValBinCompFilename, header, ccValBinCompVol);
%             end
        


            % make a histogram of the correlation coefficients of all voxels
            ccValBinFlat = reshape(ccValBinVol, prod(size(ccValBinVol)), 1);
            f = figure('Visible','off');
            set(gca, 'XScale', 'log');
            hist(ccValBinFlat,50);
            ccValFlatFilename = fullfile(niftisDir, sprintf('correlationHistogram-Val_bin%i.png',curFIRBin));
            print('-dpng', ccValFlatFilename);
        end
    else
    
        % get the variance accounted for in the validation set
        ccVal = concatSubFields('model', 'cc', 2, valDir);
        
        % allocate a matrix to store the betas for the entire volume
        ccValVol = zeros(1,volumeVoxCount);
        ccValVol(voxIdxs) = ccVal;
        
        % write out the estimation correlation coefficients
        ccValFilename = 'correlation-Val.nii';
        ccValFilename = fullfile(niftisDir, ccValFilename);
        writeVol(ccValFilename, header, ccValVol);

        % make a histogram of the correlation coefficients of all voxels
        ccValFlat = reshape(ccValVol, prod(size(ccValVol)), 1);
        f = figure('Visible','off');
        set(gca, 'XScale', 'log');
        hist(ccValFlat,100);
        ccValFlatFilename = fullfile(niftisDir, 'correlationHistogram-Val.png');
        print('-dpng', ccValFlatFilename);
        
        % write out the histogram data to a mat file
        ccValFlatMatFilename = fullfile(niftisDir, 'correlationHistogram-Val.mat');
        save(ccValFlatMatFilename, 'ccValFlat');

        % if explainable variance maps were made, save those out along with
        % the normed correlation val
        if find(ismember(mrifOptions.modes, 'snr'))
            % get the expainable variance
            expVarCC = concatSubFields('model', 'cc', 1, snrDir);
            expVarR2 = concatSubFields('model', 'r2', 1, snrDir);
            expVarMSE = concatSubFields('model', 'mse', 1, snrDir);
            
            % get the indices of the lambdas used for each voxel and then get the
            % correlations used for each voxel
            snrLambdasUsed = concatSubFields('model','lambdasUsed', 1, snrDir);
            expVarUsedCC = expVarCC(sub2ind(size(expVarCC),[1:length(snrLambdasUsed)]',snrLambdasUsed));
            expVarUsedR2 = expVarR2(sub2ind(size(expVarR2),[1:length(snrLambdasUsed)]',snrLambdasUsed));
            expVarUsedMSE = expVarMSE(sub2ind(size(expVarMSE),[1:length(snrLambdasUsed)]',snrLambdasUsed));
            
            % allocate a matrix to store the betas for the entire volume
            expVarVolCC = zeros(1,volumeVoxCount);
            expVarVolCC(voxIdxs) = expVarUsedCC;
            expVarVolR2 = zeros(1,volumeVoxCount);
            expVarVolR2(voxIdxs) = expVarUsedR2;
            expVarVolR = sqrt(expVarVolR2);
            expVarVolMSE = zeros(1,volumeVoxCount);
            expVarVolMSE(voxIdxs) = expVarUsedMSE;
            
            % write out the estimation correlation coefficients
            expVarFilenameCC = 'explainableVariance_CC.nii';
            expVarFilenameCC = fullfile(niftisDir, expVarFilenameCC);
            writeVol(expVarFilenameCC, header, expVarVolCC);
            expVarFilenameR2 = 'explainableVariance_R2.nii';
            expVarFilenameR2 = fullfile(niftisDir, expVarFilenameR2);
            writeVol(expVarFilenameR2, header, expVarVolR2);
            expVarFilenameR = 'explainableVariance_R.nii';
            expVarFilenameR = fullfile(niftisDir, expVarFilenameR);
            writeVol(expVarFilenameR, header, expVarVolR);
            expVarFilenameMSE = 'explainableVariance_MSE.nii';
            expVarFilenameMSE = fullfile(niftisDir, expVarFilenameMSE);
            writeVol(expVarFilenameMSE, header, expVarVolMSE);
        
            % get the r-value at which the noise ceiling map is considered to
            % be significant, not noise
            sigCCaArr = concatSubFields('model', 'sigCC', 1, snrDir);
            sigCC = sigCCaArr(1);
            
            % make a thresholded explainable variance map, that exclude (sets
            % to nan) everything that is less significant than a given value.
            % This is because for voxels that are essentially noise, you can
            % have very large looking normed prediction accuracy that is the
            % result of just a very small explainable variance.
            noiseCeilingMap = expVarVolCC;
            noiseCeilingMap(noiseCeilingMap<sigCC) = nan;
            
            % now scale the correlation coefficients by the variance map. Variance
            % here is explainable variance, so this is a way to normalize the
            % correlation coefficients to show how well they did compared to what
            % the data allows (i.e. how good our signal is)
            ccValNorm = ccValVol ./ noiseCeilingMap;
            ccValNormFilename = 'correlationNormed-Val.nii';
            ccValNormFilename = fullfile(niftisDir, ccValNormFilename);
            writeVol(ccValNormFilename, header, ccValNorm);
        end
    end

    % Sam ToDo - open one of the validation mat files to get the number of
    % time points in the validation set (from the pred field), then create
    % the p-value map and write it out
    
    % make a map of the p-values for the associated validation correlation
    % coefficients
%    pValVol = r2pVal(ccValVol,
    
%     
%     % get the varianece and mean images from the snr folder
%     meanMap = concatSubFields('snr', 'meanMap',2,snrDir);
%     varMap = concatSubFields('snr', 'varMap',2,snrDir);
%     snrMap = concatSubFields('snr', 'snrMap',2,snrDir);
%             
%     % allocate a matrix to store the betas for the entire volume
%     meanMapVol = zeros(1,volumeVoxCount);
%     meanMapVol(voxIdxs) = meanMap;
%     varMapVol = zeros(1,volumeVoxCount);
%     varMapVol(voxIdxs) = varMap;
%     snrMapVol = zeros(1,volumeVoxCount);
%     snrMapVol(voxIdxs) = snrMap;
% 
%     % write out the mean, var and snr maps
%     meanMapFilename = fullfile(niftisDir, 'meanMap_Val.nii');
%     writeVol(meanMapFilename, header, meanMapVol);
%     varMapFilename = fullfile(niftisDir, 'varMap_Val.nii');
%     writeVol(varMapFilename, header, varMapVol);
%     snrMapFilename = fullfile(niftisDir, 'snrMap_Val.nii');
%     writeVol(snrMapFilename, header, snrMapVol);
%     

    % write out the top voxel histogram of betas if the flag says to
    if writeTopVoxels
        plotTopVoxels(estDir, valDir, niftisDir, featureNames, header.dim, binsFIR, TOP_VOXEL_COUNT);
    end
end
end

function writeVol(filename, header, data)
    
    % update the filename in the header
    header.fname = filename;
    header.private.dat.fname = filename;
    
    % get the data and reshape it to the dimensions of the volume
    vol = reshape(data,header.dim(1),header.dim(2),header.dim(3));
    
    % write out the nifti
    spm_write_vol(header, vol);    
end