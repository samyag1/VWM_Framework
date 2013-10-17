function opt = mrifOpts(varargin);
%---------------------------------------------------------
%  FUNCTION  opt = mrifOpts(varargin);
%---------------------------------------------------------
% Defines values for MRI time series estimation.  
% Can be used in a number of ways:
%
% 1) Can be used simply to provide a template options 
%    structure, by supplying no inputs.
%
% 2) Allows manual editing of a field name in the structure.  
%    This can also happen in two ways: by either directly 
%    providing the value as in...
%
%>> optNew = defaultProcOpts(optOld,fieldname,value)
%
%    ...or, if only the fieldname is provided, you will be 
%    prompted for a a value...
%
%>> optNew = defaultProcOpts(optOld,fieldname)
%
%    The value for field '<fieldName>' is 'blah.'
%    Replace it with what?''
%
%----------------------------------------------------------
% AN EXPLANATION OF <opt>.settings FIELDS:
%
% <onSlurm>:   Flag for whether or to process jobs in parallel 
%              on the SLURM queue default = 1, means to process 
%              on slurm queue;
%
% <segChunk>:  is the size of chunks that are spread across jobs
%              default = 5000;  NOTE: depending on the number of 
%              voxels being fit, this value may or may not be 
%              compatable.
%
% <mode>:      is the current mode of fitTimeSeries.m that you 
%              will want to run. ie 'estimation', 'hrf', or 
%             'validation' default = 'estimation' is the first 
%              step of fitTimeSeries.m
%----------------------------------------------------------
% DES

%---------------------------------------------
% TIME SERIES FITTING FIELDS (opt.settings.)
%---------------------------------------------
%


defOpt =  struct('onSlurm', 1, ...
                 'queueType', 'slurm', ... % indicates which type of queue to use for paralellization. Can be either slurm (Gallant lab cluster) or SGE (NX Cluster) only applies when onSlurm is true (Didn't want to change all references of onSlurm to something like onQueue)
                 'matlabCommand', '', ... % This is the command that will execute matlab when paralellizing jobs on the queue. It should be a fully qualified path if it's not in the $PATH environment variable
                 'pathsToAdd', {{}}, ... % A cell array of paths that will be added to matlab when jobs are parallelized on the queue
                 'chunkSize', 5000, ...
                 'mode','esthrf', ...
                 'dataPre', 'RCMvolume', ...
                 'fwhm',[], ...
                 'detrendType', 'SG', ... % this can be 'SG' (Savitzky-Goly, sliding window polynomial), 'poly', or 'none'
                 'detrendPolyDeg',3, ...
                 'detrendWindow', 2, ... % in minutes
                 'zscore',true, ...
                 'slurmPartition','all', ...
                 'fitVal',0, ...
                 'processChunks',[], ...
                 'featureName','unknown', ...
                 'preprocFun', 'noPreproc', ...
                 'cacheData', 1, ...
                 'delays', 1:4, ...
                 'addBias', 0, ...
                 'estType', 'ols', ... % Currently this can be 'ols' (Ordinary Least Squares) or 'elasticNet'
                 'crossValFolds', 10, ... % does k-fold cross validation, this specifies k
                 'binsFIR', [3,4], ...     % this is list of the bin offsets to use for the FIR model
                 'predictionCompFeatures', [], ...
                 'averageValRuns', true, ...
                 'writeMeanVar', true, ...
                 'removeNuissanceFilenamesEst', {{}}, ... % filenames of nuissance regressor files (in each run folder) to regress out from data in a pre-model regression when estimating the model
                 'removeNuissanceFilenamesVal', {{}}, ... % filenames of nuissance regressor files (in each run folder) to regress out from data in a pre-model regression when validating the model
                 'modelNuissanceFilenames', {{}}, ... % filenames of nuissance regressors (in each run folder) whose values are to be entered as regressors of non-interest in the main model, and thus will share variance with regressors of interest
                 'averageCrossValFolds',false, ... % when true the prediction accuracy of the cross val folds during estimation (which also affects the noise ceiling calculation) will be averaged. When false one prediction value will calculated by concatonating all the predictions together. THe later is most likely what you wan, unless you need a variance measure over the predictions. 
                 'crossValSignificant',.05, ... % the p-value associated with the prediction accuracy (r-value) considered to be significant. This is how the neoise ceiling threhold used to determine which voxels are considered to hvae enough signal to care about.
                 'noiseCeilingIgnoreReps',[], ... % If our data has a non-stationarity then we may want to exclude certain repititions of the validation set when doing the noise ceililng calculation. This specifies the reps to withhold
                 'onlyUseStimVols',false, ... % indicates whether to only include volumes that are of a given offset (blwo) where stim are presented
                 'onlyUseStimVolsOffset',2, ...
                 'writeDesignMatrix',false, ... % whether to write out the design matrix to a file.
                 'useSingleLambda',true, ... % whether to select the ridge hyper parameter lambda per voxel, or select one for all the voexels
                 'excludeValFeatures',[]); % NOT TO BE SET, DONE INTERNALL. This is for nuissance regressors that are stim based and so have FIR binning. This tells the framework to include in the estimation model, but not the validation prediction. 
				  

%----------------------------------------
%         VARIOUS FUNCTIONALITIES
%----------------------------------------
switch nargin
    case 0 % SPIT OUT DEFAULT OPTIONS STRUCTURE
        opt = defOpt;
%      otherwise
%             error('input must be a string indicating processing mode (e.g ''train'' or ''val'')');
%         end
%          
    case 2 % REQUEST A REPLACEMENT VALUE FROM USER
    opt = varargin{1}; if ~isstruct(opt); error('<opt> must be a structure!');return;end
        field = varargin{2}; if ~isstr(field); error('<field> must be a string!');return;end
     
        fieldNames = fieldnames(opt);
        pass = 0;

        for fn = 1:length(fieldNames)

            if isfield(opt.(fieldNames{fn}),field)
                curValue = opt.(fieldNames{fn}).(field)
                if strcmp(class(curValue),'char')
                    fprintf('\nThe field ''%s'' has the current value: ''%s''\n ',field,curValue);
                else
                    fprintf('\nThe field ''%s'' has the current value: %d\n ',field,curValue);
                end
                fprintf('\nReplace this value with');
                value = input(sprintf('\n(press <Enter> to keep current value): '));

                if ~isempty(value)
                    opt.(fieldNames{fn}).(field) = value;
                end
                pass = 1;
                break
            end
        end
        if ~pass
            error(sprintf('Provided fieldname %s does not exist',field));
            return
        end
        
    case 3  % FIND AND REPLACE
        
        opt = varargin{1}; if ~isstruct(opt); error('<opt> must be a structure!');return;end
        field = varargin{2}; if ~isstr(field); error('<field> must be a string!');return;end
        value = varargin{3}; 

        fieldNames = fieldnames(opt);
        pass = 0;
        
        for fn = 1:length(fieldNames)
        
            if isfield(opt.(fieldNames{fn}),field)
                opt.(fieldNames{fn}).(field) = value;
                pass = 1;
                break
            end
        end
        if ~pass
            error(sprintf('Provided fieldname %s does not exist',field));
            return
        end
end


