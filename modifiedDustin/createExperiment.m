function experiment=createExperiment(experiment,subj,dataDirs,runIdxs,stimFiles,preFixes,preprocPostfix,expDir,expName,modelName,voxSize, voxDims, TR)
%-----------------------------------------------------------------
% experiment=createExperiment(experiment,subj,dataDirs,runIdxs,preFixes,stimFiles,expDir)
%------------------------------------------------------------------
%  Creates a data profile entry for subject <subj> in the experimental
%  datastrucure <experiment>. 
%  
%  INPUTS:
%    <experiment>:  a datastructure containing all information that
%                   would added in this function. an <experiment>
%                   datastructure is
%                   created if one is not provided. 
%
%    <subj>:        a string subject identifier whose data is added to
%                   <experiment>. e.g 'DS'
%
%    <dataDirs>:    cell array that is numSessions in length containing
%                   strings of fullfile directory paths to dicoms for each
%                   session
%
%    <runIdxs>:      cell array numSessions in length containing K vectors
%                   of indices for K different condition runs recorded
%                   in each session.
%
%    <stimFiles>:   fullpath string to location of onoff.mat matrix file
%                   for the an experimental paradigm condition
%
%    <preFixes>:    (optional) a cell array of strings indicating the scan
%                   protocol prefix.  e.g. {'ep_fatsat_32ch_'} or
%                   {'ep_fatsat_32ch_, 'ep_2d_fatsat_'}. IF not provided,
%                   the function assumes that there is only one type of
%                   protocol being used in each session's data directory.
%                   Therefore, if more than one protocol is being used, then
%                   this will throw an error.
%
%<preprocPostfix>:  (optional) a string indicating the the postfix to append
%                   to directories to indicate preprocessing. default = 
%                   '.preproc'.
%
%    <expDir>:      (optional) a fullpath to a directory for which to store
%                   all preceding processing metadata for current experimen.
%                   Default is '/auto/data/archive/mri/experiments/', which
%                   can be written by all of the gallant group members.
%  
%  OUTPUTS:
%    <experiment>:  same datastructure provided with data profile for the
%                   provided subject added. ie
%                   experiment.(subj){1:numSessions}.
%
%-----------------------------------------------------------------------
% DES  

%TODO:
% ENSURE VARIABLES DIMENSIONS ARE CONSISTENT
% -IE runIdx must be the same length as <prefix>
% REMOVE NEED TO PROVIDE info.txt
%  -mainly used to find scan protocol, but this can
%  probably be done by some clever directory parsing
if nargin == 0
    help createDataProfile
    return
end

warning('off')

if notDefined('expDir');
     expDir = '/auto/data/archive/mri/experiments/';
end

if ~isdir(expDir)
	mkdir(expDir)
end

% define the directory to put all the output as the experiment's directory
% followed by the model name then subject name
outputDir = fullfile(expDir, modelName, subj); 

if notDefined('experiment');
	experiment = {};
	experiment.expName = expName;
    experiment.modelName = modelName;
    experiment.expDir = fullfile(expDir,modelName); 
end

% allow an empty preprocPostfix for the case that preprocessing was done
% using another system
experiment.preprocPostfix = preprocPostfix;
% if notDefined('preprocPostfix')
% 	preprocPostfix = '.preproc';
%     experiment.preprocPostfix = '.preproc';
% else
%     experiment.preprocPostfix = preprocPostfix;
% end

fprintf('\nNOTE: All preprocessing will fall under a ''<DATA>''%s directory\n',experiment.preprocPostfix);

if length(runIdxs) ~= length(dataDirs)
    error('Dimensions of run indexes and dataDirs do not correspond!');
end

if ~isfield(experiment,subj)
    experiment.(subj) = {};
end

numSessions = length(dataDirs);

if length(stimFiles) == 1
    for sess = 1:numSessions
	   temp{sess} = stimFiles{:};
    end
    stimFiles = temp;    
end


% PARSE SCAN PROTOCOL FROM DATA DIRECTORY
if notDefined('preFixes')
	dirIdx = [];
	for iD = 1:numel(dataDirs)
		preFixes{iD} =repmat({identifyScanProtocol(dataDirs{iD})},numel(runIdxs{iD}),1);
	end
end

% LOOP OVER SESSION INFO
for sess = 1:numSessions
   
    if strcmp(dataDirs{sess}(end),'/')
	   dataDir{sess} = dataDirs{sess};
    else
	   dataDir{sess} = [dataDirs{sess}, '/'];
    end
	% CHECK DIMENSIONS OF PROVIDED preFixes
	if numel(preFixes{sess}) == numel(runIdxs{sess})
		prefixes{sess} = preFixes{sess};
	else
		disp('Dimension of prefix doesn''t match runIdx')
		error('May need to duplicate a prefix in setup!');
	end    
end

for sess = 1:numSessions
    exp = {};
    exp.dataDir = dataDir{sess};
    exp.prefix = prefixes{sess};
    exp.sessID = ['session',sprintf('%02d',sess)];
    exp.runIdx = runIdxs{sess};  
    
    if ~notDefined('preprocPostfix')
        exp.preprocDir = [fullfile([exp.dataDir(1:end-1),preprocPostfix]),'/'];
    end
    
    % SAM this is where the model output will be saved
    exp.outputDir = outputDir;
    sessionPath = fullfile(exp.outputDir, 'modeling', exp.sessID);
    exp.sessionDir = sessionPath;
    exp.estDir = fullfile(sessionPath, 'est');
    exp.valDir = fullfile(sessionPath, 'val');
    exp.snrDir = fullfile(sessionPath, 'snr');
    exp.niftisDir = fullfile(sessionPath, 'niftis');
    if ~exist(sessionPath,'dir')
        mkdir(sessionPath)
    end
    if ~exist(exp.estDir,'dir')
        mkdir(exp.estDir);
    end
    if ~exist(exp.valDir,'dir')
        mkdir(exp.valDir);
    end
    if ~exist(exp.snrDir,'dir')
        mkdir(exp.snrDir);
    end
    if ~exist(exp.niftisDir,'dir')
        mkdir(exp.niftisDir);
    end
   
    for cond = 1:length(exp.runIdx)
		try
			if isstr(stimFiles{sess}{cond})
				stimFiles{sess}{cond} = parseRunString(stimFiles{sess}{cond});
			end
			exp.stimFile{cond} = stimFiles{sess}{cond};
		catch
			exp.stimFile{cond} = 'unknown';
			fprintf('Warning: could not parse runfiles.\n');
        end
        if notDefined('voxDims')
            % PARSE INFORMATION FROM DICOMS
            dicoms = matchwildcards(fullfile(exp.dataDir, ...
                [exp.prefix{cond},num2str(exp.runIdx{cond}(1))],'/IM-*.dcm'));
            exp.dicom{cond}=dicominfo(dicoms{1});
            
            % VOLUME DIMENSIONS
            exp.voxDim{cond} = double([exp.dicom{cond}.PixelSpacing', exp.dicom{cond}.SpacingBetweenSlices]);
            numSlices = single(exp.dicom{cond}.Private_0019_100a);
            
            imageDim = exp.dicom{cond}.AcquisitionMatrix(find(exp.dicom{cond}.AcquisitionMatrix));
            exp.volumeDim{cond} = double([imageDim(1), imageDim(2), numSlices]);
            exp.TR{cond} = double(exp.dicom{cond}.RepetitionTime/1000);
        else
            exp.voxDim{cond} = voxSize;
            exp.volumeDim{cond} = voxDims;
            exp.TR{cond} = TR;
        end
        
		% GET # OF MEASUREMENTS (ASSUMING SAME # FOR EACH CONDITION)
		exp.nMeasure{cond} = length(dir([exp.dataDir,exp.prefix{cond},sprintf('%d',exp.runIdx{cond}(1)),'/*.dcm']));
		runFold = [];
	   
		for rn = 1:length(exp.runIdx{cond})
			runFold{rn} = fullfile(exp.dataDir,sprintf('%s%d',exp.prefix{cond},exp.runIdx{cond}(rn)));
		end
	   
		exp.runFold{cond} = [runFold(:)];
		exp.runFold0{cond} = [runFold(:)];
		if ~isempty(findstr(exp.stimFile{cond}{1},'localizer'))
			onSetType = '';
		else
			onSetType = 'Impulse';
		end
		try
			exp.paradigmFile{cond} = identifyParadigm(exp.stimFile{cond}{1},onSetType);
		catch
			exp.paradigmFile{cond} = 'unkown';
			fprintf('Warning: could not identify experiment paradigm\n');
		end
    end    
   % IDENTIFY PARADIGM FILES FROM STIMULUS FILES
    experiment.(subj){sess} = exp;
    
end
%experiment.expDir = [expDir,experiment.expName,filesep]; SAM - removed,
%don't want the experiment directory to have the experiment name appended
