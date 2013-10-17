function [opt] = mriFit(experiment,subj,sessIdx,features,featuresNuissanceEst,featuresNuissanceVal,opt,voxIdx)
%optOut= mriFit(exper,subj,sessIdx,features,optIn,voxIdx)
%----------------------------------------------------------------------------
% Models time series data for a set of voxels:
%----------------------------------------------------------------------------
% DES - stan_s_bury@berkeley.edu

if nargin == 0
    return
end

startTime = tic;

% PRELIMS...
if notDefined('experiment')
     error(sprintf('This analysis depends on an <experiment> data structure, which was not provided.\n(See createDataProfile.m'));
end

if notDefined('subj')
     error('No subject provided!!');
end

if notDefined('sessIdx')
     fprintf('\n***Note: no session index provided\nanalyzing all sessions***\n');
     sessIdx = 1:numel(experiment.(subj));
end

% CREATE OPTIONS STRUCTURE, IF NOT PROVIDED
if notDefined('opt')
     opt = mrifOpts;
     fprintf('\n***Note: no options structure provided, using default template.***\n');
end

% LOAD DEFAULT BASIS SET IF NOT SUPPLIED
if notDefined('features')
    error(sprintf('This analysis depends on a set of features corresponding to the stimuli.\n See mrifRun.m and mrifEstimate.m for details'))
end

if notDefined('voxIdx')
    voxIdx = [];
end

if opt.onSlurm
     queueDir = fullfile(experiment.expDir,subj,'queue');

	if ~exist(queueDir,'dir')
		mkdir(queueDir);
		mkdir(fullfile(queueDir,'scripts'));
		mkdir(fullfile(queueDir,'variables'));
		mkdir(fullfile(queueDir,'logs'));
	end
end

% ALLOW FOR SUBSETS OF POSSIBLE CHUNKS
try
	chunks = opt.chunks
	nChunks = numel(chunks);
catch
	chunks = [];
end
if notDefined('chunks')
    xpmt1 = experiment.(subj){sessIdx(1)};
    nVox = prod(xpmt1.volumeDim{1});
    
    % ONLY ANLAYZE VOXEL IN VOXIDX, IF PROVIDED
    if ~isempty(voxIdx)
        fprintf('\n***Note: Fitting a subset of total voxels***\n');
        chunks = segmentVector(voxIdx,opt.chunkSize);
    else
        chunks = colonWin(1,opt.chunkSize,nVox);
    end
    nChunks = numel(chunks);
end

% put the experiment's name ionto the options structure
opt.experimentName = experiment.expName;

% iterate over all the modes to run
for curModeIdx = 1:numel(opt.modes)
    
    opt.mode = lower(opt.modes{curModeIdx});
    jobIDs = [];
    sessCnt = 1;
        
    % display starting message and update cond Idx
    switch opt.mode
        case {'est'}
            fprintf('\n------------------------------------------------');
            fprintf('\n       ***Estimating Feature Weights***');
            fprintf('\n------------------------------------------------');
            opt.condIdx = 1;
        case {'val'}
            fprintf('\n--------------------------------------------------------------');
            fprintf('\n    ***Calculating validation data response amplitudes***');
            fprintf('\n--------------------------------------------------------------');
            opt.condIdx = 2;
        case {'snr'}
            fprintf('\n--------------------------------------------------------------');
            fprintf('\n    ***Creating Signal to Noise Maps***');
            fprintf('\n--------------------------------------------------------------');
            opt.condIdx = 2; %use the validation data to make the snr maps
        otherwise
            opt.condIdx = 1;
            opt.mode = 'est'
            fprintf('\nWARNING: Unkown processing mode.\n');
    end
    
    % LOOP OVER SESSIONS
    for sess = sessIdx
        fprintf('\n\nProcessing session %d (%d out of %d)\n', sess,sessCnt,numel(sessIdx));

        % GET CURRENT XPMT
        xpmt = experiment.(subj){sess};
        opt.session = sess;
        
        % LOOP OVER CHUNKS
        for chunk = 1:nChunks
            fprintf('\n-----------------------------------');
            fprintf('\n      Chunk %d out of %d ',chunk,nChunks);
            fprintf('\n-----------------------------------');

            chunkIdx = chunks{chunk};
            opt.chunkNum = chunk;
            
            if ~isempty(chunkIdx)
                if opt.onSlurm
                    % GET SLURM ID; SAVE INPUTS; SEND TO SLURM
                    queueID = ['mrifRun_',opt.mode,'_Sess0',num2str(sess),'_',sprintf('%04d',chunk)];
                    saveFile = fullfile(queueDir,'variables',[queueID,'.mat']);
                    save(saveFile,'xpmt','chunkIdx','opt','features','featuresNuissanceEst','featuresNuissanceVal');
                    jobID = mrifRunParallel(opt.queueType,queueDir,queueID,saveFile, chunk, opt.matlabCommand, opt.pathsToAdd);
                    jobIDs = [jobIDs jobID];
                else
                    opt = mrifRun(xpmt,chunkIdx,features,featuresNuissanceEst,featuresNuissanceVal,opt);
                end
            end
        end% (END CHUNK LOOP)
        
        sessCnt = sessCnt+1;
        
    end % END SESSION LOOP
    
    % wait for the jobs to finish
    if opt.onSlurm && ~isempty(jobIDs)
        disp(jobIDs);
        waitForQueue(opt.queueType,jobIDs);
    end
    
    % there are some tasks that require all the voxels to determine, such
    % as finding a single value of lambda for all voxels.
    mrifPostRun(xpmt, opt);
    
    % display starting message
    switch opt.mode
        case {'est'}
            fprintf('\n--------------------------------------');
            fprintf('\n Feature Weight Estimation Complete! ');
            fprintf('\n--------------------------------------\n');
        case {'val'}
            fprintf('\n-----------------------------------');
            fprintf('\n Response Estimation Complete!');
            fprintf('\n-----------------------------------\n');
        case {'snr'}
            fprintf('\n-----------------------------------');
            fprintf('\n Signal to Noise Map Creation Complete!');
            fprintf('\n-----------------------------------\n');
    end
end % END MODE LOOP

% DISPLAY COMPLETION
try
	endTime = toc(startTime);
	fprintf('\nFinished. Analysis duration: %1.2f minutes.\n',endTime/60);
catch
	fprintf('\nFinished.\n')
end



% SLURM-QUEUEING FUNCTIONS
%------------------------------
function [jobID, comStr] = mrifRunParallel(queueType,queueDir,queueID,saveFile, chunkNo, matlabCommand, pathsToAdd)

script = [queueID,'.script'];
scriptFile = [queueDir,'/scripts/',script];

%REMOVE THE SCRIPT IF IT ALREADY EXISTS...
if exist(scriptFile)
    system(['rm ',sprintf('%s',scriptFile)]);
end

% %CHECK PROVIDED SLURM PARTITION
% switch slurmPartition
% 	case {'all','i7','highprio','lowprio'}
% 	otherwise
% 		slurmPartition = 'all'
% end

% define the shell command that will launch the job based on the queue type
switch(queueType)
    case {'slurm'}
        MEMORY_NEEDED = '4500';
        logFilename = ['/logs/slurm_chunk-' num2str(chunkNo) '-%j.out'];
        logPath = fullfile(queueDir, logFilename);
        batchStr =  ['sbatch --mem ',MEMORY_NEEDED,' -x krokodil ',' -o ',logPath,' ',scriptFile];
        jobOptions = '';
    case {'SGE'}
        MEMORY_NEEDED = '15G';
%        batchStr =  ['qsub -j yes -i virtual_free=',MEMORY_NEEDED,' -o ',logPath,'-e ',logPath,scriptFile];
        logPath = fullfile(queueDir, 'logs', ['slurm_chunk-' num2str(chunkNo) '-$JOB_ID.out']);
        errorLogPath = fullfile(queueDir, 'logs', ['error_slurm_chunk-' num2str(chunkNo) '-$JOB_ID.out']);
        batchStr =  ['qsub ' ,scriptFile];
        jobOptions = ['#$ -o ',logPath,' -e ',logPath];
    otherwise
        error('Invalid queueType speicified while trying to parallelize. The value is: %s and the possible values are "slurm" or "SGE"',queueType);
end

% create the string that will add all the paths necessary to run this
comstr2 = '';
for curPath = 1:numel(pathsToAdd)
    comstr2 = [comstr2, sprintf('addpath(genpath(''%s''));', pathsToAdd{curPath})];
end

comstr0 = '#!/bin/bash';
comstr1 = [matlabCommand, ' -singleCompThread -r "'];
comstr3 =['load ',sprintf('%s',saveFile),';'];
comstr4 = 'opt = mrifRun(xpmt,chunkIdx,features,featuresNuissanceEst,featuresNuissanceVal,opt);exit"';

% WRITE THE SCRIPT
fid = fopen(scriptFile,'w');
fprintf(fid,'%s\n%s\n%s%s%s%s',comstr0,jobOptions,comstr1,comstr2,comstr3,comstr4);
fclose(fid);	

% MAKE SURE EVERYTING IS READY TO BE READ
pause(2); 

% RUN COMMAND AND GET JOB ID
[status, message] = system(batchStr);

switch(queueType)
    case {'slurm'}
        message = regexp(message,' ','split');
        jobID = str2double(message{end});
    case {'SGE'}
        words = strsplit(message, ' ');
        jobID = str2double(words{3});
end

%------------------------------------------------------------------------
function waitForQueue(queueType,jobIDs)
fprintf('\n\n*---------------------------------------------------------*')
fprintf('\n    Waiting for queued jobs to finish before proceeding ');
fprintf('\n*---------------------------------------------------------*\n')

checkPeriod = 30; % CHECK EVERY 30 SEC
waiting = 1;
tic
pause(1);
while waiting
    tt = mod(round(toc),checkPeriod);
    if tt == 0
        stat = checkQueue(queueType,jobIDs);
        if stat
            waiting = 0;
            break;
        end
        pause(checkPeriod*.9);  % WAIT BEFORE START CHECKING TIME AGAIN
    end
end

%---------------------------------------------------------
function status = checkQueue(queueType,jobIDs)
numJobs = length(jobIDs);
status = 0;
passCnt = 0;
for job = 1:numJobs
    
    switch(queueType)
        case {'slurm'}
            checkStr = ['squeue --jobs ',num2str(jobIDs(job))];
            [stat, msg]=system(checkStr);
            if length(msg) < 90  % KINDA HACKY, BUT WORKS
                passCnt = passCnt+1;
            end
        case {'SGE'}
            checkStr = ['qstat -j ',num2str(jobIDs(job))];
            [stat, msg]=system(checkStr);
            jobCompletedString =  'Following jobs do not exist';
            if strncmp(msg, jobCompletedString, numel(jobCompletedString)) % KINDA HACKY, BUT WORKS
                passCnt = passCnt+1;
            end
    end
end
if passCnt == numJobs
    status = 1;
end
