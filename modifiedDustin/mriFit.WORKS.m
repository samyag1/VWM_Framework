function [opt] = mriFit(experiment,subj,sessIdx,features,opt,voxIdx)
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

% IF ON SLURM, HOW ARE WE ALLOCATING RESOURCES?
try
	slurmPartition = opt.slurmPartition;
catch
	slurmPartition = 'all';
end

if opt.onSlurm
     slurmDir = fullfile(experiment.expDir,subj,'slurm');

	if ~exist(slurmDir,'dir')
		mkdir(slurmDir);
		mkdir(fullfile(slurmDir,'scripts'));
		mkdir(fullfile(slurmDir,'variables'));
		mkdir(fullfile(slurmDir,'logs'));
	end
end



% ALLOW FOR SUBSETS OF POSSIBLE CHUNKS
try
	chunks = opt.chunks
	nChunks = numel(chunks);
catch
	chunks = [];
end
opt.experimentName = experiment.expName;
% TYPE OF DATASET WE'RE WORKING WITH
% (GALLANT LAB STILL IMAGES PARADIGM OR THE LIKE)
switch lower(opt.mode)
     case {'estamp','est'}
          opt.condIdx = 1;
          opt.mode = 'est';

    case {'val','valresp'}
          opt.condIdx = 2;
          opt.mode = 'val';

     otherwise
          opt.condIdx = 1;
          opt.mode = 'est'
          fprintf('\nWARNING: Unkown processing mode.\n');   
end

sessCnt = 1;
mode0 = opt.mode;
jobIDs = [];
% LOOP OVER SESSIONS
for sess = sessIdx
	fprintf('\n\nProcessing session %d (%d out of %d)\n', sess,sessCnt,numel(sessIdx));
	% GET CURRENT XPMT
	xpmt = experiment.(subj){sess};
	opt.session = sess;
	if notDefined('chunks')
		nVox = prod(xpmt.volumeDim{opt.condIdx});

		% ONLY ANLAYZE VOXEL IN VOXIDX, IF PROVIDED
		if ~isempty(voxIdx)
			fprintf('\n***Note: Fitting a subset of total voxels***\n');
			chunks = segmentVector(voxIdx,opt.chunkSize);
		else
			chunks = colonWin(1,opt.chunkSize,nVox);
		end
		nChunks = numel(chunks);
	end
	
	if strcmp(lower(opt.mode),'est');
		if sessCnt == 1
			fprintf('\n------------------------------------------------');
			fprintf('\n       ***Estimating Feature Weights***');
			fprintf('\n------------------------------------------------');
		end
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
						slurmID = ['mrifRun_',opt.mode,'_Sess0',num2str(sess),'_',sprintf('%04d',chunk)];
						saveFile = fullfile(slurmDir,'variables',[slurmID,'.mat']);
						save(saveFile,'xpmt','chunkIdx','opt','features');
						jobID = mrifRunParallel(slurmDir,slurmID,saveFile, chunk);
						jobIDs = [jobIDs jobID];
				else
						opt = mrifRun(xpmt,chunkIdx,features,opt);
				end
			end
		end% (END CHUNK LOOP)
	end
	sessCnt = sessCnt+1;
end % END SESSION LOOP


if opt.onSlurm && ~isempty(jobIDs)
	disp(jobIDs);
	waitForQueue(jobIDs);
end

if strcmp(opt.mode,'est')
	fprintf('\n--------------------------------------');
	fprintf('\n Feature Weight Estimation Complete! ');
	fprintf('\n--------------------------------------\n');
	opt.mode = 'val';
end

% VALIDATION RESPONSE ESTIMATION
if strcmp(lower(opt.mode),'val'); 
	opt.condIdx = 2;
	jobIDs = [];
	sessCnt = 1;
	fprintf('\n--------------------------------------------------------------');
	fprintf('\n    ***Calculating validation data response amplitudes***');
	fprintf('\n--------------------------------------------------------------');
	for sess = sessIdx
		opt.session = sess;
		xpmt = experiment.(subj){sess};
		fprintf('\n\nProcessing session %d (%d out of %d)\n', sess,sessCnt,numel(sessIdx));
		% LOOP OVER RUNS
		for chunk = 1:nChunks
			fprintf('\n-----------------------------------');
			fprintf('\n      Chunk %d out of %d ',chunk,nChunks);
			fprintf('\n----------------------------------');
			chunkIdx = chunks{chunk};
			opt.chunkNum = chunk;
			if ~isempty(chunkIdx)
				if opt.onSlurm
					% GET SLURM ID; SAVE INPUTS; SEND TO SLURM
					slurmID = ['mrifRun_',opt.mode,'_Sess0',num2str(sess),'_',sprintf('%04d',chunk)];
					saveFile = fullfile(slurmDir,'variables',[slurmID,'.mat']);
					save(saveFile,'xpmt','chunkIdx','opt','features');
					jobID = mrifRunParallel(slurmDir,slurmID,saveFile, chunk);
					jobIDs = [jobIDs jobID];
				else
					opt = mrifRun(xpmt,chunkIdx,features,opt);
				end
			end
		end
		sessCnt = sessCnt+1;
	end
	
	if opt.onSlurm && ~isempty(jobIDs)
		disp(jobIDs);
		waitForQueue(jobIDs);
	end
	fprintf('\n-----------------------------------');
	fprintf('\n Response Estimation Complete!');
	fprintf('\n-----------------------------------\n');
else
	error('unknown mode!')
end

% DISPLAY COMPLETION
try
	endTime = toc(startTime);
	fprintf('\nFinished. Analysis duration: %1.2f minutes.\n',endTime/60);
catch
	fprintf('\nFinished.\n')
end


% SLURM-QUEUEING FUNCTIONS
%------------------------------
function [jobID, comStr] = mrifRunParallel(slurmDir,slurmID,saveFile, chunkNo)

script = [slurmID,'.script'];
scriptFile = [slurmDir,'/scripts/',script];

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
MEMORY_NEEDED = '4500';
logFilename = ['/logs/slurm_chunk-' num2str(chunkNo) '-%j.out '];

%DEFINE COMMANDS
batchStr =  ['sbatch --mem ',MEMORY_NEEDED,' -x krokodil ',' -o ',slurmDir,logFilename,scriptFile];
comstr0 = '#!/bin/bash';
comstr1 = ['nice -n +4 /auto/k5/matlab714/bin/matlab -singleCompThread -r "'];
comstr2 = ['addpath(genpath(''/auto/k2/spm8''));addpath(genpath(''/auto/k1/samyag1/code/''));addpath(''/auto/k1/samyag1/IAPS/analysis/scripts/modifiedDustin/'');'];
comstr3 =['load ',sprintf('%s',saveFile),';'];
comstr4 = 'opt = mrifRun(xpmt,chunkIdx,features,opt);exit"';

% WRITE THE SCRIPT
fid = fopen(scriptFile,'w');
fprintf(fid,'%s\n%s%s%s%s',comstr0,comstr1,comstr2,comstr3,comstr4);
fclose(fid);	

% MAKE SURE EVERYTING IS READY TO BE READ
pause(2); 

% RUN COMMAND AND GET JOB ID
[status, message] = system(batchStr);
message = regexp(message,' ','split');
jobID = str2double(message{end});

%------------------------------------------------------------------------
function waitForQueue(jobIDs)
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
        stat = checkQueue(jobIDs);
        pause(checkPeriod*.9);  % WAIT BEFORE START CHECKING TIME AGAIN
        if stat
            waiting = 0;
        end
    end
end

%---------------------------------------------------------
function status = checkQueue(jobIDs)
numJobs = length(jobIDs);
status = 0;
passCnt = 0;
for job = 1:numJobs
    checkStr = ['squeue --jobs ',num2str(jobIDs(job))];
    [stat, msg]=system(checkStr);
    if length(msg) < 90  % KINDA HACKY, BUT WORKS
        passCnt = passCnt+1;
    end
end
if passCnt == numJobs
    status = 1;
end
