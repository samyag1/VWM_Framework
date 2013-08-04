function stimIDMapFilename = wrapperCreateParadigms(options)

% the value in the frames vector that indicates no stimilus. make it 0 for
% ease
IGNORE_IDX = 0;
STIM_WINDOW=0; % number of frames to collapse into one stimulus

% calculate the minimum presentation frequency needed for the stimuli and
% ISI
options.presHz = 1 / min(options.stimOnTime, options.ISI);

% calculate how many frames occur per TR
framesPerTR = options.TR * options.presHz;
framesPerTrial = (options.stimOnTime + options.ISI)*options.presHz;

% create a dictionary to store the mapping between experiment relative stimIDs and filenames
stimIDMap = Dictionary();

% keep a count of how many stim have been used to make the stimIDs relative
% to the experiment and not the runs, as they are in the stimOrder files
previousStimuli = 0;

% determine the runs that we want to use
allRuns = [options.estRuns options.valRuns];

% iterate through all Sessions
newFilesWritten = false;
for curSess = 1:options.sessCount
    
    % iterate through all the runs in the sessions
    for curRunIdx = 1:numel(allRuns)
        
        curRun = allRuns(curRunIdx);
        
        % determine if this is an est or val run
        if any(options.estRuns==curRun)
            runType = options.estimationType;
            dummyScans = options.dummyScansEst;
        else
            runType = options.validationType;
            dummyScans = options.dummyScansVal;
        end
        
        % if the sessions are to be flattened, then always use session one
        % and calculate the run number
        outputSession = curSess;
        outputRun = curRun;
        if options.flattenSessions
            outputSession = 1;
            outputRun = (curSess-1)*options.runCount + curRun;
        end
        
        % create the paradigm's filename
        paradigmFilename = sprintf('paradigm_%s_Session%i_Run%02i.mat', runType, outputSession, outputRun);
        paradigmFilename = fullfile(options.paradigmFolder, paradigmFilename);

        % if the file already exists, then move on to the next one
        if exist(paradigmFilename, 'file')
            continue;
        end
        
        % create the stimOrder filename
        stimOrderFilename = sprintf(options.stimOrderFilenameTemplate, runType, curSess, curRun);
        stimOrderFilename = fullfile(options.stimFilesFolder, stimOrderFilename);
        
        % create the fileList filename
        fileListFilename = sprintf(options.fileListFilenameTemplate, runType, curSess, curRun);
        fileListFilename = fullfile(options.stimFilesFolder, fileListFilename);
        
        % read in the stimOrder file
        stimOrder = csvread(stimOrderFilename);
        
        % read in the fileList file. Use csvimport instead of csvread
        % because it is a CSV file of text strings and matlab cries with
        % anything but numbers
        fileList = csvimport(fileListFilename, 'noheader', true);
        
        % turn the list of stimIDs into a frames list at 1hz by putting in
        % the ignoreID (0) into the second second of the TR and the 2 seconds
        % of the no-stim TR
        newStimCount = 0;
        stimCount = numel(stimOrder);
        stimFrameCount = stimCount*framesPerTrial;
        preDummyFrameCount = dummyScans(1)*framesPerTR;
        postDummyFrameCount = dummyScans(2)*framesPerTR;
        frames = zeros(preDummyFrameCount + stimFrameCount + postDummyFrameCount,1);
        for curStim = 1:stimCount
            stimIdx = (curStim-1)*framesPerTrial + preDummyFrameCount + 1;
            
            % since the stimIDs in the stimOrder file are relative to the
            % fileList for that run, and we want the stimIDs to be relative
            % to the experiment for the paradigms and features matrix, add
            % the count of previous unique stimIDs to the current stimID
            runRelativeStimID = stimOrder(curStim);
            if runRelativeStimID == 0
                frames(stimIdx) = IGNORE_IDX;
            else
                % get the filename of the current stim
                curFilename = fileList{runRelativeStimID};
                
                % if this stim filename has been used previously (in
                % another run), then reuse it's experiment relative ID
                if stimIDMap.containsKey(curFilename)
                    experimentRelativeStimID  = stimIDMap(curFilename);
                else
                    experimentRelativeStimID = runRelativeStimID + previousStimuli;
                    newStimCount = newStimCount + 1;
                end
                
                % store the experiment relative StimID in the frames
                frames(stimIdx) = experimentRelativeStimID;
                
                % store the filename in the stimID map
                stimIDMap(curFilename) = experimentRelativeStimID;
            end
        end
        
        % create a paradigm for the current run
        paradigm = framesToParadigm(frames, options.presHz, options.TR, dummyScans, STIM_WINDOW, IGNORE_IDX);
        
        % write out the paradigm to file
        save(paradigmFilename, 'paradigm');
        
        % update the flad indicating a file has been written
        newFilesWritten = true;
        
        % now update the previous stim count
        previousStimuli = previousStimuli + newStimCount;
    end
end

% create the file path for the stimIDMap file
stimIDMapFilename = fullfile(options.paradigmFolder, options.stimIDMapFilename);

% only write out the stimID map if new paradigm files have been written
if newFilesWritten
    % write out the stimIDMap
    save(stimIDMapFilename, 'stimIDMap');
end

end