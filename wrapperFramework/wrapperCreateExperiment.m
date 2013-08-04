function exper = wrapperCreateExperiment(subj, options)

% create the filename for the experiment
experimentFilename = fullfile(options.outputFolder,options.modelName,sprintf('experiment-%s.mat',date));

% if it already exists, then load it and return
if exist(experimentFilename, 'file')
    load(experimentFilename, 'exper');
    return;
end

% if the sessions are flattened, then figure out the estimation and
% validation run numbers
sessCount = options.sessCount;
estRuns = options.estRuns;
valRuns = options.valRuns;
if options.flattenSessions
    sessCount = 1;
    for curSession = 1:(options.sessCount-1)
        estRuns = [estRuns, (options.estRuns + (curSession*options.runCount))];
        valRuns = [valRuns, (options.valRuns + (curSession*options.runCount))];
    end
end

stimFiles = cell(1,sessCount);
dataDirs = cell(1,sessCount);
preFixes = cell(1,sessCount);
runIdxs = cell(1,sessCount);
for curSession = 1:sessCount
    
    % DEFINE THE DIRECTORIES THAT CONTAIN THE (DICOM) DATA FOR EACH SESSION
    % This will only be used to read the one dicom I copied into run01 of
    % each Session. The createExperiment method uses info from the dicom
    % header to populate the experiment structure.
    dataDirs{curSession} = fullfile(options.dicomFolder, subj, sprintf('Session%i',curSession));
    
    %___________________________________________________
    % DEFINE THE PROTOCOL PREFIXES OF THE DICOM FOLDERS
    % FOR EACH EXPERIMENTAL CONDITION.
    % NOTE HERE THAT THERE ARE TWO DIFFERENT PROTOCOL NAMES,
    % ONE FOR TRAINING, AND ONE FOR VALIDATION CONDITIONS
    % Again, this is only used to get the header info from the first dicom
    preFixes{curSession} =  {options.estDicomFolderPrefix,options.valDicomFolderPrefix};
    
    %___________________________________________________
    % DEFINE THE RUN INDICES FOR EACH SESSION AND CONDITION
    % HERE THERE IS ONE SESSION, WITH TWO CONDITIONS
    % A TRAINING AND A VALIDATION
    % (WITHOUT PHASE VOLUMES)
    % Again, this is only used to get the header info from the first dicom
    % The 1 here is the postfix on the dicom folder, so for example I named
    % them est_IAPS_1, so I use a 1 below. I also have to specify an empty
    % array for the validation condition to make createExperiment happy.
    runIdxs{curSession} = {[1],[2]};
    
    %_______________________________________________________
    % POINTERS TO STIMULUS FILES FOR EACH OF THE CONDITIONS
    % THESE ARE NOT NECESSARY IF YOU'RE SIMPLY RUNNING
    % PREPROCESSING, BUT ARE NEEDED IF FITTING HRFs
    % OR RESPONSE AMPLITUDES; IT'S BETTER TO JUST ADD THEM
    % AT THE BEGINNING.
    % assemble the estimation run folders
    estRunParadigms = cell(1,numel(estRuns));
    for curEstRun = 1:numel(estRuns)
        % the run folder is the nifti root path, plus subject, ession and run
        estRunParadigms{curEstRun} = fullfile(options.paradigmFolder, ...
                                              sprintf('paradigm_est_Session%i_Run%02i.mat', curSession, estRuns(curEstRun)));
    end
    % assemble the estimation run folders
    valRunParadigms = cell(1,numel(valRuns));
    for curValRun = 1:numel(valRuns)
        % the run folder is the nifti root path, plus subject, ession and run
        valRunParadigms{curValRun} = fullfile(options.paradigmFolder, ...
                                              sprintf('paradigm_val_Session%i_Run%02i.mat', curSession, valRuns(curValRun)));
    end
    
    % store the estimation and validation runs for the current run
    % into the subject's session's runFold field
    stimFiles{curSession} = {estRunParadigms, valRunParadigms};
end

%____________________________________________________
% DEFINE THE POSTFIX FOR THE PREPROCESSED DATA FOLDER
% LEAVE THIS BLANK INDICATING THE RUN FOLDERS ARE PREPROCESSED
preprocPostfix = ''; 

%____________________________________
% CREATE EXPERIMENTAL DATA STRUCTURE
exper = [];
exper = createExperiment(exper, ...
                         subj, ...
                         dataDirs, ...
                         runIdxs, ...
                         stimFiles, ...
                         preFixes, ...
                         preprocPostfix, ...
                         options.outputFolder, ...
                         options.experimentName, ...
                         options.modelName, ...
                         options.voxSize, ...
                         options.voxDims, ...
                         options.TR);

%___________________________________________
% SINCE I WON'T BE USING THE PREPROCESSING PIPELINE I NEED
% TO TELL THE EXPERIMENT WHERE MY PROCESSED (NIFTI) RUN FILES ARE LOCATED
% NOTE: I use the original sessionCount, est and val runs here because the
% data is organized by sessions and runs, but i use the out session for
% storage
masterEstRunFolders = {};
masterValRunFolders = {};
for curSession = 1:options.sessCount
    
    % assemble the estimation run folders
    estRunFolders = cell(1,numel(options.estRuns));
    for curEstRun = 1:numel(options.estRuns)
        % the run folder is the nifti root path, plus subject, ession and run
        estRunFolders{curEstRun} = fullfile(options.niftiFolder, ...
                                            subj, ...
                                            sprintf('Session%i',curSession), ...
                                            sprintf('run%02i', options.estRuns(curEstRun)));
    end
    % assemble the estimation run folders
    valRunFolders = cell(1,numel(options.valRuns));
    for curValRun = 1:numel(options.valRuns)
        % the run folder is the nifti root path, plus subject, ession and run
        valRunFolders{curValRun} = fullfile(options.niftiFolder, ...
                                        subj, ...
                                        sprintf('Session%i',curSession), ...
                                        sprintf('run%02i', options.valRuns(curValRun)));
    end

    % if the sessions are to be flattened, then store the run folders into
    % master lists to be added afterwards
    if options.flattenSessions
        masterEstRunFolders = [masterEstRunFolders, estRunFolders];
        masterValRunFolders = [masterValRunFolders, valRunFolders];
    else
        % store the estimation and validation runs for the current run
        % into the subject's session's runFold field
        exper.(subj){curSession}.runFold = {estRunFolders, valRunFolders};
    end
end

if options.flattenSessions
    % store the estimation and validation runs for the current run
    % into the subject's session's runFold field
    exper.(subj){1}.runFold = {masterEstRunFolders, masterValRunFolders};
end

%_________________________
% BACK UP COPY IN ARCHIVE
save(experimentFilename, 'exper');