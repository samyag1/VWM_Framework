function [X,frames] = featuresToDesignMatrix(paradigm,features,presHz,TR,nDummy,stimWindow,ignoreIdx,binsFIR, stimAsFeatures)
%  [X,X0] = featuresToParadigm(paradigm,features,presHz,TR,nDummy,stimWindow,ignore)
%-----------------------------------------------------------------------------------
%INPUT:
% <frames>:    - an nEvents x 1 vector of stimulus frame indices
%
% <features>:  - an nStimTot x nFeatures matrix of feature values. This is model-
%                specific, where nStimTot is the total number of stimuli available
%                (across all possible paradigms).
%
% <presHz>:    - the presentation rate (i.e. the inverse duration each frame
%                in <frames>)
%
% <TR>:        - the measurement frequency (repetition time).
%
% <nDummy>     - the number of dummy TRs at the beginning and end of stimulus
%                presentation. Default = [0 0].
%
% <stimWindow> - the size of the stimulus window over which to collapse close-
%                occurring stimuli of the same type into a single event. If 
%                this number is negative, all values of the a stimulus within
%                the window are zeroed out (i.e. resulting in a single event
%                of one frame at the first onset of the stimulus).
%
% <binsFIR>    - A vector containing the bin offsets to model, starting with 
%                0 as the TR the stim appeared in. The matrix returned will have
%                (features x binsFIR) columns and thus each regressor will
%                have a copy at each binsFIR timepoint
%
% OUTPUT:
% <X>:         - an nTR x nFeatures FIR design matrix
%-------------------------------------------------------------------------------------
% Note: assumes that frames includes any dummy scans at the beginning/end of
% the experiment.
%-------------------------------------------------------------------------------------

if isstr(paradigm)
	if exist(paradigm,'file')
		p = load(paradigm);
		if isstruct(p)
			f = fields(p);
			p = p.(f{1});
		end
	end
end
frames = p.frames;

if notDefined('nDummy')
	nDummy = [0 0];
end

if notDefined('ignoreIdx')
	ignoreIdx = [0];
end

if notDefined('stimWindow')
	stimWindow = [];
end
if notDefined('binsFIR')
    binsFIR = [0];
end
if notDefined('stimAsFeatures')
    stimAsFeatures = false;
end
if notDefined('excludeValFeatures')
    excludeValFeatures = [];
end


frames0 = frames;
frames(find(frames==ignoreIdx)) = 0;
ignoreIdx = 0;

% PARSE OUT THE STIMULUS INDICES
idx = frames(find(frames));
idx = idx(find(diff([0;idx])'));
uniqueStim = unique(idx);
uniqueStim = sort(uniqueStim);
uniqueStimCount = numel(uniqueStim);
stims = idx(1:uniqueStimCount);


% IF THERE ARE DUMMY SCANS WHERE THERE IS NO
% DATA COLLECTION, REMOVE THEM
frames(1:nDummy(1)*TR*presHz) = [];
% Don't remove the dummy scans at the end since they were added in to allow
% for the HRF to settle down for the last stimuli shown
%frames(numel(frames):-1:numel(frames)+1-nDummy(2)*TR*presHz) = [];

if stimAsFeatures
    nFeatures = uniqueStimCount;
else
    nFeatures = size(features,2);
end

%  stims = unique(frames(frames~=ignoreIdx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO - add support for nuissance regressors

%I.E. TREAT AS SINGLE STIMULUS:
%....|.|.|.... = ....|||||.....
if stimWindow 	
	cnt = 1;
	while 1
		if frames(cnt) ~= ignoreIdx
			if stimWindow > 0
				frames(cnt:cnt+abs(stimWindow)) = frames(cnt);
			else
				frames(cnt+1:cnt+abs(stimWindow)) = 0;
			end
			cnt = cnt + max(1,abs(stimWindow))+1;
		else
			cnt = cnt+1;
		end
		if cnt >=numel(frames)
			break
		end
	end
end

%nEvents = numel(stims);
%nRepeats = sum(ismember(frames,stims))/nEvents; % THIS ASSUMES EQUAL # OF REPEATS PER STIM

% instead of calculating the features for all the frames and downsampling
% the resulting matrix, downsample the frames to create a smaller matrix.
% If there were some sort of interpolation goiing on in the downsampling of
% the matrix, then order would matter, but there's not, so this is more
% efficient and makes dealing with the FIR binning easier
frames = frames(1:presHz*TR:end,:);

% add all the features to the design matrix, using an FIR model that has
% one set of all features per bin specified in the binsFIR vector. This
% vector is a list of offsets (bins) from the current TR that the will have
% the features of the current stimulus added to
matWidth = nFeatures * numel(binsFIR);
nFrames = numel(frames);
%X0 = zeros(nFrames,matWidth);
repCounts = ones(uniqueStimCount,1);
X = zeros(nFrames,matWidth);
for iF = 1:nFrames
	if frames(iF) > 0
        curStimID = frames(iF);
        curUniqueStimID = find(uniqueStim == curStimID);
        curRep = repCounts(curUniqueStimID);
        for iBin = 1:numel(binsFIR)
            curBinOffset = binsFIR(iBin);
            volumeIdx = iF + curBinOffset;
            if volumeIdx <= nFrames
                if stimAsFeatures
                    stimFeatIdx = (iBin-1)*nFeatures + curUniqueStimID;
                    X(volumeIdx,stimFeatIdx) = 1;
                else
                    startFeatIdx = (iBin-1)*nFeatures + 1;
                    endFeatIdx = iBin*nFeatures;
                    X(volumeIdx,startFeatIdx:endFeatIdx) = squeeze(features(curStimID,:,curRep));
                end
            end
        end
        
        % update the counter for how many times the unique stim has been
        % encountered so far. This is used to know the rep index into the
        % features matrix.
        repCounts(curUniqueStimID) = repCounts(curUniqueStimID) + 1;
	end
end

% DOWNSAMPLE TO MEASUREMENT RATE
%X = X0(1:presHz*TR:end,:);
