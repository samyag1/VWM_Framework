function fit = estimateHRF(Y, X, TR, hrfBasis, w0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if notDefined('hrfBasis'),
	hrfBasis = canonicalBasis(TR);
end

% DISPLAY MODEL ESTIMATION
debug = 0;

[hrfLen,nBasis]=size(hrfBasis);
hrfDur = hrfLen * TR;
[nDataPoints, nVox] = size(Y);

% determine the unique combinations of regressors in the desing matrix and
% make a new design matrix that has these groups as the features
stimGroups = unique(X, 'rows');
for curGroup = 1:size(stimGroups,1)
    if sum(stimGroups(curGroup,:)) == 0
        stimGroups(curGroup,:) = [];
        break;
    end
end
nStimGroups = size(stimGroups,1);
XStimGroups = zeros(nDataPoints, nStimGroups);
for curDataPoint = 1:nDataPoints
    
    % only update the new design matrix if the current row has regressor
    % values (i.e. a stimulus was shown), otherwise the row will remain all
    % zeros
    if sum(X(curDataPoint,:)) ~= 0
        
        % determine which stim group the current row belongs to
        curGroup = find(ismember(stimGroups, X(curDataPoint,:), 'rows'));
        XStimGroups(curDataPoint,curGroup) = 1;
    end
end

tol = 1e-8;
nIters = 100;

% INITIALIZE STORAGE
basisParams = zeros(nBasis,nVox);
kernels = zeros(hrfLen,nVox);
amp = zeros(nStimGroups,nVox);
iters = zeros(1,nVox);
scale = iters; ttp = iters; fstat = iters; pval = iters;
R2 = iters; normFactor = iters;
fwhm = zeros(2,nVox);

fprintf('\r%1.2f %% complete',0);
for iV = 1:nVox
	y = Y(:,iV);
	iter = 1;
	failed = 0;
	a = 1/nStimGroups*ones(nStimGroups,1);	% ACTIVATIONS
	aOld = -inf(nStimGroups,1); 						% OLD ACTIVATIONS
	
	if notDefined('w0')
		w = ones(size(hrfBasis,2),1);						% BASIS WEIGHTS
	else
		w = w0(:,iV);
	end
	
	while 1
		% ESTIMATE EVENT AMPLITUDES USING HRF
		hrf = hrfBasis*w;
        XAMP = convolveHRF(XStimGroups, hrf);

		try
%            a = ols(XAMP,Y');
            a = gls(Y(:,iV),XAMP,4);
        catch error
            failed = 1; 
        end

		% CONVERGED?
		% (MEAN SQUARED DIFFERENCE OF ACTIVATIONS, CORRELATION
		% (OF WEIGHTS, MAX ITERATIONS, FAILURE FLAG)
		if iter > 1 && ((sum((a - aOld).^2)/(nStimGroups) <= tol) || ((1 - abs(corr2(a,aOld)))<=tol) || (iter >= nIters) || failed == 1)
			
			[fStat(iV),R2(iV),pVal(iV)]=fStats(y,XAMP,a);
			
			% CALC HRF AND SAMPLE TO TR
			basisParams(:,iV) = w;
			
			% NORMALIZE HRF TO BE UNIT LENGTH, 
			% SCALE RESPONSE AMPLITUDES ACCORDINGLY
			[hrf,flip(iV),normFact] = normFlipHRF(hrf);
			a = a.*normFact;
			
			%---------------------------------
			% LOG INFO
			try
				if failed == 1,error(' ');end
				kernels(:,iV) = hrf;
				normFactors(iV) = normFact;
				amp(:,iV) = a;
				iters(iV) = iter;
				
				% WIDTH AND SCALE (INTEGRAL OF WIDTH), TIME TO PEAK
%				[fwhm(:,iV),scale(:,iV),ttp(iV)] = calcFWHM(hrf,info.TR);
			catch
				kernels(:,iV) = NaN(size(kernels,1),1);
				amp(:,iV) = NaN(nStimGroups,1);
				normFactors(iV) = NaN;
				iters(iV) = iter;
				% HRF SCALE...
				scale(iV) = NaN;
				% TIME TO PEAK...
				ttp(iV) = NaN;
				fwhm(:,iV) = [NaN NaN];
				fStat(iV) = NaN;
				R2(iV) = NaN;
				pVal(iV) = NaN;
				failed = 1;
			end
			%----------------------------------
			if iter >= nIters & ~failed
				fprintf('\r%1.2f %% complete (max iterations)',iV/nVox*100);
			elseif ~failed
				fprintf('\r%1.2f %% complete (converged)     ',iV/nVox*100);
			else
				fprintf('\r%1.2f %% complete (failed)        ',iV/nVox*100);
			end
			break
		end

		aOld = a;

        XHRF = createBasisMatrix(XStimGroups, hrfBasis, a);
		try 
%            ols(XHRF,Y');
            w = gls(Y(:,iV),XHRF,4);
        catch error
            failed = 1; 
        end

		
		% DISPLAY
		if debug
			[f,r2,pv]=fStats(y,X1,a);
			clf(gcf)
			subplot(221);plot(hrfBasis*w);

			subplot(222);bar(a'); colormap hot;
			xlabel('Event'); ylabel('Amplitude');
			title(sprintf('Iteration %d',iter));
			xlim([0,nStimGroups+1]);
			subplot(212);plot(y,'r');
			title(sprintf('R2 = %2.2f; F-value = %1.2f; p-val = %1.5f',r2,f,pv))
			y1 = X1*a;
			y2 = X2*w;
			hold on; plot(y1,'b');
			plot(y2,'g');
			legend({ 'Actual','Fit 1','Fit 2' })
			for iR = 1:nRuns-1
				line([TR*iR,TR*iR],[-max(abs(y)) max(abs(y))],'color','k','Linewidth',2)
			end
			drawnow
		end
		iter = iter+1;
	end % END WHILE
end % END LOOP OVER VOX

fit.hrfLen = hrfLen;
fit.hrfDur = hrfDur;
fit.hrf = kernels; clear kernels;
fit.basis = hrfBasis;
fit.basisParams = basisParams;
fit.scale = scale; clear scale;
fit.amp = amp; clear amp;
fit.ttp = ttp; clear ttp;
fit.fwhm = fwhm;
fit.iters = iters;
fit.fStat = fStat;
fit.R2 = R2;
fit.pVal = pVal;
fit.flip = flip;
fit.normFactor = normFactors;
fprintf('\ndone.\n');


end

