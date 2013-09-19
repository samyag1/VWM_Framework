function mrifPostRun(xpmt, options)

% if this is the estimation mode and a single lambda value is to be used,
% then load all the beta weights and determine which value of lambda gives
% the best prediction
if options.useSingleLambda
    switch(options.mode)
        case {'est'}
            % use a helper to choose the best lambda and write out the
            % corresponding weights back to their original files
            selectBestLambda(xpmt.estDir, xpmt, options);
        case {'snr'}
            % use a helper to choose the best lambda and write out the
            % corresponding weights back to their original files
            selectBestLambda(xpmt.snrDir, xpmt, options);
    end
end