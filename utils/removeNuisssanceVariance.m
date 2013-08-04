function [ Yresid ] = removeNuisssanceVariance(Y, runPath, nuissanceFilenames )

% read the nuissance regressor files into a design matrix
X = readNuissanceFiles(runPath, nuissanceFilenames);

% if the current run doesn't have any nuissance regressors, then just
% return the Y values passed in
if isempty(X)
    Yresid = Y;
else
    assert(size(Y,1) == size(X,1));
    % do simple OLS and return the residuals
    [betas, ~, ~, Yresid,~] = ols(X,Y',[],true);
end
end

