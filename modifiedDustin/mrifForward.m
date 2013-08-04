function model = mrifForward(Y,X,weights,addBias,excludeValFeatures);
%  model = mrifForward(Y,X,weights,addBias);

if addBias
	X = [ones(size(X,1),1),X];
end

% create the indices that are to be used for prediction. THese will be used
% for cross-validation of the estimation weights, as well as determine
% which weights are stored.
useValFeatures = [1:size(X,2)];
useValFeatures(excludeValFeatures) = [];
assert(numel(useValFeatures) == size(weights,1));

XUse = X(:,useValFeatures);
model.pred = XUse*weights;
model.resid = model.pred - Y;
[~,model.r2] = calcR2(Y,model.pred);
[model.cc, model.cI] = ccMatrix(model.pred,Y,1);
model.X = X;
