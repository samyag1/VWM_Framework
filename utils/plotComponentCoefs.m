function plotComponentCoefs(componentCoefs, componentNos, coefsNum, featureNames, subject, modelName, outputFilePath)
visibleState = 'on';
writePlots = false;
if nargin > 6
%if ~isempty(outputFilePath)
    visibleState = 'off';
    writePlots = true;
end

for curComponentNoIdx = 1:numel(componentNos)
    
    % get the current component number
    curComponentNo = componentNos(curComponentNoIdx);
    
    % sort the coefficients of the component to plot
    component = componentCoefs(:,curComponentNo);
    [compSortNeg, compOrdNeg] = sort(component,1,'ascend');
    [compSortPos, compOrdPos] = sort(component,1,'descend');
    
    % get the top N coefficients, and their feature names, for display
    dispCoefsNeg = compSortNeg(1:coefsNum);
    dispFeatureNamesNeg = featureNames(compOrdNeg(1:coefsNum));
    dispCoefsPos = compSortPos(1:coefsNum);
    dispFeatureNamesPos = featureNames(compOrdPos(1:coefsNum));
    
    % display the coeficients
    f = figure('Visible', visibleState);
    bar(dispCoefsPos);
    title(sprintf('%s %s Component %i Positive Coefficients', subject, modelName, curComponentNo));
    set(gca, 'XTick', 1:coefsNum);
    set(gca, 'XTickLabel', dispFeatureNamesPos);
    rotateXLabels(gca, 45);
    if writePlots
        % write out a jpg file of the plot if a path was provided
        posOutputFilename = fullfile(outputFilePath, sprintf('%s_%s_Component%i_PositiveCoefficients',subject, modelName, curComponentNo));
        print('-dpng', posOutputFilename);        
    end
    
    f = figure('Visible', visibleState);
    bar(dispCoefsNeg);
    title(sprintf('%s %s Component %i Negative Coefficients', subject, modelName, curComponentNo));
    set(gca, 'XTick', 1:coefsNum);
    set(gca, 'XTickLabel', dispFeatureNamesNeg);
    rotateXLabels(gca, 45);
    if writePlots
        % write out a jpg file of the plot if a path was provided
        negOutputFilename = fullfile(outputFilePath, sprintf('%s_%s_Component%i_NegitiveCoefficients',subject, modelName, curComponentNo));
        print('-dpng', negOutputFilename);        
    end
end