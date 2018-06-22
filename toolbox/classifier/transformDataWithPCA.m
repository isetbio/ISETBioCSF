function [data, coeff, varianceExplained] = transformDataWithPCA(data,numPCAComponents,STANDARDIZE)
% data = transformDataWithPCA(data,numPCAComponents,[STANDARDIZE])
%
% Projects data along a specified number of principal components. PCA is
% done on the data matrix and the first numPCA components (ordered by
% variance explained) will be used to project the data into a lower
% dimensional space.
%
% By default, the data will be standardized before the PCA calculation.
% 
% Inputs:
%   data             -  A matrix containing data to project into a lower dimension.
%                       Rows represent instances of data and columns are features.
%
%   numPCAComponents -  The number of principal components to project onto.
%
%   STANDARDIZE      -  true/false: determine whether or not to standardize the
%                       data.
%
% 7/7/16  xd  wrote it

%% Set default
if (nargin < 3 || isempty(STANDARDIZE))
    STANDARDIZE = true;
end

%% Get rid of features with no variance
s = std(data,1);
index = find(s ~= 0);
data = data(:,index);
s = s(index);

%% Standardize the data
if (STANDARDIZE)
    m = mean(data,1);
    data = (data - repmat(m,size(data,1),1)) ./ repmat(s,size(data,1),1);
end

%% Do PCA and project data into new vector space
if (numPCAComponents > 0)   
    % Request only the first numPCAComponents
    [coeff,score,latent,tsquared,varianceExplained] = pca(data,'NumComponents',numPCAComponents);
    % varianceExplained is equal to latent/sum(latent) * 100
    % tsquared is a statistic measuring how far each observation is from the ?center? of the entire dataset.
    % usefull for identifying outliers
    size(coeff)
    size(score)
    size(latent)
    size(varianceExplained)
    numPCAComponents
    data = data*coeff;
end

end

