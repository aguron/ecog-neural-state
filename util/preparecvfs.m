function cvfs = preparecvfs(labels, nFolds)
%PREPARECVFS Prepares a cell array of indicator vectors for the
%   cross-validation folds while ensuring that the proportion of trial
%   types is similar across the folds
%
% INPUTS:
%
% labels      - numeric vector of trial type indices
%
% nFolds      - number of cross-validation folds
%
% OUTPUTS:
%
% cvfs        - a cell array of indicator vectors for the cross-validation
%               folds. Vector entry is true in a cross-validation fold when
%               it belongs to the training set, and false otherwise
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if ~isnumeric(labels) || ~isvector(labels)
    error('labels must be a numeric vector');
  end

  nLabels               = numel(labels);

  categories            = unique(labels);
  nCategories           = numel(categories);
  sCategories           = uhistc(labels); % size of each category

  cvfs                  = cell(nFolds,1);
  for cvf=1:nFolds
    cvfs{cvf}           = false(1,nLabels);
    
    for c=1:nCategories
      categoryInd       = find(ismember(labels,categories(c)));
      fdiv              = floor(linspace(1, sCategories(c)+1, nFolds+1));
      
      categoryIndTest   = categoryInd(fdiv(cvf):fdiv(cvf+1)-1);
      categoryIndTrain  = setdiff(categoryInd,categoryIndTest);
      cvfs{cvf}(categoryIndTrain)...
                        = true;
    end % for c=1:nCategories
    if all(cvfs{cvf})
      error('Unable to partition into valid cross-validation folds');
    end % if all(cvfs{cvf})
  end % for cvf=1:nFolds
end

