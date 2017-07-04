function cvfs = preparecvfs(labels, nFolds)
%PREPARECVFS Prepare indices cross-validation folds
%

  nLabels               = length(labels);

  categories            = unique(labels);
  nCategories           = length(categories);
  sCategories           = histc(labels,categories); % size of each category

  cvfs                  = cell(nFolds,1);
  for cvf=1:nFolds
    cvfs{cvf}           = false(1,nLabels);
    
    for c=1:nCategories
      if iscell(labels)
        categoryInd     = find(ismember(labels,categories{c}));
      else
        categoryInd     = find(ismember(labels,categories(c)));
      end % if iscell(labels)
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

