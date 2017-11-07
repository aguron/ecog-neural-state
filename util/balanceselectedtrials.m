function selectedTrials = balanceselectedtrials(labels, nTrials)
%BALANCESELECTEDTRIALS Prepares an indicator vector for selected trials
%   while attempting to equalize the numbers of the trial types in the
%   selection
%
% INPUTS:
%
% labels          - numeric vector of trial type indices
%
% nTrials         - number of trials to be selected
%
% OUTPUTS:
%
% selectedTrials	- an indicator vector for selected trials. Vector entry
%                   is true if a trial is selected, and false otherwise
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if (nTrials > numel(labels))
    error('nTrials must be less than or equal to numel(labels)');
  end

  if (nTrials == numel(labels))
    selectedTrials  = true(size(labels));
    return
  end
  
  selectedTrials    = false(size(labels));
  
  categories        = unique(labels);
  nCategories       = numel(categories);
  sCategories       = uhistc(labels); % size of each category

  while (nTrials < sum(sCategories))
    i               = maxidx(sCategories);
    sCategories(i)  = sCategories(i) - 1;
  end % while (nTrials < sum(sCategories)

  for c=1:nCategories
    if (sCategories(c))
      selectedTrials(find(labels == categories(c),sCategories(c)))...
                    = true;
    end % if (sCategories(c))
  end % for c=1:nCategories
end

