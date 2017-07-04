function selectedTrials = balanceselectedtrials(labels, nTrials)
%BALANCESELECTEDTRIALS
  
  if (nTrials > numel(labels))
    error('nTrials must be less than or equal to numel(labels)');
  end % if (nTrials >= numel(labels))
  
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

