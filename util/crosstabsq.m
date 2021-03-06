function sqTable = crosstabsq(labels,assignments,uniqass)
%CROSSTABSQ computes a square cross tabulation for numeric labels and
%   assignments even when the number of unique labels is not equal
%   to the number of unique assignments
%
% INPUTS:
%
% labels      - data labels
% assignments - data assignments
% uniqass   	- specified unique data assignments
%
% OUTPUTS:
%
% sqTable   	- square cross tabulation
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if (nunique(labels) < nunique(assignments))
    error(['Number of unique labels must be greater',...
           ' than the number of unique assignments']);
  end % if (nunique(labels) < nunique(assignments))

  if (nunique(labels) ~= uniqass)
    error(['Number of unique labels must equal the',...
           ' number of specified unique assignments']);
  end % if (nunique(labels) ~= uniqass)
  
  if ~all(ismember(unique(assignments),uniqass))
    error('Invalid unique assignment specification');
  end % if ~all(ismember(unique(assignments),uniqass))

  ind               = ismember(sort(uniqass),unique(assignments));

  table             = crosstab(labels,assignments);
  sqTable           = zeros(size(table,1),numel(ind));
  j                 = 0;
  for i=1:numel(ind)
    if ind(i)
      j             = j + 1;
      sqTable(:,i)	= table(:,j);
    end % if ind(i)
  end % for i=1:numel(ind)
end

