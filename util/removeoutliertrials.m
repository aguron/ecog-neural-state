function L = removeoutliertrials(L, thr)
%REMOVEOUTLIERTRIALS
%   

  nModels               = size(L,1);

  if (nargin == 1)
    thr                 = nModels - 1;
  end % if (nargin == 1)

  if (median([0 thr nModels]) ~= thr)
    error('Invalid outlier threshold specification')
  end % if (median([0 thr nModels]) ~= thr)

  % Least probable model-sequence pair
  [~, lpp]              = min(L,[],2);
  lpp(isnan(diag(L)))   = [];
  lppCounts             = uhistc(lpp);
  lppModels             = unique(lpp);
  idx                   = find(lppCounts == max(lppCounts));
  
  thr2                  = thr;
  for i=idx'
    if (getargout(2,@min,L(lppModels(i),:)) == lppModels(i))
      lppCounts(i)      = lppCounts(i) - 1;
    end % if (getargout(2,@min,L(lppModels(i),:)) == lppModels(i))
    
    if (lppCounts(i) >= thr)
      L(:,lppModels(i))	= NaN;
      L(lppModels(i),:)	= NaN;
      thr2              = thr2 - 1;
    end % if (lppCounts(i) >= thr)
  end % for i=idx'
  
  if (thr2 == thr)
    return
  else
    L                   = removeoutliertrials(L, thr2);
  end % if (thr2 == thr)
end

