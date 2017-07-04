function [accuracy, idx] = crosstabaccuracy(A)
%CROSSTABACCURACY
%   
  if ~issqmat(A)
    error('Cross tabulation must be a square matrix');
  end % if ~issqmat(A)
  d   = size(A,1);
  
  p   = perms(1:d);
  
  nP  = size(p,1);
  w   = nan(1,nP);
  for i=1:nP
    w(i)	=...
      sum(A(sub2ind(size(A),1:d,p(i,:))));
  end % for i=1:nP
  
  [a, b]	= max(w);
  accuracy	= a/sum(A(:));
  idx       = [1:d; p(b,:)]';
end

