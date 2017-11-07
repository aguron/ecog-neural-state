function [accuracy, idx] = crosstabaccuracy(A)
%CROSSTABACCURACY computes classification accuracy from a square
%   cross tabulation
%
% INPUTS:
%
% A           - square cross tabulation
%
% OUTPUTS:
%
% accuracy   	- classification (based on an optimal assignment ordering)
% idx         - ordering of assignments for the best class separation
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if ~issqmat(A)
    error('Cross tabulation must be a square matrix');
  end % if ~issqmat(A)
  d         = size(A,1);

  p         = perms(1:d);

  nP        = size(p,1);
  w         = nan(1,nP);
  for i=1:nP
    w(i)    = sum(A(sub2ind(size(A),1:d,p(i,:))));
  end % for i=1:nP
  
  [a, b]    = max(w);
  accuracy	= a/sum(A(:));
  idx       = [1:d; p(b,:)]';
end