function [loglik, alpha] = hmmFilter2(initDist, transmat, softev)
% Calculate p(S(t)=i| y(1:t))
% INPUT:
% initDist(i) = p(S(1) = i)
% transmat(i,j) = p(S(t) = j | S(t-1)=i)
% softev(i,t) = p(y(t)| S(t)=i)
%
% OUTPUT
% loglik = log p(y(1:T))
% alpha(i,t)  = p(S(t)=i| y(1:t))
%
% Code adapted from hmmFilter.m in pmtk3.
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu



[K T]                          	= size(softev);
scale                           = zeros(T,1);
AT                              = transmat';
if nargout >= 2
    alpha                       = zeros(K,T);
    [alpha(:,1), scale(1)]      =...
      normalize(initDist(:) .* softev(:,1) + eps);
    for t=2:T
        [alpha(:,t), scale(t)]  =...
          normalize((AT * alpha(:,t-1)) .* softev(:,t) + eps);
    end
else
    % save some memory
    [alpha, scale(1)]           =...
      normalize(initDist(:) .* softev(:,1) + eps);
    for t=2:T
     	[alpha, scale(t)]         =...
        normalize((AT * alpha) .* softev(:,t) + eps);
    end
end
loglik                          = sum(log(scale+eps));

end
