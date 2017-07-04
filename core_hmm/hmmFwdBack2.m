function [gamma, alpha, beta, loglik] = hmmFwdBack2(initDist, transmat, softev)
% Calculate p(S(t)=i | y(1:T))
% INPUT:
% initDist(i) = p(S(1) = i)
% transmat(i,j) = p(S(t) = j | S(t-1)=i)
% softev(i,t) = p(y(t)| S(t)=i)
%
% OUTPUT
% gamma(i,t) = p(S(t)=i | y(1:T))
% alpha(i,t)  = p(S(t)=i| y(1:t))
% beta(i,t) propto p(y(t+1:T) | S(t)=i)
% loglik = log p(y(1:T))

% This file is from pmtk3.googlecode.com


% Matlab Version by Kevin Murphy
% C Version by Guillaume Alain
%PMTKauthor Guillaume Alain
%PMTKmex


[loglik, alpha] = hmmFilter2(initDist, transmat, softev);
if any(abs(sum(alpha,1) - 1) > 10*eps)
  disp('alpha (forward probability) columns must sum to 1');
  keyboard
end % if any(abs(sum(alpha,1) - 1) > 10*eps)

beta            = hmmBackwards(transmat, softev);
betaPos         = all(beta < eps);
if any(betaPos)
  disp('beta (backward probability) columns should sum to nonzero value');
	keyboard
end % if any(betaPos)

                % make each column sum to 1
gamma           = normalize(alpha .* beta, 1);

end

function [beta] = hmmBackwards(transmat, softev)
[K T]           = size(softev);
beta            = zeros(K,T);
beta(:,T)       = ones(K,1);
for t=T-1:-1:1
    beta(:,t)   =...
      normalize(transmat * (beta(:,t+1) .* softev(:,t+1)) + eps);
end

end
