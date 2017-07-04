function [BIC, nParams] = bic(modelType, model, LL, nSamples)
%AIC Akaike Information Criterion
  switch(modelType)
    case 'chmfa'
      nParams   = [];
      BIC       = [];
    case 'phmfa'
      nParams   = [];
      BIC       = [];
    case 'hmfa'
      nParams   = length(nonzero(model.pi)) +...
                  length(nonzero(model.trans)) +...
                  numel(model.d) +...
                  numel(model.C) +...
                  length(nonzero(model.R));
      BIC       = -2*(LL - (1/2)*nParams*log(nSamples));
    case 'mfa'
      nParams   = length(nonzero(model.Pi)) +...
                  numel(model.d) +...
                  numel(model.C) +...
                  length(nonzero(model.R));
      BIC       = -2*(LL - (1/2)*nParams*log(nSamples));
    otherwise
      error('invalid specification of neural trajectory extraction method');
  end % switch(modelType)
end