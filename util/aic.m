function [AIC, nParams] = aic(modelType, model, LL)
%AIC Akaike Information Criterion
  switch(modelType)
    case 'chmfa'
      nParams   = [];
      AIC       = [];
    case 'phmfa'
      nParams   = [];
      AIC       = [];
    case 'hmfa'
      nParams   = numel(nonzero(model.pi)) +...
                  numel(nonzero(model.trans)) +...
                  numel(model.d) +...
                  numel(model.C) +...
                  numel(nonzero(model.R));
      AIC       = -2*(LL - nParams);
    case 'mfa'
      nParams   = numel(nonzero(model.Pi)) +...
                  numel(model.d) +...
                  numel(model.C) +...
                  numel(nonzero(model.R));
      AIC       = -2*(LL - nParams);
    otherwise
     error('invalid specification of neural trajectory extraction method');
  end % switch(modelType)
end

