function plotmetric(ranking, selected, pSpecs, varargin)
%PLOTMETRIC plots line graphs comparing latent variable models with respect
%   to a specified metric
%
% INPUTS:
%
% ranking       - data for plotting is extracted from a structure RANKING
%                 (which is the output argument of RANKMODELS) or a
%                 structure array of rankings when the independent plot
%                 variable changes in moving from one structure array entry
%                 to another
% selected     	- structure whose i-th entry (corresponding to the i-th
%                 data point) has fields
%                   method        -- latent variable model ('gmm', 'hmm',
%                                    'mfa', 'hmfa')
%                   covType       -- covariance type: 'full', 'diagonal',
%                                    or not applicable to the model ('')
%                   sharedCov     -- covariance tied (true), untied (false)
%                                    or not applicable to the model (nan)
%                   faType        -- factor analyzers specification
%                                    ([nan nan nan] when not applicable
%                                    to the model)
% pSpecs        - structure array (with the i-th entry corresponding to the
%                 i-th plot line) with fields
%                   LineStyle     -- plot line style
%                   Marker        -- plot line marker
%                   Color         -- plot line color
%                   LineWidth     -- plot line width
%                   legendEntry  	-- plot line legend entry
%
% OPTIONAL ARGUMENTS:
%
% newPlot      	- if true, a new plot figure is generated; otherwise the
%                 line plots are made in an existing figure (default: true)
% makeLeg      	- if true, a plot legend is generated (default: true)
% metric      	- 'mInf' or 'predError' (default: 'mInf')
% nMixComp     	- number of mixture components (default: 3)
% nStates     	- number of states (default: 3)
% xDim          - state dimensionality (default: 3)
% indepVar     	- specified independent variable that is not xDim or
%                 nStates/nMixComp (default: [])
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  newPlot                         = true;
  makeLeg                         = true;
  metric                          = 'mInf';
  nMixComp                        = 3;
  nStates                         = 3;
  xDim                            = 3;
  indepVar                        = [];
  assignopts(who, varargin);

  if ~ismember(metric, {'predError', 'mInf'})
    error('Invalid or unsupported metric specification');
  end
  
  if ~any(ismember({selected.method}, {'gmm', 'hmm', 'mfa', 'hmfa'}))
    error('Invalid or unsupported method specified in input');
  else
    if ~isequal(nMixComp, nStates)
      error('nMixComp and nStates must be equal');
    end
  end

  if (numel(ranking) > 1)
    if (numel(indepVar) ~= numel(ranking))
      error('Invalid independent variable specification');
    end
    indepVarPts                   = indepVar;
  else
    if ((numel(nStates) > 1) && (numel(xDim) > 1))
     error(['Only one of nMixComp and/or nStates and ',...
            'xDim can be the independent variable']);
    end
    if (numel(nStates) > 1)
      indepVarPts                	= nStates;
    else
      indepVarPts               	= xDim;
    end
  end
  nIndepVarPts                    = numel(indepVarPts);

  % Reorganize data
  if (numel(ranking) > 1)
    for iIndepVarPts=1:nIndepVarPts
     models                       = ranking(iIndepVarPts).name;
     for i=1:numel(models)
      s.(models{i})(iIndepVarPts)	= ranking(iIndepVarPts).metric(i);
     end
    end
    models                        = fieldnames(s);
  else
    models                        = ranking.name;
    for i=1:numel(models)
     s.(models{i})                = ranking.metric(i);
    end
  end

  switch(metric)
    case 'predError'
      for i=1:numel(models)
        s.(models{i})(s.(models{i}) == 0)...
                                  = nan;
      end % for i=1:numel(models)
    case 'mInf'
      % do nothing
   otherwise
     error('Invalid or unsupported metric specification');
  end % switch(metric)

  % Plot data
  nModels                         = numel(models);
  for i=1:nModels
   M(i)                           = processmodelname(models{i});
  end % for i=1:nModels
  
  if (newPlot)
   show;
  end
  hold on
  
  for i=1:numel(selected)
   sModels    = true(1,nModels);
   sModels    = sModels & ismember({M.method}, {selected(i).method});
   switch(selected(i).method)
    case 'gmm'
     sModels	= sModels & ismember({M.covType}, {selected(i).covType});
     sModels	= sModels & ismember([M.sharedCov], selected(i).sharedCov);
     sModels 	= sModels & ismember([M.nMixComp], nMixComp);
    case 'hmm'
     sModels	= sModels & ismember({M.covType}, {selected(i).covType});
     sModels	= sModels & ismember([M.sharedCov], selected(i).sharedCov);
     sModels  = sModels & ismember([M.nStates], nStates);
    case 'mfa'
     sModels	= sModels &...
      all(bsxfun(@eq,vertcat(M.faType),selected(i).faType),2)';
     sModels	= sModels & ismember([M.xDim], xDim);
     sModels 	= sModels & ismember([M.nMixComp], nMixComp);
    case 'hmfa'
     sModels	= sModels &...
      all(bsxfun(@eq,vertcat(M.faType),selected(i).faType),2)';
     sModels	= sModels & ismember([M.xDim], xDim);
     sModels  = sModels & ismember([M.nStates], nStates);
    otherwise
     error('Invalid or unsupported method');
   end
   idx                            = find(sModels);
   switch(selected(i).method)
    case {'gmm', 'mfa'}
     idx                          = idx(sortidx([M(idx).nMixComp]));
     if isequal(selected(i).method, 'mfa')
       for nS=unique([M(idx).nMixComp])
        idx([M(idx).nMixComp] == nS)...
                                  =...
         arrayaccess(idx([M(idx).nMixComp] == nS),...
                     ['([',...
                      num2str(sortidx([M(idx([M(idx).nMixComp] == nS)).xDim])),...
                      '])']);
       end
     end
    case {'hmm', 'hmfa'}
     idx                          = idx(sortidx([M(idx).nStates]));
     if isequal(selected(i).method, 'hmfa')
       for nS=unique([M(idx).nStates])
        idx([M(idx).nStates] == nS)...
                                  =...
         arrayaccess(idx([M(idx).nStates] == nS),...
                     ['([',...
                      num2str(sortidx([M(idx([M(idx).nStates] == nS)).xDim])),...
                      '])']);
       end
     end
   end

   if (numel(ranking) > 1)
     plot(indepVarPts, s.(models{idx}),...
          'Color', pSpecs(i).Color,...
          'LineWidth', pSpecs(i).LineWidth,...
          'LineStyle', pSpecs(i).LineStyle)    
   else
     for j=idx
      depVarPts(end+1)            = s.(models{j});
     end % for j=idx
     plot(indepVarPts, depVarPts,...
          'Color', pSpecs(i).Color,...
          'LineWidth', pSpecs(i).LineWidth,...
          'LineStyle', pSpecs(i).LineStyle)    
   end
  end % for i=1:numel(selected)

  if (makeLeg)
    makelegend(pSpecs)
  end
end