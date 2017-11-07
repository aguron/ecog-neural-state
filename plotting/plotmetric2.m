function data = plotmetric2(ranking, selected, bSpecs, varargin)
%PLOTMETRIC2 appends data for plotting or plots a bar graph comparing
%   latent variable models with respect to a specified metric
%
% INPUTS:
%
% ranking       - data for plotting is extracted from RANKING (which is
%                 the output argument of RANKMODELS)
% selected     	- structure whose i-th entry (corresponding to the i-th
%                 data point) has fields
%                   method        -- latent variable model ('gmm', 'hmm',
%                                    'mhmm', 'mfa', 'hmfa', 'mhmfa')
%                   covType       -- covariance type: 'full', 'diagonal',
%                                    or not applicable to the model ('')
%                   sharedCov     -- covariance tied (true), untied (false)
%                                    or not applicable to the model (nan)
%                   faType        -- factor analyzers specification
%                                    ([nan nan nan] when not applicable
%                                    to the model)
% bSpecs        - structure array (with the i-th entry corresponding to the
%                 bar of the i-th data point while barWidth and
%                 groupSpacing are only specified in the 1st entry)
%                   FaceColor     -- bar face color
%                   EdgeColor     -- bar edge color
%                   LineStyle     -- bar edge line style
%                   LineWidth     -- bar edge line width
%                   legendEntry  	-- bar legend entry
%                   barWidth      -- bar width
%                   groupSpacing  -- spacing between bar groups
%
% OUTPUTS:
%
% data          - data to be plotted or appended to
%
% OPTIONAL ARGUMENTS:
%
% mode          - plot ('genFigs') or append to ('addData') data
%                 (default: 'genFigs')
% data          - data to be appended to (default: [])
% nMixComp      - number of mixture components (default: 3)
% nStates      	- number of states (default: 3)
% xDim          - state dimensionality (default: 3)
% FontSize     	- font size of plot axes handle (default: 30)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  mode                = 'genFigs';
  data                = [];
  nMixComp            = 3;
  nStates            	= 3;
  xDim               	= 3;
  FontSize            = 30;
  assignopts(who, varargin);
  
  if ~ismember(mode, {'genFigs', 'addData'})
    error('Invalid function mode');
  end

  models              = ranking.name;
  nModels             = numel(models);
  for i=1:nModels
   M(i)               = processmodelname(models{i});
  end
  
  data                = [data, nan(1,bSpecs(1).groupSpacing)];
  
  for i=1:numel(selected)
   sModels           	= true(1,nModels);
   sModels            = sModels & ismember({M.method},{selected(i).method});
   switch(selected(i).method)
    case 'gmm'
     sModels	= sModels & ismember({M.covType}, {selected(i).covType});
     sModels	= sModels & ismember([M.sharedCov], selected(i).sharedCov);
     sModels 	= sModels & ismember([M.nMixComp], nMixComp);
    case 'hmm'
     sModels	= sModels & ismember({M.covType}, {selected(i).covType});
     sModels	= sModels & ismember([M.sharedCov], selected(i).sharedCov);
     sModels  = sModels & ismember([M.nStates], nStates);
    case 'mhmm'
     sModels	= sModels & ismember({M.covType}, {selected(i).covType});
     sModels	= sModels & ismember([M.sharedCov], selected(i).sharedCov);
     sModels 	= sModels & ismember([M.nMixComp], nMixComp);
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
     sModels 	= sModels & ismember([M.nStates], nStates);
    case 'mhmfa'
     sModels	= sModels &...
      all(bsxfun(@eq,vertcat(M.faType),selected(i).faType),2)';
     sModels	= sModels & ismember([M.xDim], xDim);
     sModels 	= sModels & ismember([M.nMixComp], nMixComp);
     sModels 	= sModels & ismember([M.nStates], nStates);
    otherwise
     error('Invalid or unsupported specification of neural state method');
   end % switch(selected(i).method)

   if ~any(sModels)
    data(1,end+1)     = -inf;
   else % if any(sModels)
    data(1,end+1)     = ranking.metric(sModels);
   end   
  end % for i=1:numel(selected)

  if isequal(mode,'genFigs')
   show;
   b                 	= bar([data; nan(1,numel(data))], bSpecs(1).barWidth);
   xlim([0.5 1.5])
   
   groupOffsets       = isnan(data);
   groupOffsets       = diff(groupOffsets);
   
   if (groupOffsets(find(groupOffsets,1)) > 0)
    groupOffsets      = [-1, groupOffsets];
   else
    groupOffsets      = [0, groupOffsets];
   end
   groupOffsets       = (groupOffsets < 0);
   for j=find(groupOffsets)
    for i=1:numel(selected)
     set(b(j+i-1),...
         'FaceColor',bSpecs(i).FaceColor,...
         'LineStyle',bSpecs(i).LineStyle,...
         'LineWidth',bSpecs(i).LineWidth)     
    end % for i=1:numel(selected)
   end % for j=find(groupOffsets)
   set(gca, 'Fontsize', FontSize, 'XTickLabel',{});

   [bSpecs.barWidth]	= deal(bSpecs(1).barWidth);
   makelegend(bSpecs,'FontSize',FontSize)
  end % if isequal(mode,'genFigs')
end