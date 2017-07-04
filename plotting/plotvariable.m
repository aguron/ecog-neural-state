function results = plotvariable(seq, varargin)
%PLOTVARIABLE plots observed or latent variables
%   
%   yDim: number of electrodes for each process
%   xDim: state dimensionality for each process
%   nProcesses: number of Markov processes
%
% INPUTS:
%
% seq             - data structure with extracted trajectories
%                   	trialId                 -- unique trial identifier
%                     trialType               -- index from 1 to the number
%                                                of trial types inclusive
%                     y (# electrodes x T)    -- neural data
%                     T (1 x 1)               -- number of timesteps
%                     state (mixComp)         -- state (mixComp) at each
%                     (nProcesses x T)           time point
%                     x {nProcesses}          -- state (mixComp)
%                       (xDim(p) x T             at each time point
%                       x nVar(p))
%                     p {nProcesses}          -- state (mixComp)
%                       (nVar(p) x T)         probabilities at each
%                                                time point
% OPTIONAL ARGUMENTS:
%
% variable        - field name of variable in seq to be plotted
%                   'y', 'state' (default), 'mixComp', 'x' or 'p'
% selectVar       - selected variable slice(s) (default: 1)
% process         - selected process (default: 1)
% selected       	- selected state for variable 'x' (default: 1)
% plotRef         - 'type' (default), 'variable', or 'trial'
% selectRef       - selected variable slice(s)
%                   (default: determined by function)
% subPlotRef      - 'variable' (default), 'type', or 'trial'
% nCols           - number of subplot columns (default: 4)
% nRows           - number of subplot rows (default: 4)
% styles          - cell array of graph styles
% LineWidth       - width of plot lines (default: 0.05)
% FontSize        - font size for axis information (default: 24)
%
% @ 2015 Akinyinka Omigbodun -- aomigbod@ucsd.edu

  variable                            = 'state';
  selectVar                           = 1;
  process                             = 1;
  state                               = 1;
  plotRef                             = 'type';
  selectRef                           = [];
  subPlotRef                          = 'variable';
  statespace                          = false;
  discreteInfo                        = false;
  nCols                               = 4;
  nRows                               = 4;
  styles                              =...                     	%line
    {{'r','b','g','m','y','c','k'};                             %COLORS
     {'-','-.', '--', ':'}                                      %STYLES
     {'','o','.','+','*','<','>','^','v','x','s','d','p','h'}};	%MARKERS
	LineWidth                           = 0.05;
  MarkerSize                          = 10;
  FontSize                            = 24;
  
  assignopts(who, varargin);
  
  switch(plotRef)
   case 'type'
    if ~ismember(subPlotRef,{'variable'})
     error('plotRef and subPlotRef incompatible');
    end % if ~ismember(subPlotRef,{'variable'})
    if (statespace) && ~ismember(variable, {'y'})
     error('statespace and variable incompatible');
    end % if (statespace) && ~ismember(variable, {'y'})
    if (discreteInfo)
     error('discreteInfo == true and plotRef incompatible');
    end % if (discreteInfo)
   case 'variable'
    if ~ismember(subPlotRef,{'type', 'trial'})
     error('plotRef and subPlotRef incompatible');
    end % if ~ismember(subPlotRef,{'type', 'trial'})
    if ismember(subPlotRef, {'trial'}) && ~ismember(variable, {'y'})
     error('plotRef and subPlotRef incompatible');
    end % if ismember(subPlotRef, {'trial'}) && ~ismember(variable, {'y'})
    if (statespace)
      error('plotRef cannot be ''variable'' if statespace == true');
    end % if (statespace)
    if (discreteInfo)
     error('discreteInfo == true and plotRef incompatible');
    end % if (discreteInfo)
   case 'trial'
    if ~ismember(subPlotRef,{'variable'})
     error('plotRef and subPlotRef incompatible');
    end % if ~ismember(subPlotRef,{'variable'})
    if (statespace) && ~ismember(variable, {'y'})
     error('statespace and variable incompatible');
    end % if (statespace) && ~ismember(variable, {'y'})
   otherwise
    error('Invalid specification of plotRef');
  end % switch(plotRef)
  
  if ~any(ismember({plotRef, subPlotRef},{'trial'})) &&...
     (numunique([seq.T]) > 1)
    error(['Sequences must have the same length ',...
           'for plot averages or probability plots']);
  end

  switch(variable)
   case {'state', 'mixComp', 'y'}
    switch(subPlotRef)
     case {'variable', 'type', 'trial'}
      switch(plotRef)
       case {'type', 'variable', 'trial'}
        T                             = seq(1).T;
        switch(variable)
         case {'state', 'mixComp'}
          nVar                        =...
           numunique(cell2mat(modifycells({seq.(variable)},...
                                          @(x) x(process,:))'));
         case 'y'
          nVar                        = size(seq(1).(variable),1);
          nDiscrete                   = [];
          if (discreteInfo)
           if isfield(seq, 'mixComp')
            discrete                  = 'mixComp';
           elseif isfield(seq, 'state')
            discrete                  = 'state';
           end
           nDiscrete                  =...
            numunique(cell2mat(modifycells({seq.(discrete)},...
                                          @(x) x(process,:))'));
          end % if (discreteInfo)
         otherwise
          error('Invalid variable specification');
        end % switch(variable)

        if isempty(selectRef)
         if ismember(plotRef, {'type'})
          selectRef                   = unique([seq.trialType]);
         elseif ismember(plotRef, {'variable'})
          selectRef                  	= 1:nVar;
         elseif ismember(plotRef, {'trial'})
          selectRef                   = 1:numel(seq);
         end
        end % if isempty(selectRef)

        if ismember(plotRef, {'type'})
         if ~isempty(setdiff(selectRef, unique([seq.trialType])))
          error('Invalid specification of selectRef');
         end
        elseif ismember(plotRef, {'variable'})
         if ~isempty(setdiff(selectRef, 1:nVar))
          error('Invalid specification of selectRef');
         end
        elseif ismember(plotRef, {'trial'})
         if ~isempty(setdiff(selectRef, 1:numel(seq)))
          error('Invalid specification of selectRef');
         end
        end

        if ismember(subPlotRef, {'type'})
         if ~isempty(setdiff(selectVar, unique([seq.trialType])))
          error('Invalid specification of selectVar');
         end
        elseif ismember(subPlotRef, {'variable'})
         if ~isempty(setdiff(selectVar, 1:nVar))
          error('Invalid specification of selectVar');
         end
         if (statespace)
          nSelectVar                  = numel(selectVar);
          if ~ismember(nSelectVar, [2 3])
           error('numel(selectVar) must be 2 or 3 when statespace == true');
          end % if ~ismember(nSelectVar, [2 3])
         end % if (statespace)
        elseif ismember(subPlotRef, {'trial'})
         if ~isempty(setdiff(selectVar, 1:numel(seq)))
          error('Invalid specification of selectVar');
         end
        end
        
        if isempty(selectVar)
         switch(variable)
          case {'state', 'mixComp', 'y'}
           if ismember(plotRef, {'type','trial'})
             selectVar               	= 1:nVar;
           elseif ismember(plotRef, {'variable'})
             selectVar              	= unique([seq.trialType]);
           end
          otherwise
           error('Invalid variable specification');
         end % switch(variable)
        end % if isempty(selectVar)

        sp                            = 0;
        for i=selectRef
         sp                         	= sp + 1;
         if (statespace)
          if (sp == 1)
           show;
          end % if (sp == 1)
         else % if (~statespace)
          p                          	= mod(sp-1, nRows*nCols);
          if (~p)
            show;
          end % if (~p)
          subplot(nRows, nCols, p+1);
         end
         hold on

         switch(variable)
          case {'state', 'mixComp'}
           if ismember(plotRef, {'type'})
            varVal                   	=...
             cell2mat(modifycells({seq([seq.trialType] == i).(variable)},...
                                       @(x) x(process,:))');
            varVal                   	=...
             histc(varVal,1:nVar,1)/size(varVal,1);
           elseif ismember(plotRef, {'variable'})
            varSeq                    =...
             cell2mat(modifycells({seq.(variable)},...
                                  @(x) x(process,:))');
            varSeq                    = (varSeq == i);
            uTrialTypes               = unique([seq.trialType]);
            nTrialTypes               = numunique([seq.trialType]);
            varVal                    = nan(nTrialTypes,T);
            for tt=1:nTrialTypes
             varVal(tt,:)             =...
              sum(varSeq([seq.trialType]==uTrialTypes(tt),:))./sum(varSeq);
            end % for tt=1:nTrialTypes
           elseif ismember(plotRef, {'trial'})
            varVal                    = seq(i).(variable)(process,:);
           end
          case 'y'
           if ismember(plotRef, {'type'})
            varVal                    =...
             [seq([seq.trialType] == i).(variable)];
            varVal                    = varVal(selectVar,:);
            varVal                    =...
             reshape(varVal,size(varVal,1),T,[]);
            varVal                    = mean(varVal,3);
           elseif ismember(plotRef, {'variable'})
            switch(subPlotRef)
             case 'type'
              uTrialTypes           	= unique([seq.trialType]);
              nTrialTypes            	= numunique([seq.trialType]);
              varVal                 	= nan(nTrialTypes,T);
              for tt=1:nTrialTypes
               varSeq                	=...
                [seq([seq.trialType]==uTrialTypes(tt)).(variable)];
               varSeq                	= varSeq(i,:);
               varSeq                	=...
                reshape(varSeq,size(varSeq,1),T,[]);
               varVal(tt,:)          	= mean(varSeq,3);
              end % for tt=1:nTrialTypes
             case 'trial'
              nTrials                 = numel(selectVar);
              varVal                 	= nan(nTrials,T);
              for tr=1:nTrials
               varVal(tr,:)          	= seq(selectVar(tr)).(variable)(i,:);
              end % for tr=1:nTrials
             otherwise
              error('Invalid variable specification');
            end % switch(subPlotRef)
           elseif ismember(plotRef, {'trial'})
             varVal                  	= seq(i).(variable)(selectVar,:);
           end
          otherwise
           error('Invalid variable specification');
         end % switch(variable)
         results.variable             = variable;
         results.plotRef              = plotRef;
         results.subPlotRef           = subPlotRef;
         results.plot(sp).ref         = i;
         results.plot(sp).var         = varVal;
         leg                          = {};
         if (statespace)
          if (discreteInfo)
           k                         	=...
            ((find(i==selectRef)-1)*nDiscrete) + 1;
          else % if (~discreteInfo)
           k                         	= find(i==selectRef);
          end
          
          for s=1:max([nDiscrete 1])
           nlineStyle                	= cellfun(@numel, styles);
           lsIdx                      =...
            dec2dvb(mod(k+s-2, prod(nlineStyle)), rowvec(nlineStyle),...
                    'endianness', 'little');
           lsIdx                     	= lsIdx + 1;
           lineStyle                 	= '';
           for l=1:numel(lsIdx)
            lineStyle                 = [lineStyle, styles{l}{lsIdx(l)}];
           end % for l=1:numel(lsIdx)
           
           plotData                  	= {};
           for l=1:nSelectVar
            if (discreteInfo)
             plotData                	=...
              [plotData, varVal(l,seq(i).(discrete)(process,:) == s)];
            else % if (~discreteInfo)
             plotData                	= [plotData, varVal(l,:)];
            end
           end % for l=1:nSelectVar
           if isempty(plotData)
            continue
           end % if isempty(plotData)
           plotData                  	= [plotData, lineStyle,...
                                                   'LineWidth', LineWidth];
           switch(nSelectVar)
            case 2
             plot(plotData{:},'Markers', MarkerSize);
            case 3
             plot3(plotData{:},'Markers', MarkerSize);
            otherwise
             error(['Invalid specification of ',...
                    'selectVar for statespace == true']);
           end % switch(nSelectVar)
          end % for s=1:max([nDiscrete 1])
         else % if (~statespace)
          for j=selectVar
           k                         	= find(j==selectVar);

           nlineStyle                  = cellfun(@numel, styles);
           lsIdx                       =...
            dec2dvb(mod(k-1, prod(nlineStyle)), rowvec(nlineStyle),...
                    'endianness', 'little');
           if ismember(plotRef, {'type','variable'})
            lsIdx                    	= lsIdx + 1;
           elseif ismember(plotRef, {'trial'})
            lsIdx(1)                 	= lsIdx(1) + 1;
            lsIdx(2)                 	= lsIdx(2) + 1;
            lsIdx(3)                 	= lsIdx(3) + 2;
           end
           lineStyle                  = '';
           for l=1:numel(lsIdx)
            lineStyle                 = [lineStyle, styles{l}{lsIdx(l)}];
           end % for l=1:numel(lsIdx)

           if ismember(plotRef, {'type'})
            plot(1:T, varVal(j==selectVar,:),...
                 lineStyle, 'LineWidth', LineWidth);
            leg                      	=...
             [leg, sprintf('%s %d', variable, j)];
           elseif ismember(plotRef, {'variable'})
            plot(1:T, varVal(j==selectVar,:),...
                 lineStyle, 'LineWidth', LineWidth);
            leg                      	=...
             [leg, sprintf('%s %d', subPlotRef, j)];
           elseif ismember(plotRef, {'trial'})
            switch(variable)
             case {'state', 'mixComp'}
              stem(find(varVal==j), varVal(varVal==j), lineStyle)
             case 'y'
              plot(1:T, varVal(j==selectVar,:),...
                   lineStyle, 'LineWidth', LineWidth);
              leg                     =...
               [leg, sprintf('%s %d', subPlotRef, j)];
             otherwise
              error('Invalid variable specification');
            end % switch(variable)
           end
          end % for j=selectVar
         end

         if (statespace)
          switch(nSelectVar)
           case 2
            xlabel(sprintf('%s %d (arb. unit)', variable, selectVar(1)));
            ylabel(sprintf('%s %d (arb. unit)', variable, selectVar(2)));
           case 3
            xlabel(sprintf('%s %d (arb. unit)', variable, selectVar(1)));
            ylabel(sprintf('%s %d (arb. unit)', variable, selectVar(2)));
            zlabel(sprintf('%s %d (arb. unit)', variable, selectVar(3)));
           otherwise
            error(['Invalid specification of ',...
                   'selectVar for statespace = true']);
          end % switch(numel(selectVar))
          title(plotRef);
          set(gca,'FontSize',FontSize)
         else
          if ismember(plotRef, {'type'})
           if (find(i==selectRef) == numel(selectRef))
            legend(leg)
           elseif (sp == 1)
            switch(variable)
             case 'state'
              ylabel('P(state)');
             case 'mixComp'
              ylabel('P(mixComp)');
             case 'y'
              ylabel('arb. unit');
             otherwise
              error('Invalid variable specification');
            end
           end
           if ismember(variable,{'state','mixComp'})
            ylim([0 1])
           end % if ismember(variable,{'state','mixComp'})
           title(sprintf('%s %d', plotRef, i));
          elseif ismember(plotRef, {'variable'})
           if (find(i==selectRef) == numel(selectRef))
            legend(leg)
           elseif (sp == 1)
            switch(variable)
             case {'state', 'mixComp'}
              ylabel('P(type)');
             case 'y'
              ylabel('arb. unit');
             otherwise
              error('Invalid variable specification');
            end
           end
           if ismember(variable,{'state','mixComp'})
            ylim([0 1])
           end % if ismember(variable,{'state','mixComp'})
           title(sprintf('%s %d', variable, i));
          elseif ismember(plotRef, {'trial'})
           if (sp == 1)
            switch(variable)
             case 'state'
              ylabel('state');
             case 'mixComp'
              ylabel('mixComp');
             case 'y'
              ylabel('arb. unit');
             otherwise
              error('Invalid variable specification');
            end
           end % if (sp == 1)
           if ismember(variable,{'state','mixComp'})
            ylim([0 nVar])
            set(gca,'YTick',1:nVar,'YTicklabel',1:nVar)
           end % if ismember(variable,{'state','mixComp'})
           title(sprintf('%s %d', plotRef, i));
          end
          set(gca,'FontSize',FontSize)
          xlim([0 T])
         end
        end % for i=selectRef
       otherwise
        error('Invalid plot reference specification');
      end % switch(plotRef)
     otherwise
      error('Invalid subplot reference specification');
    end % switch(subPlotRef)
   otherwise
    error('Invalid variable specification');
  end % switch(variable)
end