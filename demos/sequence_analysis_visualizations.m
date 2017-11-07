%%
% First rank models according to mutual information and prediction error
if (~synthetic)
 str                        = 'real_data';
else % if (synthetic)
 str                        = 'synthetic_data';
end

% Metrics
metrics                     = {'mInf', 'predError'};

% Statistics computed within or across cross-validation folds
stats                       = {'mean', 'median', 'sum'};

for runIdx=pa_select
 if exist(sprintf('results/sequence_analysis/%s/r.mat', str), 'file')
   break
 end
 for metric=metrics
  for stat=stats
   args.(metric{1})         = {'method', {'gmm', 'hmm', 'mfa', 'hmfa'},...
                               'metric', metric{1},...
                               'stat', stat{1},...
                               'nFolds', 4,...
                              };

   for iF=1:numel(fracTrainData)
    for iB=1:numel(binWidth)
     iD                     = 1; % models directory index
     modelDir               =...
       sprintf(['results/sequence_analysis/%s/%g',...
                '/mat_results/run%03d/binWidth_%g'],...
               str, fracTrainData(iF), runIdx, binWidth(iB));
     r(runIdx,iF,iB,iD).(stat{1}).(metric{1})...
                            = rankmodels(modelDir, args.(metric{1}){:});

     iD                     = 2; % models directory index
     modelDir               =...
       sprintf(['results/sequence_analysis/%s/%g',...
                '/mfa-mm/mat_results/run%03d/binWidth_%g'],...
               str, fracTrainData(iF), runIdx, binWidth(iB));
     r(runIdx,iF,iB,iD).(stat{1}).(metric{1})...
                            = rankmodels(modelDir, args.(metric{1}){:});

     iD                     = 3; % models directory index
     modelDir               =...
       sprintf(['results/sequence_analysis/%s/%g',...
                '/gmm-mm/mat_results/run%03d/binWidth_%g'],...
               str, fracTrainData(iF), runIdx, binWidth(iB));
     r(runIdx,iF,iB,iD).(stat{1}).(metric{1})...
                            = rankmodels(modelDir, args.(metric{1}){:});
    end % for iB=1:numel(binWidth)
   end % for iF=1:numel(fracTrainData)
  end % for stat=stats
 end % for metric=metrics
end % for runIdx=pa_select
if ~exist(sprintf('results/sequence_analysis/%s/r.mat', str), 'file')
  save(sprintf('results/sequence_analysis/%s/r.mat', str), 'r');
end
%%
% Plot specified data within ranking structure, r
load(sprintf('results/sequence_analysis/%s/r.mat', str), 'r');

selectB                     = 1;    % bin width
selectF                     = 1:3;  % fraction of training data

nStates                     = 3;
nStatesMax                  = 4;
xDim                        = 1;

% Generate figure comparing prediction error and mutual information
% metrics in the manuscript (Figure 6)
if (~synthetic)
 pa                        	= 2;
else % if (synthetic)
 pa                        	= 1;
end
for plotGroup=4
 pSpecs                     = [];
 selected                   = [];
 switch(plotGroup)
  case 1
   selected(end+1).method   = 'gmm';
   selected(end).covType    = 'diagonal';
   selected(end).sharedCov  = false;
   selected(end).faType     = nan(1,3);

   selected(end+1).method   = 'gmm';
   selected(end).covType    = 'full';
   selected(end).sharedCov  = false;
   selected(end).faType     = nan(1,3);

   selected(end+1).method   = 'mfa';
   selected(end).covType    = '';
   selected(end).sharedCov  = nan;
   selected(end).faType     = [1 1 0];

   selected(end+1).method   = 'mfa';
   selected(end).covType    = '';
   selected(end).sharedCov  = nan;
   selected(end).faType     = [1 1 1];

   LineStyle                = '-';

   iD                       = 1;
  case {2, 3, 4}
   if ismember(plotGroup, [3 4])
    selected(end+1).method  = 'hmm';
    selected(end).covType   = 'diagonal';
    selected(end).sharedCov = false;
    selected(end).faType    = nan(1,3);

    selected(end+1).method	= 'hmm';
    selected(end).covType   = 'full';
    selected(end).sharedCov = false;
    selected(end).faType    = nan(1,3);
   end % if ismember(plotGroup, [3 4])

   if ismember(plotGroup, [2 4])
    selected(end+1).method	= 'hmfa';
    selected(end).covType 	= '';
    selected(end).sharedCov = nan;
    selected(end).faType  	= [1 1 0];

    selected(end+1).method  = 'hmfa';
    selected(end).covType 	= '';
    selected(end).sharedCov = nan;
    selected(end).faType  	= [1 1 1];
   end % if ismember(plotGroup, [2 4])

   switch(plotGroup)
    case {2, 3}
     LineStyle              = '-';

     iD                     = plotGroup;
    case 4
     LineStyle              = '-';

     iD                     = 1;
   end % switch(plotGroup)
 end % switch(plotGroup)
 LineWidth                  = 3;
 
 if ismember(plotGroup, [1 3 4])
  pSpecs(end+1).LineStyle 	= LineStyle;
  pSpecs(end).Marker      	= 'none';
  pSpecs(end).Color       	= (nStates/nStatesMax)*[1 0 0];   % RED
  pSpecs(end).LineWidth   	= LineWidth;
  if (plotGroup == 1)
   pSpecs(end).legendEntry	= sprintf('GMM_{diag}');
  elseif (plotGroup == 3)
   pSpecs(end).legendEntry	= sprintf('GMM_{diag}-MM');
  else
   pSpecs(end).legendEntry	= sprintf('HMM_{diag}');
  end

  pSpecs(end+1).LineStyle 	= LineStyle;
  pSpecs(end).Marker      	= 'none';
  pSpecs(end).Color       	= (nStates/nStatesMax)*[0 1 0];   % GREEN
  pSpecs(end).LineWidth   	= LineWidth;

  if (plotGroup == 1)
   pSpecs(end).legendEntry	= sprintf('GMM_{full}');
  elseif (plotGroup == 3)
   pSpecs(end).legendEntry	= sprintf('GMM_{full}-MM');
  else
   pSpecs(end).legendEntry	= sprintf('HMM_{full}');
  end
 end % if ismember(plotGroup, [1 3 4])

 if ismember(plotGroup, [1 2 4])
  pSpecs(end+1).LineStyle 	= LineStyle;
  pSpecs(end).Marker      	= 'none';
  pSpecs(end).Color       	= (nStates/nStatesMax)*[0 0 1];   % BLUE
  pSpecs(end).LineWidth   	= LineWidth;
  if (plotGroup == 1)
   pSpecs(end).legendEntry	= sprintf('MFA_{CIPN}');
  elseif (plotGroup == 2)
   pSpecs(end).legendEntry	= sprintf('MFA_{CIPN}-MM');
  else
   pSpecs(end).legendEntry	= sprintf('HMFA_{CIPN}');
  end

  pSpecs(end+1).LineStyle 	= LineStyle;
  pSpecs(end).Marker      	= 'none';
  pSpecs(end).Color       	= (nStates/nStatesMax)*[1 0 1];   % MAGENTA
  pSpecs(end).LineWidth   	= LineWidth;
  if (plotGroup == 1)
   pSpecs(end).legendEntry	= sprintf('MFA_{CVPN}');
  elseif (plotGroup == 2)
   pSpecs(end).legendEntry	= sprintf('MFA_{CVPN}-MM');
  else
   pSpecs(end).legendEntry	= sprintf('HMFA_{CVPN}');
  end
 end % if ismember(plotGroup, [1 2 4])

 if (plotGroup == 2)
  newPlot                   = true;
  makeLeg                   = true;
 elseif (plotGroup == 3)
  newPlot                   = false;
  gcf
  makeLeg                   = true;
 else
  newPlot                   = true;
  makeLeg                   = true;
 end

 for metric=metrics
  args                      = {'newPlot', newPlot,...
                               'makeLeg', makeLeg,...
                               'metric', metric{1},...
                               'nMixComp', nStates,...
                               'nStates', nStates,...
                               'xDim', xDim,...
                               'indepVar', fracTrainData(selectF),...
                               };

  for iB=selectB
   stat                     = 'median';
   for iF=selectF
    p(iF)                   = r(pa,iF,iB,iD).(stat).(metric{1});
   end
   plotmetric(p, selected, pSpecs, args{:})

   switch(metric{1})
    case 'predError'
     figure(3)
    case 'mInf'
     figure(1)
    otherwise
     error('Invalid metric specification');
   end
   set(gca, 'XTick',[0 0.5 1], 'XTickLabel',{0,0.5,1});
   xlim([0 1])
   switch(metric{1})
    case 'predError'
     ylim([3800 4300])
    case 'mInf'
     ylim([0 0.01])
    otherwise
     error('Invalid metric specification');
   end
   set(gca,'FontSize',24)
  end % for iB=selectB
 end % for metric=metrics
end

% Generate figure about the Markov constraint (Figure 7) and about the
% Markov constraint and model covariance (Figure 8)
selectF                         = 3;  % fraction of training data
for model_select={1:7, 3:4}
 if isequal(model_select{1}, 1:7)
  if (~synthetic)
   pa                          	= 5;
  else % if (synthetic)
   pa                           = 1;
  end
 else % if isequal(model_select{1}, 3:4)
  pa                           	= 1;
 end
 for iB=selectB
  for iF=selectF
   for plotGroup=[1 3 2 4]
    bSpecs                     	= [];
    selected                    = [];
    switch(plotGroup)
     case 1
      if ismember(1, model_select{1})
       selected(end+1).method 	= 'gmm';
       selected(end).covType  	= 'diagonal';
       selected(end).sharedCov	= false;
       selected(end).faType   	= nan(1,3);
      end

      if ismember(2, model_select{1})
       selected(end+1).method 	= 'gmm';
       selected(end).covType  	= 'full';
       selected(end).sharedCov	= false;
       selected(end).faType   	= nan(1,3);
      end

      if ismember(3, model_select{1})
       selected(end+1).method 	= 'gmm';
       selected(end).covType  	= 'diagonal';
       selected(end).sharedCov	= true;
       selected(end).faType   	= nan(1,3);
      end

      if ismember(4, model_select{1})
       selected(end+1).method 	= 'gmm';
       selected(end).covType  	= 'full';
       selected(end).sharedCov	= true;
       selected(end).faType   	= nan(1,3);
      end

      if ismember(5, model_select{1})
       selected(end+1).method 	= 'mfa';
       selected(end).covType  	= '';
       selected(end).sharedCov	= nan;
       selected(end).faType   	= [1 1 0];
      end

      if ismember(6, model_select{1})
       selected(end+1).method 	= 'mfa';
       selected(end).covType  	= '';
       selected(end).sharedCov	= nan;
       selected(end).faType   	= [1 1 1];
      end

      if ismember(7, model_select{1})
       selected(end+1).method   = 'mfa';
       selected(end).covType  	= '';
       selected(end).sharedCov	= nan;
       selected(end).faType   	= [1 0 0];
      end

      iD                        = 1;
     case {2, 3, 4}
      if ismember(plotGroup, [3 4])
       if ismember(1, model_select{1})
        selected(end+1).method	= 'hmm';
        selected(end).covType 	= 'diagonal';
        selected(end).sharedCov = false;
        selected(end).faType  	= nan(1,3);
       end

       if ismember(2, model_select{1})
        selected(end+1).method	= 'hmm';
        selected(end).covType 	= 'full';
        selected(end).sharedCov = false;
        selected(end).faType  	= nan(1,3);
       end

       if ismember(3, model_select{1})
        selected(end+1).method	= 'hmm';
        selected(end).covType 	= 'diagonal';
        selected(end).sharedCov = true;
        selected(end).faType  	= nan(1,3);
       end

       if ismember(4, model_select{1})
        selected(end+1).method	= 'hmm';
        selected(end).covType 	= 'full';
        selected(end).sharedCov = true;
        selected(end).faType  	= nan(1,3);
       end
      end % if ismember(plotGroup, [3 4])

      if ismember(plotGroup, [2 4])
       if ismember(5, model_select{1})
        selected(end+1).method	= 'hmfa';
        selected(end).covType 	= '';
        selected(end).sharedCov = nan;
        selected(end).faType  	= [1 1 0];
       end

       if ismember(6, model_select{1})
        selected(end+1).method	= 'hmfa';
        selected(end).covType 	= '';
        selected(end).sharedCov = nan;
        selected(end).faType  	= [1 1 1];
       end

       if ismember(7, model_select{1})
        selected(end+1).method	= 'hmfa';
        selected(end).covType 	= '';
        selected(end).sharedCov = nan;
        selected(end).faType  	= [1 0 0];
       end
      end % if ismember(plotGroup, [2 4])

      switch(plotGroup)
       case {2, 3}
        iD                      = plotGroup;
       case 4
        iD                      = 1;
      end % switch(plotGroup)
    end % switch(plotGroup)

    LineStyle                   = '-';
    LineWidth                   = 3;
    if (plotGroup == 4)
      if ismember(1, model_select{1})
       bSpecs(end+1).FaceColor 	= [1 0 0];	% RED
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'Diagonal';
      end

      if ismember(2, model_select{1})
       bSpecs(end+1).FaceColor 	= [0 1 0];	% GREEN
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'Full';
      end

      if ismember(3, model_select{1})
       bSpecs(end+1).FaceColor 	= [0.5 0 0];	% DARK RED 
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'Diagonal-tied';
      end

      if ismember(4, model_select{1})
       bSpecs(end+1).FaceColor 	= [0 0.5 0];	% DARK GREEN
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'Full-tied';
      end

      if ismember(5, model_select{1})
       bSpecs(end+1).FaceColor 	= [0 0 1];	% BLUE
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'FA-CIPN';
      end

      if ismember(6, model_select{1})
       bSpecs(end+1).FaceColor 	= [1 0 1];	% MAGENTA
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'FA-CVPN';
      end

      if ismember(7, model_select{1})
       bSpecs(end+1).FaceColor 	= [0 1 1];	% CYAN
       bSpecs(end).EdgeColor   	= [0 0 0];
       bSpecs(end).LineStyle   	= LineStyle;
       bSpecs(end).LineWidth   	= LineWidth;
       bSpecs(end).legendEntry 	= 'FA-tied';
      end
    end % if (plotGroup == 4)

    if (plotGroup == 4)
     mode                       = 'genFigs';
    else
     mode                       = 'addData';
    end
    if ismember(plotGroup, [1 2])
     groupSpacing               = 0;
    else
     groupSpacing               = 4;
    end
    if (plotGroup == 1)
     data                       = [];
    end

    bSpecs(1).barWidth          = 1;
    bSpecs(1).groupSpacing      = groupSpacing;
    metric                      = 'mInf';
    stat                        = 'median';
    args                        = {'mode', mode,...
                                   'data', data,...
                                   'nMixComp', nStates,...
                                   'nStates', nStates,...
                                   'xDim', xDim,...
                                   'FontSize', 24,...
                                  };
    data                        =...
    plotmetric2(r(pa,iF,iB,iD).(stat).(metric), selected, bSpecs, args{:});
   end

   if isequal(model_select{1}, 1:7)
    figure(5)
   else % if isequal(model_select{1}, 3:4)
    figure(7)
   end
   switch(metric)
    case 'predError'
     switch(pa)
      case 1
       ylim([1600 1800])
      case 2
       ylim([3800 4200])
      case 3
       ylim([4200 4500])
      case 4
       ylim([2800 3100])
      case 5
       ylim([2600 2800])
      case 6
       ylim([12000 20000])
      otherwise
       error('Invalid patient specification');
     end % switch(pa)
    case 'mInf'
     switch(pa)
      case 1
       ylim([0 0.016])
      case 2
       ylim([0 0.012])
      case 3
       ylim([0 0.02])
      case 4
       ylim([0 0.02])
      case 5
       ylim([0 0.012])
      case 6
       ylim([0 0.01])
      otherwise
       error('Invalid patient specification');
     end % switch(pa)
    otherwise
     error('Invalid metric specification');
   end % switch(metric)
  end % for iF=selectF
 end % for iB=selectB
end
%%