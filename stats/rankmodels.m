function ranking = rankmodels(modelDir, varargin)
%RANKMODELS ranks models (in descending order) in a specified directory
%   according to prediction error, mutual information, or accuracy
%
% INPUTS:
%
% modelDir        - specified directory
%
% OUTPUTS:
%
% ranking         - model ranking (in descending order) with fields
%                     fold      -- structure with fields
%                                    metric   -- 2-D (or 3-D) array of
%                                                metric values separated
%                                                according to
%                                                cross-validation folds
%                                                (and specified channel
%                                                groupings)
%                     metric    -- array of metric values
%                     name      -- list of model names
%                     std       -- array of metric value standard
%                                  deviations (only for stat 'mean')
%
% OPTIONAL ARGUMENTS:
%
% methods         - methods to be included in ranking
% metric          - 'LLtest', 'predError', 'mInf', or 'accuracy'
%                   (default: 'predError')
% stat            - statistic computed within (for prediction error) or
%                   across (for all the metrics) cross-validation folds
%                   ('mean', 'sum', 'median') (default: 'mean')
% nFolds          - number of cross-validation folds (default: 4)
% xDim            - state dimensionality (default: [] (unspecified))
% nMixComp        - number of mixture components
%                   (default: [] (unspecified))
% nStates         - number of states (default: [] (unspecified))
% faType          - factor analyzers specification
%                   (default: [] (unspecified))
% covType         - covariance type (default: '' (unspecified))
% sharedCov       - covariance tied (true), untied (false), or
%                   unspecified ([]). (default: [])
% trialTypeMap   	- cell array of length 2 with a numeric array of trial
%                   type indices as the first entry and a numeric array of
%                   remapped trial type indices as the second entry
%                   (default: [] (unspecified so that there is no mapping))
% mInfLatVar     	- selected latent variable ('mixComp' or 'state') in
%                   computing the mutual information for MHMM and MHMFA
%                   (default: 'mixComp')
% nMInf           - if true, mutual information is normalized by entropy
%                   of trial type indices (default: true)
% predErrorChnAss - cell array of length numel(METHODS) with numeric
%                   arrays of values indicating channel groupings in
%                   prediction error computations for each of the methods
%                   included in the ranking (default: [] (unspecified so
%                   that all channels are in a single group for each
%                   method; a subset of the methods selected for the
%                   ranking can be similarly unspecified))
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  methods                           = {'pca', 'ppca', 'fa', 'gpfa',...
                                       'gmm', 'hmm',...
                                       'mfa', 'hmfa',...
                                       'mhmfa', 'mhmm'};
  metric                            = 'predError';
  stat                              = 'mean';
  nFolds                            = 4;
  xDim                              = [];
  nMixComp                          = [];
  nStates                           = [];
  faType                            = [1 1 1];
  covType                           = '';
  sharedCov                        	= [];
  trialTypeMap                      = [];
  mInfLatVar                        = 'mixComp';
  nMInf                             = true;
  predErrorChnAss                   = [];
  assignopts(who, varargin);

  if ~ismember(metric, {'LLtest', 'predError', 'mInf', 'accuracy'})
    error('Invalid metric specification');
  else
   fprintf('metric: %s\n', metric);
  end
  
  if ~ismember(stat, {'mean', 'sum', 'median'})
    error('Invalid stat specification');
  end
  
  switch(stat)
    case {'mean', 'median'}
      % do nothing
    case 'sum'
      if ismember(metric,{'accuracy'})
        error('Incompatible metric and stat');
      end
  end
  
  if isequal(mInfLatVar, 'mixComp') &&...
     ~any(methods, {'mhmfa', 'mhmm'})
    error(['mixComp variable in MHMFA and MHMM methods is not',...
           'equivalent to the mixComp variable in the other methods']);
  end
  
  if ~isempty(predErrorChnAss) &&...
     ~(numel(predErrorChnAss) ~= numel(methods))
    error('predErrorChnSel and methods must have the same length');
  end

  if ~isdir(modelDir)
    fprintf('ERROR: %s does not exist.  Exiting...\n', modelDir);
    ranking                         = [];
    return
  else % if isdir(modelDir)
    D                               = dir([modelDir '/*.mat']);
  end

  if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
    ranking                         = [];
    return
  end % if isempty(D)

  for i=1:numel(D)
    try
      P                            	= processfilename(D(i).name);

      D(i).method                 	= P.method;
      D(i).nMixComp                	= P.nMixComp;
      D(i).nStates                  = P.nStates;
      D(i).xDim                    	= P.xDim;
      D(i).faType                  	= P.faType;
      D(i).covType                 	= P.covType;
      D(i).sharedCov               	= P.sharedCov;
      D(i).cvf                     	= P.cvf;
      if ~isnan(P.ncvf)
        D(i).ncvf                  	= P.ncvf;
      else % if isnan(P.ncvf)
         if ~isnatnum(nFolds)
           error('nFolds must be a positive integer');
         end
        D(i).ncvf                 	= nFolds;
      end
    catch err
      D(i).method                 	= '';
      D(i).nMixComp                	= nan;
      D(i).nStates                 	= nan;
      D(i).xDim                    	= nan;
      D(i).faType                  	= nan(1,3);
      D(i).covType                 	= '';
      D(i).sharedCov               	= nan;
      D(i).cvf                     	= nan;
      D(i).ncvf                     = nan;

      programcontrol
    end % try
  end

  % Specified methods
  if isempty(cell2mat(methods))
    error('methods should be a cell array with valid entries');
  end % if isempty(cell2mat(methods))
  crit                              = ismember({D.methods}, methods);

  % Specified number of cross validation folds
  crit                              = crit & ([D.ncvf] == nFolds);

  % Specified number of mixture components
  if ~isempty(nMixComp)
    crit                = crit &...
      cellfun(@(x)isequal(x, nMixComp) || all(isnan(x)), {D.nMixComp});
  end % if ~isempty(nMixComp)

  % Specified number of states
  if ~isempty(nStates)
    crit                            = crit &...
      cellfun(@(x)isequal(x, nStates) || all(isnan(x)), {D.nStates});
  end % if ~isempty(nStates)

  % Specified latent dimensionality
  if ~isempty(xDim)
    crit                            = crit &...
      cellfun(@(x)isequal(x, xDim) || all(isnan(x)), {D.xDim});
  end % if ~isempty(xDim)

  % Specified faType
  if ~isempty(faType)
    crit                           	= crit &...
      cellfun(@(x)isequal(x, faType) || all(isnan(class2mat(x))),...
              {D.faType});
  end % if ~isempty(faType)

  % Specified covType
  if ~isempty(covType)
    if any(ismember({'gmm', 'hmm', 'mhmm'}, methods))
      crit                          = crit &...
       (ismember({D.covType},covType) | ismember({D.covType},''));
    end % if any(ismember({'gmm', 'hmm', 'mhmm'}, methods))
  end % if ~isempty(covType)

  % Specified sharedCov
  if ~isempty(sharedCov)
    if any(ismember({'gmm', 'hmm', 'mhmm'}, methods))
     crit                           = crit &...
      cellfun(@(x)isequal(x, sharedCov) || all(isnan(x)), {D.sharedCov});
    end % if any(ismember({'gmm', 'hmm', 'mhmm'}, methods))
  end % if ~isempty(sharedCov)

  % Only continue processing files that satisfy criteria
  D                                 = D(crit);

  nFiles                            = numel(D);
  nameList                          = {};
  metricList                        = [];
  numTrialsList                     = [];
  for i=1:nFiles
    fprintf('Loading %s/%s...\n', modelDir, D(i).name);
    ws                              =...
        load(sprintf('%s/%s', modelDir, D(i).name));

    % Compute metric
    switch(metric)
      case 'LLtest'
        if ~isfield(ws, metric)
          fprintf(['%s/%s does not have metric. ',...
                   'Moving to next file...\n'],...
                  modelDir, D(i).name);
          continue
        end
      case 'predError'
        switch(D(i).method)
          case {'pca', 'ppca', 'fa'}
            if ~isfield(ws.kern(1).seqTest, 'ycs')
              fprintf(['%s/%s does not have metric. ',...
                       'Moving to next file...\n'],...
                      modelDir, D(i).name);
              continue
            end
          case 'gpfa'
            if ~isfield(ws.seqTest, 'ycsOrth01')
              fprintf(['%s/%s does not have metric. ',...
                       'Moving to next file...\n'],...
                      modelDir, D(i).name);
              continue
            end
          case {'gmm', 'hmm', 'mfa', 'hmfa'}
            if ~isfield(ws.seqTest, 'ycs')
              fprintf(['%s/%s does not have metric. ',...
                       'Moving to next file...\n'],...
                      modelDir, D(i).name);
              continue
            end
          otherwise
            error('Invalid method specification');
        end % switch(D(i).method)
      case 'mInf'
        switch(D(i).method)
          case {'pca', 'ppca', 'fa', 'gpfa'}
            fprintf(['%s/%s does not have metric. ',...
                     'Moving to next file...\n'],...
                    modelDir, D(i).name);
          case {'gmm', 'hmm', 'mfa', 'hmfa', 'mhmfa', 'mhmm'}
            % do nothing
          otherwise
            error('Invalid method specification');
        end % switch(D(i).method)
      case 'accuracy'
        switch(D(i).method)
          case {'pca', 'ppca', 'fa', 'gpfa',...
                'gmm', 'hmm', 'mfa', 'hmfa'}
            fprintf(['%s/%s does not have metric. ',...
                     'Moving to next file...\n'],...
                    modelDir, D(i).name);
          case {'mhmfa', 'mhmm'}
            % do nothing
          otherwise
            error('Invalid methods');
        end % switch(D(i).method)
      otherwise
        error('Invalid metric specification');
    end % switch(metric)

    switch(metric)
      case 'LLtest'
        nKernSD                    	= 0;
        D(i).(metric)               = ws.(metric);
        D(i).numTrials              = numel(ws.seqTest);
        clear ws
      case 'predError'
        switch(D(i).method)
          case {'pca', 'ppca', 'fa'}
            nKernSD                	= numel(ws.kern);
            D(i).(metric)          	= nan(nKernSD,1);
            D(i).numTrials          = nan(nKernSD,1);
            for kidx=1:nKernSD
              ws.kern(kidx).seqTest	=...
               segmentByTrial(ws.kern(kidx).seqTest, ws.YtestRaw, 'yOrig');
              D(i).(metric)(kidx)   =...
              eval([stat,...
                    '(cell2mat(cellfun(@(x,y)ndarraysum((x-y).^2),',...
                                       '{ws.kern(kidx).seqTest.ycs},',...
                                       '{ws.kern(kidx).seqTest.yOrig},',...
                                      '''UniformOutput'',false)));']);
              D(i).numTrials(kidx) 	= numel(ws.kern(kidx).seqTest);
            end
          case 'gpfa'
            nKernSD               	= 0;
            D(i).numTrials          = numel(ws.seqTest);
            fn                      = sprintf('ycsOrth%02d', D(i).xDim);
            D(i).(metric)           =...
            eval([stat,...
                  '(cell2mat(cellfun(@(x,y)ndarraysum((x-y).^2),',...
                                     '{ws.seqTest.(''',fn,''')},',...
                                     '{ws.seqTest.y},',...
                                    '''UniformOutput'',false)));']);
          case {'gmm', 'hmm', 'mfa', 'hmfa'}
            nKernSD                	= 0;
            D(i).numTrials          = numel(ws.seqTest);

            if ~isfield(ws, 'seqTestOrig')
             str                    = 'seqTest';
            else % if isfield(ws, 'seqTestOrig')
             str                    = 'seqTestOrig';
            end

            m                       = ismember(methods, D(i).method);
            if isempty(predErrorChnAss) || isempty(predErrorChnAss{m})
              D(i).(metric)        	=...
               eval([stat,...
                     '(cell2mat(cellfun(@(x,y)ndarraysum((x-y).^2),',...
                                        '{ws.seqTest.ycs},',...
                                        '{ws.(''',str,''').y},',...
                                       '''UniformOutput'',false)));']);
            else
              for g=rowvec(unique(predErrorChnAss{m}))
                D(i).(metric)(g)   	=...
                  eval([stat,...
                       '(cell2mat(cellfun(@(x,y)ndarraysum((x-y).^2),',...
                                          'modifycells({ws.seqTest.ycs},',...
                                             '@(x)x(predErrorChnAss{m}==g,:)),',...
                                          'modifycells({ws.(str).y},',...
                                             '@(x)x(predErrorChnAss{m}==g,:)),',...
                                         '''UniformOutput'',false)));']);
              end
            end
          otherwise
            error('Invalid methods for metric');
        end % switch(D(i).method)
      case 'mInf'
        nKernSD                     = 0;
        if ~isfield(ws.seqTest,'trialType')
          error('trial type information must be specified')
        end

        switch(D(i).method)
          case {'gmm', 'mfa'}
            D(i).state              = [ws.seqTest.mixComp]';
          case {'hmm', 'hmfa'}
            D(i).state              = [ws.seqTest.state]';
          case {'mhmfa', 'mhmm'}
            D(i).state              = [ws.seqTest.(mInfLatVar)]';
          otherwise
            error('Invalid methods');
        end % switch(D(i).method)
        if ismember(D(i).method, {'mhmfa', 'mhmm'}) &&...
           isequal(mInfLatVar, 'mixComp')
          D(i).trialType            = [ws.seqTest.trialType]';
        else
          D(i).T                   	= [ws.seqTest.T];
          D(i).trialType           	=...
            cell2mat(cellfun(@(x,y)x*ones(1,y),...
                             num2cell([ws.seqTest.trialType]),...
                             num2cell([ws.seqTest.T]),...
                             'UniformOutput', false))';
        end
        if ~isempty(trialTypeMap)
          D(i).trialType           	= maparray(D(i).trialType,....
                                               trialTypeMap{1},...
                                               trialTypeMap{2});
        end % if ~isempty(trialTypeMap)
        X                           = [D(i).state, D(i).trialType];
        mVar                      	= [false(1,size(D(i).state,2)), true];
        D(i).(metric)             	=...
         minformation(X, mVar, 'inputType', 'variable');
        if (nMInf)
         D(i).(metric)             	=...
          D(i).(metric)/entropy2(X(:,end), 'inputType', 'variable');
        end % if (nMInf)
        D(i).numTrials              = numel(ws.seqTest);
      case 'accuracy'
        nKernSD                     = 0;
        D(i).trialType              = [ws.seqTest.trialType];
        D(i).mixComp                = [ws.seqTest.mixComp];
        if ~isempty(trialTypeMap)
         D(i).trialType             = maparray(D(i).trialType,....
                                               trialTypeMap{1},...
                                               trialTypeMap{2});
        end
       	D(i).crossTab               =...
          crosstab(D(i).trialType, D(i).mixComp);
        if (size(D(i).crossTab,1) ~= D(i).nMixComp)
          error('Cross tabulation is incompatible with nMixComp');
        end
        if ~issqmat(D(i).crossTab)
          fprintf('Make cross tabulation a square matrix\n');
          D(i).crossTab             =...
            crosstabsq(D(i).trialType, D(i).mixComp, 1:D(i).nMixComp);
        end % if ~issqmat(D(i).crossTab)
        D(i).(metric)               = crosstabaccuracy(D(i).crossTab);
        D(i).numTrials              = numel(D(i).mixComp);
      otherwise
        error('Invalid metric specification');
    end % switch(metric)

    % Across cross-validation folds
    for kidx=1:max(nKernSD,1)
      undi                          = find(D(i).name == '_');
      if (nKernSD) && (ws.kern(kidx).kernSD)
        name                        = sprintf('%s_kernSD%g',...
                                              D(i).name(1:undi(end)-1),...
                                              ws.kern(kidx).kernSD);
      else
        name                        = D(i).name(1:undi(end)-1);
      end

      [~, pos]                    	= ismember(name,nameList);
      switch(metric)
        case {'LLtest', 'predError', 'mInf', 'accuracy'}
          if (pos)
            pos2                   	= find(isnan(metricList(:,pos,1)),1);
            if (numel(D(i).(metric)) > 1) && (~nKernSD)
              metricList(pos2,pos,:)= D(i).(metric)(:);
            else
              metricList(pos2,pos)  = D(i).(metric)(kidx);
              numTrialsList(pos)    =...
                numTrialsList(pos) + D(i).numTrials(kidx);
            end
          else % if ~(pos)
            nameList                = [nameList, name];
            numTrialsList           =...
              [numTrialsList, D(i).numTrials(kidx)];

            pos2                    = 1;
            pos                     = numel(nameList);

            if (numel(D(i).(metric)) > 1) && (nKernSD == 0)
              nChnGroups            = numel(D(i).(metric));
              if isempty(metricList)
                metricList          =...
                  nan(max([[D.cvf] nFolds 1]),1,nChnGroups);
              elseif (ndims(metricList) == 2)
                metricList(:,:,2:nChnGroups)...
                                    = nan;
              end
              metricList(:,pos,:)   = nan;
              metricList(pos2,pos,:)= D(i).(metric)(:);
            else
              metricList(:,pos,1)   = nan(max([[D.cvf] nFolds 1]),1);
              metricList(:,pos,:)   = nan;
              metricList(pos2,pos,1)= D(i).(metric)(kidx);
            end
          end
        otherwise
          error('Invalid metric specification');
      end % switch(metric)
    end % for kidx=1:max(nKernSD,1)
  end % for i=1:nFiles
  fprintf('\n');

  ranking.fold.metric               = metricList;
  if isempty(metricList)
   metricStat                       = metricList;
  else % if ~isempty(metricList)
   metricList(:,any(isnan(metricList)))...
                                    = nan;
   eval(['metricStat = ','nan',stat,'(matflat(metricList,[3 1]));']);
  end
  [ranking.metric, ordering]        = sort(metricStat,'descend');
  switch(stat)
   case 'mean'
    metricStd                       =...
     std(ndarrayflatten(metricList,[3 1]));
    ranking.std                     = metricStd(ordering);
   case {'median','sum'}
    % do nothing
  end % switch(stat)
  ranking.fold.metric               = ranking.fold.metric(:,ordering,:);
  ranking.name                      = nameList(ordering);
end