function hyptest = comparemetric(results, varargin)
%COMPAREMETRIC compares the prediction error, mutual information or
%   accuracy of two models with the Wilcoxon signed rank test
%
% INPUTS:
%
% results       - structure of two entries (with information about the
%                 model files for the two models) with fields
%                   dir               -- model directory information
%                   runIdx           	-- model directory information
%                                        (please see RUNIDX in NEURALSTATE)
%                   binWidth          -- ECoG window width in seconds
%                   method            -- method for model fitting,
%                                        inference, and/or prediction
%                   nFolds            -- number of cross-validation folds
%                   xDim              -- state dimensionality (nan when
%                                        not applicable to model)
%                   nMixComp         	-- number of mixture components (nan
%                                        when not applicable to model)
%                   nStates           -- number of states (nan when not
%                                        applicable to model)
%                   faType            -- factor analyzers specification
%                                        ([nan nan nan] when not applicable
%                                        to the model)
%                   kernSD            -- Gaussian smoothing kernel width
%                                        in seconds (for compatibility
%                                        with TWOSTAGEENGINE2 from
%                                        ecog-gpfa-functionality)
%                   covType          	-- covariance type: 'full',
%                                        'diagonal', or not applicable
%                                        to the model ('')
%                   sharedCov         -- covariance tied (true), 
%                                        untied (false), or not applicable
%                                        to the model (nan)
%
% OUTPUTS:
%
% hyptest       - structure with fields
%                   model             -- struct array of size 2 (with
%                                        information on the 2 models in the
%                                        comparison) with fields
%                                          name          -- model name
%                                          numTrials     -- number of
%                                                           trials in the
%                                                           cross-validation
%                                                           folds
%                                          nMeas         -- number of
%                                                           independent
%                                                           measurements in
%                                                           the
%                                                           cross-validation
%                                                           folds
%                                   EITHER:         ]    -- all independent
%                                          predError]       metric
%                                       OR:         ]       measurements
%                                          mInf     ]       from the
%                                       OR:         ]       cross-validation
%                                          accuracy ]       folds
%                                   EITHER:         ]    -- metric
%                                          predError]       measurements
%                                       OR:         ]       across the
%                                          mInf     ]Fold   cross-validation
%                                       OR:         ]       folds
%                                          accuracy ]
%                   rslt              -- struct array of size 3 (with
%                                        information on two-tailed,
%                                        right-tailed, and left-tailed
%                                        tests respectively) with fields
%                                          tail          -- 'both', 'right'
%                                                           or 'left'
%                                          p             --] Please see
%                                          h             --] MATLAB function
%                                          stats         --] SIGNRANK
%
% OPTIONAL ARGUMENTS:
%
% metric          - 'predError', 'mInf', or 'accuracy'
%                   (default: 'predError')
% stat            - structure with fields
%                     method          -- 'exact', 'approximate', or []
%                                        (unspecified in use of MATLAB
%                                        function SIGNRANK below). Please
%                                        see SIGNRANK optional argument
%                                        METHOD for more information
%                                        (default: [])
%                     alpha           -- scalar value in the range 0 to 1,
%                                        or [] (unspecified in use of
%                                        MATLAB function SIGNRANK below).
%                                        Please see SIGNRANK optional
%                                        argument ALPHA for more
%                                        information (default: [])
%                     maxNumMeasFold  -- maximum number of independent
%                                        mutual information or accuracy
%                                        measurements for each
%                                        cross-validation fold (default: [])
% trialTypeMap   	- cell array of length 2 with a numeric array of trial
%                   type indices as the first entry and a numeric array of
%                   remapped trial type indices as the second entry
%                   (default: [] (unspecified so that there is no mapping))
% mInfLatVar     	- selected latent variable ('mixComp' or 'state') in
%                   computing the mutual information for MHMM and MHMFA
%                   (default: 'mixComp')
% nMInf           - if true, mutual information is normalized by entropy
%                   of trial type indices (default: true)
% predErrorStat  	- computed statistic ('mean', 'median', or 'sum') of
%                   trial prediction error values in a cross-validation
%                   fold (default: 'median')
% predErrorChnSel - cell array of length 2 with numeric arrays with the
%                   indices of subsets of channels to be selected in
%                   prediction error computations for the 2 models in the
%                   comparison. (default: [] (unspecified so that all
%                   channels are selected for each method))
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  metric              	= 'predError';
  stat                  = struct('method',[],...
                                 'alpha',[],...
                                 'maxNumMeasFold',[]);
  trialTypeMap        	= [];
  mInfLatVar           	= 'mixComp';
  nMInf                	= true;
  predErrorStat         = 'median';
  predErrorChnSel       = [];
  assignopts(who, varargin);

  if ~ismember(metric, {'predError', 'mInf', 'accuracy'})
    error('Invalid metric specification');
  else
   fprintf('metric: %s\n', metric);
  end

  nResults              = numel(results);
  if ~isstruct(results) || (nResults ~= 2)
    error('results must be a struct of size 2');
  end
  
  hyptest.model         = struct(metric,cell(1,nResults),...
                                 'name',cell(1,nResults));

  for j=1:nResults
    method              = {results(j).method};
    nFolds             	= results(j).nFolds;
    xDim                = results(j).xDim;
    nMixComp            = results(j).nMixComp;
    nStates             = results(j).nStates;
    faType              = results(j).faType;
    kernSD              = results(j).kernSD;
    covType             = {results(j).covType};
    sharedCov           = results(j).sharedCov;

    directory           =...
     sprintf('%s/mat_results/run%03d/binWidth_%g',...
             results(j).dir, results(j).runIdx, results(j).binWidth);
    if ~isdir(directory)
      fprintf('ERROR: %s does not exist.  Exiting...\n', directory)
      return
    else % if isdir(directory)
      D                 = dir([directory '/*.mat']);
    end

    if isempty(D)
      fprintf('ERROR: No valid files.  Exiting...\n');
      return
    end % if isempty(D)

    for i=1:numel(D)
     try
       P                = processfilename(D(i).name);

       D(i).method      = P.method;
       D(i).nMixComp    = P.nMixComp;
       D(i).nStates     = P.nStates;
       D(i).xDim        = P.xDim;
       D(i).faType      = P.faType;
       D(i).covType     = P.covType;
       D(i).sharedCov  	= P.sharedCov;
       D(i).cvf         = P.cvf;
       if ~isnan(P.ncvf)
         D(i).ncvf      = P.ncvf;
       else % if isnan(P.ncvf)
         if ~isnatnum(nFolds)
           error('results(%d).nFolds must be a positive integer', j);
         end
         D(i).ncvf      = nFolds;
       end
     catch err
       D(i).method      = '';
       D(i).nMixComp    = nan;
       D(i).nStates     = nan;
       D(i).xDim        = nan;
       D(i).faType      = nan(1,3);
       D(i).covType     = '';
       D(i).sharedCov  	= nan;
       D(i).cvf         = nan;
       D(i).ncvf        = nan;

       programcontrol
     end % try
    end

    % Specified method
    crit                = ismember({D.method}, method);
    
    % Specified number of cross-validation folds
    crit                = crit & ([D.ncvf] == nFolds);

    % Specified number of mixture components
    crit                = crit &...
      cellfun(@(x)isequal(x, nMixComp) || all(isnan(x)), {D.nMixComp});
    
    % Specified number of states
    crit                = crit &...
      cellfun(@(x)isequal(x, nStates) || all(isnan(x)), {D.nStates});

    % Specified latent dimensionality
    crit                = crit &...
      cellfun(@(x)isequal(x, xDim) || all(isnan(x)), {D.xDim});

    % Specified faType
    if any(ismember({'mfa', 'hmfa', 'mhmfa'}, method))
      crit              = crit &...
        cellfun(@(x)isequal(x, faType) || all(isnan(class2mat(x))),...
                {D.faType});
    end

    % Specified covType
    if ~isempty(covType)
      if any(ismember({'gmm', 'hmm', 'mhmm'}, method))
        crit          	= crit &...
          ismember({D.covType},covType);
      end
    end % if ~isempty(covType)

    % Specified sharedCov
    if any(ismember({'gmm', 'hmm', 'mhmm'}, method))
     crit               = crit &...
      cellfun(@(x)isequal(x, sharedCov) || all(isnan(x)), {D.sharedCov});
    end

    % Only continue processing files that satisfy criteria
    D                   = D(crit);

    for i=1:numel(D)
      fprintf('Loading %s/%s...\n', directory, D(i).name);
      ws               	=...
          load(sprintf('%s/%s', directory, D(i).name));

      % Compute metric
      switch(metric)
        case 'predError'
          switch(D(i).method)
            case {'pca', 'ppca', 'fa'}
              if ~isfield(ws.kern(1).seqTest, 'ycs')
                error('%s/%s does not have metric.', directory, D(i).name);
              end
            case 'gpfa'
              if ~isfield(ws.seqTest, 'ycsOrth01')
                error('%s/%s does not have metric.', directory, D(i).name);
              end
            case {'gmm', 'hmm', 'mfa', 'hmfa'}
              if ~isfield(ws.seqTest, 'ycs')
                error('%s/%s does not have metric.', directory, D(i).name);
              end
            otherwise
              error('Invalid method');
          end % switch(D(i).method)
        case 'mInf'
          switch(D(i).method)
            case {'pca', 'ppca', 'fa', 'gpfa'}
              error('%s/%s does not have metric.', directory, D(i).name);
            case {'gmm', 'hmm', 'mfa', 'hmfa', 'mhmfa', 'mhmm'}
              % do nothing
            otherwise
              error('Invalid method');
          end % switch(D(i).method)
        case 'accuracy'
          switch(D(i).method)
            case {'pca', 'ppca', 'fa', 'gpfa', 'gmm', 'hmm', 'mfa', 'hmfa'}
              fprintf(['%s/%s does not have metric. ',...
                       'Moving to next file...\n'],...
                      directory, D(i).name);
            case {'mhmfa', 'mhmm'}
              % do nothing
            otherwise
              error('Invalid method');
          end % switch(D(i).method)
        otherwise
          error('Invalid metric specification');
      end % switch(metric)

      switch(metric)
        case 'predError'
          switch(D(i).method)
            case {'pca', 'ppca', 'fa'}
              kidx      = find([ws.kern.kernSD] == kernSD);
              if isempty(kidx)
                error('Metric does not exist for the specified model')
              end % if isempty(kidx)

              D(i).numTrials...
                        = numel(ws.kern(kidx).seqTest);
              ws.kern(kidx).seqTest...
                        = segmentByTrial(ws.kern(kidx).seqTest,...
                                         ws.YtestRaw, 'yOrig');
              for n=1:D(i).numTrials
                if isempty(predErrorChnSel)
                  D(i).(metric)(n)...
                        = ndarraysum((ws.kern(kidx).seqTest(n).ycs -...
                                      ws.kern(kidx).seqTest(n).yOrig).^2);
                else % if ~isempty(predErrorChnSel)
                  D(i).(metric)(n)...
                        = ndarraysum((ws.kern(kidx).seqTest(n).ycs...
                                       (predErrorChnSel{j}) -...
                                      ws.kern(kidx).seqTest(n).yOrig...
                                       (predErrorChnSel{j})).^2);
                end
              end % for n=1:D(i).numTrials
            case 'gpfa'
              D(i).numTrials...
                        = numel(ws.seqTest);
              fn      	= sprintf('ycsOrth%02d', D(i).xDim);
              for n=1:D(i).numTrials
                if isempty(predErrorChnSel)
                  D(i).(metric)(n)...
                        = ndarraysum((ws.seqTest(n).(fn) -...
                                      ws.seqTest(n).y).^2);
                else % if ~isempty(predErrorChnSel)
                  D(i).(metric)(n)...
                        = ndarraysum((ws.seqTest(n).(fn)...
                                       (predErrorChnSel{j}) -...
                                      ws.seqTest(n).y...
                                       (predErrorChnSel{j})).^2);
                end
              end % for n=1:D(i).numTrials
            case {'gmm', 'hmm', 'mfa', 'hmfa'}
              D(i).numTrials...
                        = numel(ws.seqTest);
              for n=1:D(i).numTrials
                if isempty(predErrorChnSel) || isempty(predErrorChnSel{j})
                  if ~isfield(ws, 'seqTestOrig')
                    D(i).(metric)(n)...
                        = ndarraysum((ws.seqTest(n).ycs -...
                                      ws.seqTest(n).y).^2);
                  else % if isfield(ws, 'seqTestOrig')
                    D(i).(metric)(n)...
                        = ndarraysum((ws.seqTest(n).ycs -...
                                      ws.seqTestOrig(n).y).^2);
                  end
                else
                  if ~isfield(ws, 'seqTestOrig')
                    D(i).(metric)(n)...
                        = ndarraysum((ws.seqTest(n).ycs...
                                       (predErrorChnSel{j},:) -...
                                      ws.seqTest(n).y...
                                       (predErrorChnSel{j},:)).^2);
                  else % if isfield(ws, 'seqTestOrig')
                    D(i).(metric)(n)...
                        = ndarraysum((ws.seqTest(n).ycs...
                                       (predErrorChnSel{j},:) -...
                                      ws.seqTestOrig(n).y...
                                       (predErrorChnSel{j},:)).^2);
                  end
                end
              end % for n=1:D(i).numTrials
            otherwise
              error('Invalid method for metric');
          end % switch(D(i).method)
          D(i).([metric,'Fold'])...
                        = eval([predErrorStat,'(D(i).(metric))']);
          D(i).nMeas    = D(i).numTrials;
        case 'mInf'
          if ~isfield(ws.seqTest,'trialType')
            error('trial type information must be specified')
          end
          switch(D(i).method)
           case {'gmm', 'mfa'}
            D(i).state  = [ws.seqTest.mixComp]';
           case {'hmm', 'hmfa'}
            D(i).state	= [ws.seqTest.state]';
           case {'mhmfa', 'mhmm'}
            D(i).state 	= [ws.seqTest.(mInfLatVar)]';
           otherwise
            error('Invalid method');
          end % switch(D(i).method)
          D(i).T       	= [ws.seqTest.T];
          D(i).trialType= [ws.seqTest.trialType];

          % Determine the smallest trial type set and the cardinality
          D(i).nMeas    = min(uhistc(D(i).trialType));
          D(i).nMeas    = min([stat.maxNumMeasFold, D(i).nMeas]);
          D(i).trialSets= preparecvfs(D(i).trialType, D(i).nMeas);

          if ~isempty(trialTypeMap)
            D(i).trialType...
                        =...
              maparray(D(i).trialType, trialTypeMap{1}, trialTypeMap{2});
          end % if ~isempty(trialTypeMap)

          for n=1:D(i).nMeas
           if ismember(D(i).method, {'mhmfa', 'mhmm'}) &&...
              isequal(mInfLatVar, 'mixComp')
            X           =...
             [D(i).state(~D(i).trialSets{n},:),...
              D(i).trialType(~D(i).trialSets{n})'];
           else
          	X          	=...
             [D(i).state(subvecind([],D(i).T,find(~D(i).trialSets{n})),...
                         :),...
              cell2mat(cellfun(@(x,y)x*ones(1,y),...
                               num2cell(D(i).trialType(~D(i)...
                                                       .trialSets{n})),...
                               num2cell(D(i).T(~D(i).trialSets{n})),...
                               'UniformOutput', false))'];
           end
           mVar       	= [false(1,size(D(i).state,2)), true];
           D(i).(metric)(n)...
                        = minformation(X, mVar, 'inputType','variable');
           if (nMInf)
            D(i).(metric)(n)...
                        = D(i).(metric)(n)/...
                          entropy2(X(:,end), 'inputType', 'variable');
           end % if (nMInf)
          end
          if ismember(D(i).method, {'mhmfa', 'mhmm'}) &&...
             isequal(mInfLatVar, 'mixComp')
           X            = [D(i).state, D(i).trialType'];            
          else
           X            =...
           	[D(i).state, cell2mat(cellfun(@(x,y)x*ones(1,y),...
                                          num2cell(D(i).trialType),...
                                          num2cell(D(i).T),...
                                          'UniformOutput', false))'];
          end
         	D(i).([metric,'Fold'])...
                        = minformation(X, mVar, 'inputType', 'variable');
          if (nMInf)
           D(i).([metric,'Fold'])...
                        = D(i).([metric,'Fold'])/...
                          entropy2(X(:,end), 'inputType', 'variable');
          end % if (nMInf)
          D(i).numTrials= numel(ws.seqTest);
        case 'accuracy'
          if ~isfield(ws.seqTest,'trialType')
            error('trial type information must be specified')
          end
          D(i).mixComp  = [ws.seqTest.mixComp];
          D(i).trialType= [ws.seqTest.trialType];

          % determine the smallest trial type set and the cardinality
          D(i).nMeas    = min(uhistc(D(i).trialType));
          D(i).nMeas    = min([stat.maxNumMeasFold, D(i).nMeas]);
          D(i).trialSets= preparecvfs(D(i).trialType, D(i).nMeas);
          
          if ~isempty(trialTypeMap)
            D(i).trialType...
                        =...
              maparray(D(i).trialType, trialTypeMap{1}, trialTypeMap{2});
          end % if ~isempty(trialTypeMap)

          for n=1:D(i).nMeas
            D(i).crossTab...
                        = crosstab(D(i).trialType(~D(i).trialSets{n}),...
                                   D(i).mixComp(~D(i).trialSets{n}));
            if (size(D(i).crossTab,1) ~= D(i).nMixComp)
              error('Cross tabulation is incompatible with nMixComp');
            end
            if ~issqmat(D(i).crossTab)
              fprintf('Make cross tabulation a square matrix\n');
              D(i).crossTab...
                        =...
                crosstabsq(D(i).trialType(~D(i).trialSets{n}),...
                           D(i).mixComp(~D(i).trialSets{n}),...
                           1:D(i).nMixComp);
            end
            D(i).(metric)(n)...
                        = crosstabaccuracy(D(i).crossTab);
          end
          D(i).crossTab = crosstab(D(i).trialType, D(i).mixComp);
          if (size(D(i).crossTab,1) ~= D(i).nMixComp)
            error('Cross tabulation is incompatible with nMixComp');
          end
          if ~issqmat(D(i).crossTab)
            fprintf('Make cross tabulation a square matrix\n');
            D(i).crossTab...
                        =...
              crosstabsq(D(i).trialType, D(i).mixComp, 1:D(i).nMixComp);
          end % if ~issqmat(D(i).crossTab)
          D(i).([metric,'Fold'])...
                        = crosstabaccuracy(D(i).crossTab);
          D(i).numTrials= numel(D(i).mixComp);
        otherwise
          error('Invalid metric specification');
      end % switch(metric)

      % Across cross-validation folds
      undi              = find(D(i).name == '_');
      if exist('kidx','var') && (ws.kern(kidx).kernSD)
        name            = sprintf('%s_kernSD%g',...
                                  D(i).name(1:undi(end)-1),...
                                  ws.kern(kidx).kernSD);
      else
        name            = D(i).name(1:undi(end)-1);
      end

      if isempty(hyptest.model(j).name)
        hyptest.model(j).name...
                        = name;
      elseif ~isequal(hyptest.model(j).name, name)
        error('Inconsistent model names');
      end
    end % for i=1:numel(D)
    
    hyptest.model(j).(metric)...
                        = [D.(metric)];
    hyptest.model(j).numTrials...
                        = [D.numTrials];
    if ismember(metric, {'predError', 'mInf', 'accuracy'})
      hyptest.model(j).nMeas...
                        = [D.nMeas];
      hyptest.model(j).([metric,'Fold'])...
                        = [D.([metric,'Fold'])];
    end % if ismember(metric, {'predError', 'mInf', 'accuracy'})
  end % for j=1:nResults
  fprintf('\n');

  tail                  = {'both', 'right', 'left'};
  for t=1:numel(tail)
   args                 = {'tail', tail{t}};
   if ~isempty(stat.method)
     args               = [args, 'method', stat.method];
   end % if ~isempty(statMethod)
   if ~isempty(stat.alpha)
     args               = [args, 'alpha', stat.alpha];
   end % if ~isempty(statAlpha)   
   hyptest.rslt(t).tail = tail{t};
   [hyptest.rslt(t).p, hyptest.rslt(t).h, hyptest.rslt(t).stats]...
                        = signrank(hyptest.model(1).(metric)',...
                                   hyptest.model(2).(metric)',...
                                   args{:});
  end % for t=1:numel(tail)
end