function mhmmEngine(seqTrain, seqTest, fname, varargin)
%
% mhmmEngine(seqTrain, seqTest, fname, ...) 
%
% Model fitting and inference with MHMM
%
%   yDim: number of electrodes or channels
%
% INPUTS:
%
% seqTrain    - training data structure, whose n-th entry (corresponding
%               to the n-th experimental trial) has fields
%                 trialId           -- unique trial identifier
%                 trialType (1 x 1)	-- trial type index (Optional)
%                 fs (1 x 1)       	-- sampling frequency of ECoG data
%                 T (1 x 1)         -- number of timesteps
%                 y (yDim x T)      -- neural data
% seqTest     - test data structure (same format as seqTrain)
% fname       - model filename
%
% OPTIONAL ARGUMENTS:
%
% nStates    	- number of Markov (HMM) states (default: 3)
% nMixComp   	- number of MHMM mixture components (default: 3)
%
% binWidth  	- ECoG window width in seconds (default: 0.2)
% kernSD      - Gaussian smoothing kernel width in seconds. 0 corresponds
%               to no smoothing (default: 0)
%
% d {nMixComp x 1} (yDim x nStates)               - cell array of component
%                                                   HMM observation means
% R {nMixComp x 1} (yDim x yDim x nStates (or 1))	- cell array of component
%                                                   HMM covariances
%
% pi0 {nMixComp x 1} (1 x nStates)               	- cell array of component
%                                                   HMM initial start
%                                                   probabilities
% piPrior {nMixComp x 1} (1 x nStates)           	- cell array of pseudo
%                                                   counts for component
%                                                   HMM start probabilities
% trans0 {nMixComp x 1} (nStates x nStates)     	- cell array of component
%                                                   HMM initial transition
%                                                   probabilities
% transPrior {nMixComp x 1} (nStates x nStates) 	- cell array of pseudo
%                                                   counts for component
%                                                   HMM transition
%                                                   probabilities
% Pi (1 x nMixComp)                               - component HMM priors
%
% stateGuess                                      - cell array of initial
%                                                   state guesses for time
%                                                 	points of trials for
%                                                  	training data
%
% mixCompGuess                                    - vector of initial
%                                                   component HMM guesses
%                                                   for training trials.
%                                                   Use an index of 0 if
%                                                   there is no guess for a
%                                                   trial; an index between
%                                                   1 and nMixComp
%                                                   (inclusive), if the
%                                                   trial is to be used in
%                                                   initialization; and an
%                                                   index between -nMixComp
%                                                   and -1 (inclusive), if
%                                                  	the trial is not to be
%                                                   used in model fitting
%	removeOutliers                                 	- if true, outlier trials
%                                                   identified during the
%                                                   initialization stage
%                                                   are eliminated from
%                                                   model fitting
%                                                   (default: true)
%
% init                                            - function handle for
%                                                   initialization:
%                                                   @init_mhmm or
%                                                   @init_mhmfa
%                                                   (default: @init_mhmm)
%
% learning                                        - indicates whether model
%                                                   fitting should be
%                                                   carried out
%                                                   (default: true)
% inference                                       - indicates whether
%                                                   latent variables should
%                                                   be inferred for seqTest
%                                                   (default: true)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu
  
  nStates                                   = 3;
  nMixComp                                  = 3;
  
  binWidth                                  = 0.2;
  kernSD                                    = 0;

  d                                         = [];
  R                                         = [];
  pi0                                       = [];
  piPrior                                   = [];
  trans0                                    = [];
  transPrior                                = [];
  Pi                                        = [];

  stateGuess                                = [];
  
  mixCompGuess                              = [];
  removeOutliers                            = true;
  
  init                                      = @init_mhmm;
  
  learning                                  = true;
  inference                                 = true;
  
  extraOpts                                 = assignopts(who, varargin);

  if (learning)
    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and MHMM...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig                          = seqTrain;
      for n=1:numel(seqTrain)
        seqTrain(n).y                       =...
          smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:numel(seqTrain)

      % Test data
      seqTestOrig                           = seqTest;
      for n=1:numel(seqTest)
        seqTest(n).y                        =...
          smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:numel(seqTest)
    end % if (kernSD)

    % Resegment trials for training
    [seqTrainCut, resegTrlGuess]           	=...
      resegmenttrials(seqTrain,...
                      'method', 'mhmm',...
                      'stateGuess', stateGuess,...
                      'mixCompGuess', mixCompGuess,...
                      extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      [seqTrainCut, resegTrlGuess]         	=...
        resegmenttrials(seqTrain,...
                        'method', 'mhmm',...
                        'stateGuess', stateGuess,...
                        'mixCompGuess', mixCompGuess,...
                        'segLength', Inf);
    end

    yAll                                    = [seqTrainCut.y];
    if ~isempty(yAll)
      yDim                                  = size(yAll, 1);
    else % if isempty(yAll)
      yDim                                  = size([seqTest.y], 1);
    end

    initStatus                              =...
      ~[isempty(d) isempty(R),...
        isempty(pi0) isempty(piPrior),...
        isempty(trans0) isempty(transPrior),...
        isempty(Pi)];

    if all(~initStatus) && ~isempty(seqTrainCut)
      if removeOutliers
        argsOut                             = cell(3,1);
      else % if ~removeOutliers
        argsOut                             = cell(2,1);
      end
      [argsOut{:}]                          =...
        init(seqTrainCut, fname,...
             'nStates', nStates, 'nMixComp', nMixComp,...
             'stateGuess', resegTrlGuess.state,...
             'mixCompGuess', resegTrlGuess.mixComp,...
             extraOpts{:});
      startParams                           = argsOut{1};
      if isequal(init,@init_mhmm)
        % do nothing
      else % if isequal(init,@init_mhmfa)
        for k=1:nMixComp
          if (size(startParams(k).R, 3) > 1)
            startParams(k).R               	=...
              apwfun(@(x,y)y*y' + x, startParams(k).R, startParams(k).C);
          else % if (size(startParams(k).R, 3) == 1)
            startParams(k).R               	=...
              bsxfun(@plus,...
                     apwfun(@(y)y*y', startParams(k).C),...
                     startParams(k).R);
          end

          CovType                           = [];
          SharedCov                         = [];
          assignopts(who, varargin);
          if isequal(CovType, 'full')
            % do nothing
          elseif isequal(CovType, 'diagonal')
            startParams(k).R              	=...
              nddiag(nddiag(startParams(k).R));
          end
        end % for k=1:nMixComp
        startParams(1).covType             	= CovType;
        startParams(1).sharedCov           	= SharedCov;
        startParams                        	=...
          rmfield(startParams,{'C','faType'});
      end
      resegTrlGuess.mixComp               	= argsOut{2};
      if removeOutliers
        outliers                            = argsOut{3};
      end % if removeOutliers

      if (numunique(resegTrlGuess.mixComp(abs(resegTrlGuess.mixComp)>0))...
          ~= nMixComp)
        error('Invalid mixCompGuess specification');
      end
    elseif ~any(~initStatus)
      fprintf('Initializing parameters with input...\n');
      % ==================================
      % Initialize state model parameters
      % ==================================
      startParams(1).nStates               	= nStates;
      
      CovType                               = [];
      SharedCov                             = [];
      assignopts(who, varargin);
      startParams(1).covType               	= CovType;
      startParams(1).sharedCov              = SharedCov;
      
      startParams(1).nMixComp              	= nMixComp;

      startParams(nMixComp).Pi              = num2cell(Pi);
      [startParams.Pi]                      =...
        deal(startParams(nMixComp).Pi{:});
      startParams(nMixComp).pi              = class2cell(pi0);
      [startParams.pi]                      =...
        deal(startParams(nMixComp).pi{:});
      startParams(nMixComp).piPrior        	= class2cell(piPrior);
      [startParams.piPrior]                	=...
        deal(startParams(nMixComp).piPrior{:});
      startParams(nMixComp).trans           = class2cell(trans0);
      [startParams.trans]                   =...
        deal(startParams(nMixComp).trans{:});
      startParams(nMixComp).transPrior     	= class2cell(transPrior);
      [startParams.transPrior]             	=...
        deal(startParams(nMixComp).transPrior{:});

      % ========================================
      % Initialize observation model parameters
      % ========================================
      startParams(nMixComp).d               = class2cell(d);
      [startParams.d]                       =...
        deal(startParams(nMixComp).d{:});
      startParams(nMixComp).R               = class2cell(R);
      [startParams.R]                       =...
        deal(startParams(nMixComp).R{:});
    else
      error(['d, R, Pi, pi, piPrior, trans and transPrior must',...
             ' all be specified if there is no data for training.']);
    end

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting MHMM model...\n');

      [estParams, seqTrainCut, LL, iterTime]=...
        em_mhmm(startParams, seqTrainCut,...
               	'outliers', outliers, extraOpts{:});

      % Inference for original, unsegmented trials using learned parameters
      [seqTrain, ess, LLorig]               =...
        exactInferenceWithLL_mhmm(seqTrain, estParams,...
                                 	'getSeq', true,...
                                  'getLL', true,...
                                  extraOpts{:});
    else % if isempty(seqTrainCut)
      estParams                             = startParams;
    end
  else % if (~learning)
    estParams                               = loadvars(fname, 'estParams');
    if (inference)
      if ~isempty(seqTrain)
        % Inference for original, unsegmented trials using learned
        % parameters
        [seqTrain, ess, LLorig]            	=...
          exactInferenceWithLL_mhmm(seqTrain, estParams,...
                                    'getSeq', true,...
                                    'getLL', true,...
                                    extraOpts{:});
      end % if ~isempty(seqTrain)
    end % if (inference)
  end

  if ~isempty(seqTest) % check if there are any test trials
    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      [seqTest, ~, LLtest]                  =...
        exactInferenceWithLL_mhmm(seqTest, estParams,...
                                  'getSeq', true,...
                                  'getLL', true,...
                                  extraOpts{:});
    end % if (inference)
  end % if ~isempty(seqTest)

  % =============
  % Save results
  % =============
  if (learning)
    vars                                  	= who;
    fprintf('Saving %s...\n', fname);
    save(fname, vars{~ismember(vars, {'yAll'})});
  elseif (inference)
    fprintf('Saving inferred latent variables in %s...\n', fname);
    if ~isempty(seqTrain)
      save(fname, 'seqTrain', 'ess', 'LLorig', '-append');
    end % if ~isempty(seqTrain)
    if ~isempty(seqTest)
      save(fname, 'seqTest', 'LLtest', '-append');
    end % if ~isempty(seqTest)
    save(fname, 'inference', '-append');
  end
end