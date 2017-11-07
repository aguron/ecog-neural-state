function hmmEngine(seqTrain, seqTest, fname, varargin)
%
% hmmEngine(seqTrain, seqTest, fname, ...) 
%
% Model fitting and inference with HMM
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
% nStates                           - number of Markov states (default: 3)
% CovType                           - HMM covariance type: 'full' or
%                                     'diagonal' (default: 'diagonal')
% SharedCov                         - HMM covariance tied (true) or
%                                     untied (false) (default: false)
%
% binWidth                          - ECoG window width in seconds
%                                     (default: 0.2)
% kernSD                            - Gaussian smoothing kernel width in
%                                    	seconds. 0 corresponds to no
%                                     smoothing (default: 0)
%
% Replicates                       ]  parameters for
% Regularize                       ]- MATLAB function
% Options                          ]  GMDISTRIBUTION.FIT
%
% d (yDim x nStates)                - observation means
% R (yDim x yDim x nStates (or 1))	- covariance(s)
%
% pi0 (1 x nStates)               	- HMM initial start probabilities
% piPrior (1 x nStates)           	- pseudo counts for HMM start 
%                                     probabilities
% trans0 (nStates x nStates)        - HMM initial transition probabilities
% transPrior (nStates x nStates)    - pseudo counts for HMM transition
%                                     probabilities
%
% stateGuess                      	- cell array of initial state guesses
%                                     for time points of trials for
%                                     training data
%
% learning                          - indicates whether model fitting
%                                     should be carried out (default: true)
% inference                       	- indicates whether latent variables
%                                     should be inferred for seqTest
%                                     (default: true)
% prediction                        - specifies whether
%                                     leave-one-channel-out prediction
%                                    	should be carried out
%                                     (default: false)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  nStates                                   = 3;
  CovType                                   = 'diagonal';
  SharedCov                                 = false;

  binWidth                                  = 0.2;
  kernSD                                    = 0;

  Replicates                                = 1;
  Regularize                                = 0;
  Options                                   = [];
  
  d                                         = [];
  R                                         = [];

  pi0                                       = [];
  piPrior                                   = [];
  trans0                                    = [];
  transPrior                                = [];

  stateGuess                                = [];

  learning                                  = true;
  inference                                 = true;
  prediction                                = false;

  extraOpts                                 = assignopts(who, varargin);
  
  if (learning)
    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and HMM...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig                          = seqTrain;
      for n=1:length(seqTrain)
        seqTrain(n).y                       =...
          smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:length(seqTrain)

      % Test data
      seqTestOrig                           = seqTest;
      for n=1:length(seqTest)
        seqTest(n).y                        =...
          smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:length(seqTest)
    end % if (kernSD)

    % Resegment trials for training
    [seqTrainCut, resegTrlGuess]           	=...
      resegmenttrials(seqTrain,...
                      'method', 'hmm',...
                      'stateGuess', stateGuess,...
                      extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      [seqTrainCut, resegTrlGuess]         	=...
        resegmenttrials(seqTrain,...
                        'method', 'hmm',...
                        'stateGuess', stateGuess,...
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
        isempty(pi0) isempty(piPrior) isempty(trans0) isempty(transPrior)];

    if all(~initStatus) && ~isempty(seqTrainCut)
      fprintf('Initializing parameters with GMM...\n');
      args                                  =...
        {'Replicates', Replicates, 'CovType', CovType,...
         'SharedCov', SharedCov, 'Regularize', Regularize,...
         'Options', Options};
      if ~isempty(resegTrlGuess.state)
        args                                =...
          [args, 'Start', cell2mat(resegTrlGuess.state)];
      end % if ~isempty(resegTrlGuess.state)
      obj                                   =...
        gmdistribution.fit(yAll', nStates, args{:});

     	temp                                  = cluster(obj, yAll');
     	temp                                  =...
        segmentByTrial(seqTrainCut, temp', 'mixComp');
      mixComp                               = cell(1,numel(temp));
      [mixComp{:}]                          = deal(temp.mixComp);

      startParams                           =...
        hmmestimate2(mixComp, mixComp, nStates, nStates,...
                     'Pseudostarts', ones(1,nStates),...
                   	 'Pseudotransitions', ones(nStates));
       
      % ==================================
      % Initialize state model parameters
      % ==================================
      startParams.nStates                   = nStates;
      startParams                           =...
        renamefield(startParams, {'ST', 'TR'}, {'pi', 'trans'});
      startParams                           = rmfield(startParams,'E');
      startParams.piPrior                   = ones(1, startParams.nStates);
      startParams.transPrior                = ones(startParams.nStates);

      % ========================================
      % Initialize observation model parameters
      % ========================================
      startParams.d                         = obj.mu';
      if isequal(CovType, 'full')
        startParams.R                      	= obj.Sigma;
      elseif isequal(CovType, 'diagonal')
        startParams.R                      	= nddiag(obj.Sigma);
      end
    elseif ~any(~initStatus)
      fprintf('Initializing parameters with input...\n');
      % ==================================
      % Initialize state model parameters
      % ==================================
      startParams.nStates                   = nStates;
      startParams.pi                        = pi0;
      startParams.piPrior                   = piPrior;
      startParams.trans                     = trans0;
      startParams.transPrior                = transPrior;

      % ========================================
      % Initialize observation model parameters
      % ========================================
      startParams.d                         = d;
      startParams.R                         = R;
    else
      error(['d, R, pi, piPrior, trans and transPrior must',...
             ' all be specified if there is no data for training.']);
    end
    
    % Define parameter constraints
    startParams.covType                     = CovType;
    startParams.sharedCov                   = SharedCov;

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting HMM model...\n');

      [estParams, seqTrainCut, LL, iterTime]=...
        em_hmm(startParams, seqTrainCut, extraOpts{:});

      % Extract neural trajectories for original, unsegmented trials
      % using learned parameters
      [seqTrain, ess, LLorig]               =...
        exactInferenceWithLL_hmm(seqTrain, estParams,...
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
        % Extract neural trajectories for original, unsegmented trials
        % using learned parameters
        [seqTrain, ess, LLorig]            	=...
          exactInferenceWithLL_hmm(seqTrain, estParams,...
                                   'getSeq', true,...
                                   'getLL', true,...
                                   extraOpts{:});
      end % if ~isempty(seqTrain)
    end % if (inference)
  end

  if ~isempty(seqTest) % check if there are any test trials
    % =============================================
    % Leave-one-channel-out prediction on test data
    % =============================================
    if (prediction)
      seqTest                               =...
        predict_hmm(seqTest, estParams, varargin{:});
    end % if (prediction)
    
    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      [seqTest, ~, LLtest]                  =...
        exactInferenceWithLL_hmm(seqTest, estParams,...
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
    save(fname, vars{~ismember(vars, {'yAll', 'mixComp', 'temp'})});
  elseif inference || prediction
    if (inference)
      fprintf('Saving inferred latent variables in %s...\n', fname);
      if ~isempty(seqTrain)
        save(fname, 'seqTrain', 'ess', 'LLorig', '-append');
      end % if ~isempty(seqTrain)
      if ~isempty(seqTest)
        save(fname, 'seqTest', 'LLtest', '-append');
      end % if ~isempty(seqTest)
      save(fname, 'inference', '-append');
    end % if (inference)

    if (prediction) && ~isempty(seqTest)
      fprintf('Saving leave-one-channel-out prediction in %s...\n', fname);
      if (inference)
        save(fname, 'prediction', '-append');
      else % if (~inference)
        save(fname, 'seqTest', 'prediction', '-append');
      end
    end
  end
end