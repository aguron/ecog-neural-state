function [startParams, mixCompGuess, outliers] =...
      init_mhmm(seq, fname, varargin)
%
% init_mhmm(seq, fname, ...)
%
% MHMM initialization
%
%   yDim: number of electrodes or channels
%
% INPUTS:
%
% seq           - data structure, whose n-th entry (corresponding
%                 to the n-th experimental trial) has fields
%                   trialId           -- unique trial identifier
%                   trialType (1 x 1)	-- trial type index (Optional)
%                   fs (1 x 1)       	-- sampling frequency of ECoG data
%                   T (1 x 1)         -- number of timesteps
%                   y (yDim x T)      -- neural data
% fname         - model filename
%
% OUTPUTS:
%
% startParams   - parameters from initialization
% mixCompGuess 	- vector of initial component HMM guesses for training
%                 trials. Use an index of 0 if there is no guess for a 
%                 trial; an index between 1 and nMixComp (inclusive), if
%                 the trial is to be used in initialization; and an index
%                	between -nMixComp and -1 (inclusive), if the trial is not 
%                 to be used in model fitting
% outliers      - trials identified as outliers. Index corresponds to
%                 position in seq
%
% OPTIONAL ARGUMENTS:
%
% nMixComp      - number of MHMM mixture components (default: 3)
% nStates       - number of Markov (HMM) states (default: 3)
% CovType       - component HMM covariance type: 'full' or
%                 'diagonal' (default: 'diagonal')
% SharedCov     - component HMM covariance tied (true) or untied (false)
%                 (default: false)
%
% Replicates   ]  
% Regularize   ]- parameters for MATLAB function GMDISTRIBUTION.FIT
% Options      ]  
%
% stateGuess    - cell array of initial state guesses for time points of
%                 trials for training data
%
% mixCompGuess 	- same as mixCompGuess above
% outlierThr    - if the fraction of trials that a trial is the most
%                 dissimilar from is below this threshold, that trial
%                 is not designated as an outlier
% skipSeq       - if true, only one attempt is made to fit an HMM to
%                 each trial, and trials with unsuccessful model fitting
%                 first attempts are skipped; if false, any unsuccessful
%                 model fit results in an error
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  nMixComp                                  = 3;
  nStates                                   = 3;
  CovType                                   = 'diagonal';
  SharedCov                                 = false;

  Replicates                                = 1;
  Regularize                               	= 0;
  Options                                   = [];

  stateGuess                                = [];

  mixCompGuess                              = [];
  outlierThr                                = 0.5;
  skipSeq                                   = true;

  extraOpts                                 = assignopts(who, varargin);

  if isempty(seq)
    error('No trials for initialization');
  end % if isempty(seq)
  yDim                                      = size(seq(1).y, 1);
  nSeq                                      = numel(seq);

  [~, str]                                  = rstrtok(fname,'_');
  fname                                     =...
   sprintf('%s/init_mhmm_nStates%02d_%s',...
           rstrtok(fname,'/'), nStates, CovType);
  if ~isnan(str2double(str(end)))
    fname                                   = sprintf('%s%s', fname, str);
  end % if ~isnan(str2double(str(end)))

  if exist([fname '.mat'], 'file')
    if isempty(mixCompGuess)
      [L_symm, obj, seqGMM, omitted]        =...
        loadvars(fname, 'L_symm', 'obj', 'seqGMM', 'omitted');
    else % if ~isempty(mixCompGuess)
      [obj, seqGMM, omitted]                =...
        loadvars(fname, 'obj', 'seqGMM', 'omitted');
    end
  else % if ~exist([fname '.mat'], 'file')
    fprintf('Initializing parameters...\n');

    fprintf('First step: GMM...\n');
    args                                    =...
      {'Replicates', Replicates,...
       'CovType', CovType, 'SharedCov', SharedCov,...
       'Regularize', Regularize, 'Options', Options};
    if ~isempty(stateGuess)
      args                                  =...
        [args, 'Start', cell2mat(stateGuess)];
    end % if ~isempty(stateGuess)
    obj                                     =...
      gmdistribution.fit([seq.y]', nStates, args{:});
    seqGMM                                  =...
      segmentByTrial(seq, cluster(obj,[seq.y]')', 'mixComp');

    fprintf(['HMM for computation of ',...
             'loglikelihood "distance" matrix...\n']);
    mixCompGMM                            	= cell(1,numel(seqGMM));
    [mixCompGMM{:}]                        	= deal(seqGMM.mixComp);
    hmm                                     =...
      hmmestimate2(mixCompGMM, mixCompGMM, nStates, nStates,...
                   'Pseudostarts', ones(1,nStates),...
                   'Pseudotransitions', ones(nStates));
    startParamsHMM.nStates                  = nStates;
    startParamsHMM.covType                 	= CovType;
    startParamsHMM.sharedCov               	= SharedCov;
    startParamsHMM.pi                       = hmm.ST;
    startParamsHMM.trans                    = hmm.TR;
    startParamsHMM.piPrior                  = ones(1, nStates);
    startParamsHMM.transPrior               = ones(nStates);

    startParamsHMM.d                      	= obj.mu';
    if isequal(CovType, 'full')
      startParamsHMM.R                    	= obj.Sigma;
    elseif isequal(CovType, 'diagonal')
      startParamsHMM.R                    	= nddiag(obj.Sigma);
    end

    L                                       = nan(nSeq);
    for i=1:nSeq
      try
        fprintf('Sequence %d...\n',i);
        estParamsHMM                       	=...
          em_hmm(startParamsHMM, seq(i), extraOpts{:});

        % Compute unsymmetrized loglikelihood "distance" matrix
        for j=1:nSeq
          [~, ~, L(i,j)]                      =...
            exactInferenceWithLL_hmm(seq(j), estParamsHMM,...
                                     'getLL', true,...
                                     extraOpts{:});
        end % for j=1:nSeq
      catch err
        displayerror(err)
        if (~skipSeq)
          rethrow(err)
        end % if (~skipSeq)
      end % try
    end % for i=1:nSeq

    % Symmetrize L
    L_symm                                  = symm(L);
    omitted                                 = find(isnan(diag(L_symm)))';

    % Saving variables
    fprintf('Saving %s...\n\n', fname);
    vars                                    = who;
    save(fname, vars{~ismember(vars, {'hmm'})});
  end

  if isempty(mixCompGuess)
    mixCompGuess                            = zeros(1,nSeq);
    % Remove outliers
    L_symm                                  =...
      removeoutliertrials(L_symm, ceil(outlierThr*nSeq));

    % Select trials for initialization
    % row - model
    % col - probability
    [m p]                                  	=...
     find(L_symm == min(L_symm(:)));
    idx(1:2)                              	= [m(1) p(1)];
    for k=3:nMixComp
      d                                   	= nan(1,nSeq);
      for i=setdiff(1:nSeq,idx)
        d(i)                               	=...
          sum(L_symm(sub2ind(size(L_symm),idx,i*ones(1,numel(idx)))));
      end % for i=setdiff(1:nSeq,idx)
      [~, idx(k)]                           = min(d);
    end % for k=3:nMixComp

    for k=1:nMixComp
      mixCompGuess(idx(k))                  = k;
    end % for k=1:nMixComp
  end % if isempty(mixCompGuess)

  if (numel(mixCompGuess) ~= nSeq)
    error('Invalid mixCompGuess specification');
  end % if (numel(mixCompGuess) ~= nSeq)

  outliers                                  = [];
  if (nargout == 3)
    if ~exist('L_symm','var')
      L_symm                                = loadvars(fname, 'L_symm');
      % Remove outliers
      L_symm                               	= removeoutliertrials(L_symm);
    end % if ~exist('L_symm','var')
    
    outliers                                = find(isnan(diag(L_symm)))';
    outliers                                = setdiff(outliers,omitted);
  end % if (nargout == 3)

  if exist('L_symm','var') && (numel(mixCompGuess) ~= size(L_symm, 1))
    error('mixCompGuess and L_symm are incompatible');
  end

  if ~exist('startParamsHMM', 'var')
    mixCompGMM                            	= cell(1,numel(seqGMM));
    [mixCompGMM{:}]                        	= deal(seqGMM.mixComp);
    hmm                                     =...
      hmmestimate2(mixCompGMM, mixCompGMM, nStates, nStates,...
                   'Pseudostarts', ones(1,nStates),...
                   'Pseudotransitions', ones(nStates));
    startParamsHMM.nStates                  = nStates;
    startParamsHMM.covType                 	= CovType;
    startParamsHMM.sharedCov                = SharedCov;
    startParamsHMM.pi                       = hmm.ST;
    startParamsHMM.trans                    = hmm.TR;
    startParamsHMM.piPrior                  = ones(1, nStates);
    startParamsHMM.transPrior               = ones(nStates);

    startParamsHMM.d                      	= obj.mu';
    if isequal(CovType, 'full')
      startParamsHMM.R                    	= obj.Sigma;
    elseif isequal(CovType, 'diagonal')
      startParamsHMM.R                    	= nddiag(obj.Sigma);
    end
  end % if ~exist('startParamsHMM', 'var')

  mixCompGuess(outliers)                    = -abs(mixCompGuess(outliers));
  for k=1:nMixComp
    fprintf('MHMM Component %d Parameter Initializations...\n', k);
    % ==================================
    % Initialize state model parameters
    % ==================================
    if (k == 1)
      startParams(1).nStates               	= nStates;
      startParams(1).nMixComp               = nMixComp;
    end % if (k == 1)

    estParamsHMM                            =...
      em_hmm(startParamsHMM, seq(mixCompGuess==k), extraOpts{:});

    startParams(k).Pi                       =...
      (sum(mixCompGuess==k) + 1)/(sum(mixCompGuess>0) + nMixComp);
    startParams(k).pi                       = estParamsHMM.pi;
    startParams(k).trans                    = estParamsHMM.trans;
    startParams(k).piPrior                  = ones(1, nStates);
    startParams(k).transPrior               = ones(nStates);

    % ========================================
    % Initialize observation model parameters
    % ========================================
    startParams(k).d                      	= estParamsHMM.d;
    startParams(k).R                        = estParamsHMM.R;

    % Define parameter constraints
    if (k == 1)
      startParams(1).covType               	= CovType;
      startParams(1).sharedCov             	= SharedCov;
    end % if (k == 1)
  end % for k=1:nMixComp
end