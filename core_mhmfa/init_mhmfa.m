function [startParams, mixCompGuess, outliers] =...
      init_mhmfa(seq, fname, varargin)
%
% init_mhmfa(seq, fname, ...)
%
% Initialization for MHMFA.
%
%   yDim: number of electrodes
%
% INPUTS:
%
% seq         - data structure, whose nth entry (corresponding
%               to the nth experimental trial) has fields
%                 trialId (1 x 1)         -- unique trial identifier
%                 y (# electrodes x T)    -- neural data
%                 T (1 x 1)               -- number of timesteps
% fname       - model filename
%
% OPTIONAL ARGUMENTS:
%
% xDim        - state dimensionality (default: 3)
% nStates     - number of HMFA states (default: 3)
% faType      - HMFA factor analyzer(s) (with tied (0) or untied (1)
%               mean, factor loading matrix, and covariance parameters)
%               (default: [1 1 1])
% nMixComp    - number of MHMFA mixture components (default: 3)

% stateGuess                      - initial state guesses for time points
%                                 	for training data
% mixCompGuess                    - initial MHMFA component guesses
%                                   for training data: use an index of
%                                   0 if there is no guess for a trial
%                                   mixComp; an index between 1 and
%                                   nMixComp inclusive, if the mixComp
%                                   guess is to be used in the
%                                   initialization; and an index between
%                                   -1 and -nMixComp inclusive, if the
%                                   mixComp guess is not to be used in the
%                                   initialization
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  xDim                                      = 3;
  nStates                                   = 3;
  faType                                    = [1 1 1];
  nMixComp                                  = 3;

  Replicates                                = 1;
  Regularize                               	= 0;
  Options                                   = [];

  stateGuess                                = [];
  mixCompGuess                              = [];
  outlierThr                                = 0.5;
  skipSeq                                   = true;

  extraOpts                                 = assignopts(who, varargin);

  if ~any(faType) && (nStates > 1)
    error(['At least one of the factor analyzer parameters must ',...
           'be untied if nStates > 1']);
  end % if ~any(faType) && (nStates > 1)

  if ~isequal(faType,[1 1 1]) &&...
     ~isequal(faType,[1 1 0])
    fprintf('Does not support faType = [%d %d %d]\n', faType);
  end
  
  if isempty(seq)
    error('No trials for initialization');
  end % if isempty(seq)
  yDim                                      = size(seq(1).y, 1);
  nSeq                                      = numel(seq);

  faTypeSpec                                = 'tu'; % t - tied; u - untied
  [~, str]                                  = rstrtok(fname,'_');
  fname                                     =...
   sprintf('%s/init_mhmfa_xDim%02d_nStates%02d_MFV%c%c%c',...
           rstrtok(fname,'/'), xDim, nStates, faTypeSpec(faType+1));
  if ~isnan(str2double(str(end)))
    fname                                   = sprintf('%s%s', fname, str);
  end % if ~isnan(str2double(str(end)))

  if exist([fname '.mat'], 'file')
    if isempty(mixCompGuess)
      [L_symm,estParamsMFA,seqMFA,omitted]  =...
        loadvars(fname, 'L_symm', 'estParamsMFA', 'seqMFA', 'omitted');
    else % if ~isempty(mixCompGuess)
      [estParamsMFA,seqMFA,omitted]        	=...
        loadvars(fname, 'estParamsMFA', 'seqMFA', 'omitted');
    end
  else % if ~exist([fname '.mat'], 'file')
    fprintf('Initializing parameters...\n');

    fprintf('First step: GMM...\n');
    args                                    =...
      {'Replicates', Replicates,...
       'CovType', 'diagonal', 'SharedCov', ~faType(3),...
       'Regularize', Regularize, 'Options', Options};
    if ~isempty(stateGuess)
      args                                  =...
        [args, 'Start', cell2mat(stateGuess)];
    end % if ~isempty(stateGuess)
    obj                                     =...
      gmdistribution.fit([seq.y]', nStates, args{:});

    fprintf('Second step: MFA...\n');
    startParamsMFA.nMixComp                 = nStates;
    startParamsMFA.faType                   = faType;
    startParamsMFA.Pi                       = obj.PComponents;
    startParamsMFA.d                        = obj.mu';
    for j=1:nStates
      jC                                    = j*faType(2)+(1 - faType(2));
      jR                                    = j*faType(3)+(1 - faType(3));

      startParamsMFA.C(:,:,jC)              = eye(yDim,xDim);
      startParamsMFA.R(:,:,jR)              = diag(obj.Sigma(:,:,jR));
    end % for j=1:nStates

    startParamsMFA.notes.RforceDiagonal     = true;

    [estParamsMFA, seqMFA]                  =...
      em_mfa(startParamsMFA, seq, extraOpts{:});

    fprintf(['HMFA for computation of ',...
             'loglikelihood "distance" matrix...\n']);
    mixCompMFA                              =	cell(1,numel(seqMFA));
    [mixCompMFA{:}]                         = deal(seqMFA.mixComp);

    hmm                                     =...
      hmmestimate2(mixCompMFA, mixCompMFA, nStates, nStates,...
                   'Pseudostarts', ones(1,nStates),...
                   'Pseudotransitions', ones(nStates));
    startParamsHMFA.nStates                 = nStates;
    startParamsHMFA.faType                  = faType;
    startParamsHMFA.pi                      = hmm.ST;
    startParamsHMFA.trans                   = hmm.TR;
    startParamsHMFA.piPrior                 = ones(1, nStates);
    startParamsHMFA.transPrior              = ones(nStates);

    startParamsHMFA.d                       = estParamsMFA.d;
    startParamsHMFA.C                       = estParamsMFA.C;
    startParamsHMFA.R                       = estParamsMFA.R;

    startParamsHMFA.notes.RforceDiagonal    = true;

    L                                       = nan(nSeq);
    for i=1:nSeq
      try
        fprintf('Sequence %d...\n',i);
        estParamsHMFA                      	=...
          em_hmfa(startParamsHMFA, seq(i), extraOpts{:});
        % Compute unsymmetrized loglikelihood "distance" matrix
        for j=1:nSeq
          [~, ~, L(i,j)]                  	=...
            exactInferenceWithLL_hmfa(seq(j), estParamsHMFA,...
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
    fprintf('Saving %s...\n', fname);
    vars                                    = who;
    save(fname, vars{~ismember(vars, {'faTypeSpec', 'mixCompMFA',...
                                      'hmm'})});
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

  if ~exist('startParamsHMFA', 'var')
    mixCompMFA                              =	cell(1,numel(seqMFA));
    [mixCompMFA{:}]                         = deal(seqMFA.mixComp);

    hmm                                     =...
      hmmestimate2(mixCompMFA, mixCompMFA, nStates, nStates,...
                   'Pseudostarts', ones(1,nStates),...
                   'Pseudotransitions', ones(nStates));
    startParamsHMFA.nStates                 = nStates;
    startParamsHMFA.faType                  = faType;
    startParamsHMFA.pi                      = hmm.ST;
    startParamsHMFA.trans                   = hmm.TR;
    startParamsHMFA.piPrior                 = ones(1, nStates);
    startParamsHMFA.transPrior              = ones(nStates);

    startParamsHMFA.d                       = estParamsMFA.d;
    startParamsHMFA.C                       = estParamsMFA.C;
    startParamsHMFA.R                       = estParamsMFA.R;

    startParamsHMFA.notes.RforceDiagonal    = true;
  end % if ~exist('startParamsHMFA', 'var')


  mixCompGuess(outliers)                    = -abs(mixCompGuess(outliers));
  for k=1:nMixComp
    fprintf('MHMFA Component %d Parameter Initializations...\n', k);
    % ==================================
    % Initialize state model parameters
    % ==================================
    if (k == 1)
      startParams(1).nStates               	= nStates;
      startParams(1).nMixComp               = nMixComp;
    end % if (k == 1)

    estParamsHMFA                           =...
      em_hmfa(startParamsHMFA, seq(mixCompGuess==k), extraOpts{:});

    startParams(k).Pi                       =...
      (sum(mixCompGuess==k) + 1)/(sum(mixCompGuess>0) + nMixComp);
    startParams(k).pi                       = estParamsHMFA.pi;
    startParams(k).trans                    = estParamsHMFA.trans;
    startParams(k).piPrior                  = ones(1, nStates);
    startParams(k).transPrior               = ones(nStates);

    % ========================================
    % Initialize observation model parameters
    % ========================================
    startParams(k).d                        = estParamsHMFA.d;
    startParams(k).C                        = estParamsHMFA.C;
    startParams(k).R                        = estParamsHMFA.R;
    
    % Define parameter constraints
    if (k == 1)
      startParams(1).faType                	= faType;
      startParams(1).notes.RforceDiagonal  	= true;
    end % if (k == 1)
  end % for k=1:nMixComp
end