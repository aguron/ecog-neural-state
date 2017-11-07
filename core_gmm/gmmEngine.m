function gmmEngine(seqTrain, seqTest, fname, varargin)
%
% gmmEngine(seqTrain, seqTest, fname, ...) 
%
% Model fitting and inference with GMM
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
% nMixComp                          - number of mixture components
%                                     (default: 3)
% CovType                           - GMM covariance type: 'full' or
%                                     'diagonal' (default: 'diagonal')
% SharedCov                         - GMM covariance tied (true) or
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
% d (yDim x nMixComp)              	- observation means
% R (yDim x yDim x nMixComp (or 1))	- covariance(s)
%
% Pi (1 x nMixComp)                	- mixture component priors
%
% mixCompGuess                    	- cell array of initial component
%                                     guesses for time points of trials for
%                                     training data
%
% Start                             - 'randSample' (default): please see
%                                     information on 'randSample' value of
%                                     'Start' parameter for MATLAB function
%                                     GMDISTRIBUTION.FIT
%
%                                     'initParams': initialization with
%                                     specified parameters or parameters
%                                     derived from k-means
%
%                                     'InitComp': initialization with
%                                     specified initial component guesses
%                                     for time points of trials for
%                                     training data or with initial
%                                     component guesses obtained with
%                                     k-means
% nInitMax                          - maximum number of times k-means is
%                                     run attempting to satisfy a minimum
%                                     time point cluster size threshold
%                                     (default: 5)
% clSzThr                           - minimum time point cluster size
%                                     threshold (set to yDim if not
%                                     specified as an optional argument)
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
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  nMixComp                          = 3;
  CovType                           = 'diagonal';
  SharedCov                         = false;

  binWidth                          = 0.2;
  kernSD                            = 0;

  Replicates                        = 1;
  Regularize                        = 0;
  Options                           = [];

  d                                 = [];
  R                                 = [];

  Pi                                = [];

  mixCompGuess                      = [];

  Start                             = 'randSample';
  nInitMax                          = 5;
  clSzThr                           = [];

  learning                          = true;
  inference                         = true;
  prediction                        = false;

  extraOpts                         = assignopts(who, varargin);

  if (learning)
    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and GMM...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig                  = seqTrain;
      for n=1:length(seqTrain)
        seqTrain(n).y               =...
         smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:length(seqTrain)

      % Test data
      seqTestOrig                   = seqTest;
      for n=1:length(seqTest)
        seqTest(n).y                =...
         smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:length(seqTest)
    end % if (kernSD)

    % Resegment trials for training
    [seqTrainCut, resegTrlGuess]    =...
      resegmenttrials(seqTrain,...
                      'method', 'gmm',...
                      'mixCompGuess', mixCompGuess,...
                      extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      [seqTrainCut, resegTrlGuess]	=...
        resegmenttrials(seqTrain,...
                        'method', 'gmm',...
                        'mixCompGuess', mixCompGuess,...
                        'segLength', Inf);
    end

    yAll                            = [seqTrainCut.y];
    if ~isempty(yAll)
      yDim                          = size(yAll, 1);
    else % if isempty(yAll)
      yDim                          = size([seqTest.y], 1);
    end

    switch(Start)
      case 'randSample'
        if isempty(seqTrainCut)
          error('No training data.');
        end % if isempty(seqTrainCut)
        start                       = Start;
      case 'initParams'
        initStatus                  = ~[isempty(d) isempty(R) isempty(Pi)];
        if all(~initStatus) && ~isempty(seqTrainCut)
          fprintf('Initializing parameters using K-means...\n');
          if isempty(clSzThr)
            clSzThr                 = yDim;
          end % if isempty(clSzThr)
          [temp2, temp]             =...
            kmeansclszthr(yAll',nMixComp,clSzThr,nInitMax);
          start.mu                  = temp;

          if strcmp(CovType, 'diagonal')
           start.Sigma              = nan(1,yDim,nMixComp);
           for j=1:nMixComp
            start.Sigma(:,:,j)      = var(yAll(:,temp2==j),1,2);
           end % for j=1:nMixComp
          elseif strcmp(CovType, 'full')
           start.Sigma              = nan(yDim,yDim,nMixComp);
           for j=1:nMixComp
            start.Sigma(:,:,j)      = cov(yAll(:,temp2==j)',1);
           end % for j=1:nMixComp
          end
          if (SharedCov)
            clsz(1,1,:)             = uhistc(temp2);
            start.Sigma(:,:,1)      =...
              sum(bsxfun(@times,start.Sigma,clsz),3)/numel(temp2);
            start.Sigma(:,:,2:end)  = [];
          end % if (SharedCov)
          start.PComponents         = uhistc(temp2)/numel(temp2);
        elseif ~any(~initStatus)
          fprintf('Initializing parameters with input...\n');
          start.mu                  = d';

          if strcmp(CovType, 'diagonal')
           start.Sigma              =...
            nan(1,yDim,max(1,nMixComp*~SharedCov));
           for j=1:max(1,nMixComp*~SharedCov)
            start.Sigma(:,:,j)      = diag(R(:,:,j));
           end % for j=1:max(1,nMixComp*SharedCov)
          elseif strcmp(CovType, 'full')
            start.Sigma             = R;
          end
          start.PComponents         = Pi;
        else
          error(['d, R, and Pi must all be specified',...
                 ' if there is no data for training.']);
        end % if all(~initStatus) && ~isempty(seqTrainCut)
      case 'InitComp'
        if isempty(seqTrainCut)
          error('No training data.');
        end % if isempty(seqTrainCut)
        if isempty(resegTrlGuess.mixComp)
          if isempty(clSzThr)
            clSzThr                 = yDim;
          end % if isempty(clSzThr)
          start                     =...
            kmeansclszthr(yAll',nMixComp,clSzThr,nInitMax);
        else
          start                     = cell2mat(resegTrlGuess.mixComp);
        end
      otherwise
        error('Invalid Start specification');
    end % switch(Start)

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting GMM model...\n');

      obj                           =...
        gmdistribution.fit(yAll', nMixComp,...
                           'Start', start,...
                           'Replicates', Replicates,...
                           'CovType', CovType,...
                           'SharedCov', SharedCov,...
                           'Regularize', Regularize,...
                           'Options', Options);

      % Extract neural trajectories for original, unsegmented trials
      % using learned parameters
      temp                          = cluster(obj,[seqTrain.y]');
      seqTrain                      =...
        segmentByTrial(seqTrain, temp', 'mixComp');
      [~, LLorig]                   = posterior(obj,[seqTrain.y]');
      LLorig                        = -LLorig;
    else % if isempty(seqTrainCut)
      obj                           =...
        gmdistribution(start.mu, start.Sigma, start.PComponents);
    end
  else % if (~learning)
    obj                             = loadvars(fname, 'obj');
    if (inference)
      if ~isempty(seqTrain)
        % Extract neural trajectories for original, unsegmented trials
        % using learned parameters
        temp                        = cluster(obj,[seqTrain.y]');
        seqTrain                    =...
          segmentByTrial(seqTrain, temp', 'mixComp');
        [~, LLorig]                 = posterior(obj,[seqTrain.y]');
        LLorig                      = -LLorig;
      end % if ~isempty(seqTrain)
    end % if (inference)
  end

  if ~isempty(seqTest) % check if there are any test trials
    % =========================================
    % Leave-channel-out prediction on test data
    % =========================================
    if (prediction)
      seqTest                       =...
       predict_gmm(seqTest, obj, varargin{:});
    end % if (prediction)

    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      temp                          = cluster(obj,[seqTest.y]');
      seqTest                       =...
        segmentByTrial(seqTest, temp', 'mixComp');
      [~, LLtest]                   = posterior(obj,[seqTest.y]');
      LLtest                        = -LLtest;
    end % if (inference)
  end % if ~isempty(seqTest)

  % =============
  % Save results
  % =============
  if (learning)
    vars                            = who;
    fprintf('Saving %s...\n', fname);
    save(fname, vars{~ismember(vars, {'yAll', 'temp', 'temp2'})});
  elseif inference || prediction
    if (inference)
      fprintf('Saving inferred latent variables in %s...\n', fname);
      if ~isempty(seqTrain)
        save(fname, 'seqTrain', 'LLorig', '-append');
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
