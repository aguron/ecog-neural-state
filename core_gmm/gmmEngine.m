function gmmEngine(seqTrain, seqTest, fname, varargin)
%
% gmmEngine(seqTrain, seqTest, fname, ...) 
%
% Extract neural trajectories using GMM.
%
%   yDim: number of electrodes
%
% INPUTS:
%
% seqTrain              - training data structure, whose nth
%                         entry (corresponding to the nth
%                         experimental trial) has fields
%                           trialId (1 x 1)       -- unique trial
%                                                    identifier
%                           y (# electrodes x T)	-- neural data
%                           T (1 x 1)            	-- number of timesteps
% seqTest               - test data structure (same format as seqTrain)
% fname                 - filename of where results are saved
%
% OPTIONAL ARGUMENTS:
%
% nMixComp             	- number of mixture components (default: 3)
% binWidth            	- ECoG window width in seconds (default: 0.2)
% kernSD                - Gaussian smoothing kernel width in seconds
%                         (default: 0)

% d (yDim x nMixComp)  	- observation mean(s)
% R (yDim x yDim x    	- observation noise covariance(s)
%    nMixComp (or 1))
% Pi (1 x nMixComp)   	- mixture component priors
% mixCompGuess         	- initial component guesses for time points
%                         for training data

% nInitMax            	- number of times k-means is run to satisfy
%                         minimum cluster size threshold (default: 5)
% clSzThr              	- minimum cluster size threshold (default: yDim)

% Start                	- 'randSample' (default), 'initParams', 'InitComp'
% Replicates            - number of times to repeat EM algorithm
%                         (default: 1)
% CovType             	- 'full' or 'diagonal' (default)
% SharedCov            	- true or false (default)
% Regularize          	- nonnegative number added to diagonals of
%                         covariance matrices to make them
%                         positive-definite (default: 0)
% Options               - Please see GMDISTRIBUTION.FIT MATLAB function

% prediction            - specifies whether leave-channel-out prediction
%                         should be carried out (default: false)
%
% @ 2016 Akinyinka Omigbodun    aomigbod@ucsd.edu

  nMixComp                    = 3;
  binWidth                    = 0.2;	% in sec
  kernSD                      = 0;    % in sec

  d                           = [];
  R                           = [];
  Pi                          = [];
  mixCompGuess                = [];
  
  nInitMax                    = 5;
  clSzThr                     = [];

  Start                       = 'randSample';
  Replicates                  = 1;
  CovType                     = 'diagonal';
  SharedCov                   = false;
  Regularize                  = 0;
  Options                     = [];

  learning                   	= true;
  inference                  	= true;
  prediction                 	= false;

  extraOpts                   = assignopts(who, varargin);

  if (learning)
    if (kernSD)
      fprintf('Performing presmoothing with kernSD=%g and GMM...\n',...
              kernSD);
      % ========================
      % Presmooth data over time
      % ========================

      % Training data
      seqTrainOrig            = seqTrain;
      for n=1:length(seqTrain)
        seqTrain(n).y         = smoother(seqTrain(n).y, kernSD, binWidth);
      end % for n=1:length(seqTrain)

      % Test data
      seqTestOrig             = seqTest;
      for n=1:length(seqTest)
        seqTest(n).y          = smoother(seqTest(n).y, kernSD, binWidth);
      end % for n=1:length(seqTest)
    end % if (kernSD)

    % For compute efficiency, train on equal-length segments of trials
    seqTrainCut               = cutTrials(seqTrain, extraOpts{:});
    if isempty(seqTrain)
      fprintf('No segments extracted for training.\n');
    elseif isempty(seqTrainCut)
      fprintf(['WARNING: no segments extracted for training.',...
               ' Defaulting to segLength=Inf.\n']);
      seqTrainCut             = cutTrials(seqTrain, 'segLength', Inf);
    end % if isempty(seqTrainCut)

    yAll                      = [seqTrainCut.y];
    if ~isempty(yAll)
      yDim                   	= size(yAll, 1);
    else % if isempty(yAll)
      yDim                    = size([seqTest.y], 1);
    end

    switch(Start)
      case 'randSample'
        if isempty(seqTrainCut)
          error('No training data.');
        end % if isempty(seqTrainCut)
        start                 = Start;
      case 'initParams'
        initStatus            = ~[isempty(d) isempty(R) isempty(Pi)];
        if all(~initStatus) && ~isempty(seqTrainCut)
          fprintf('Initializing parameters using K-means...\n');
          if isempty(clSzThr)
            clSzThr           = yDim;
          end % if isempty(clSzThr)
          [temp2, temp]       =...
            kmeansclszthr(yAll',nMixComp,clSzThr,nInitMax);
          start.mu            = temp;

          if strcmp(CovType, 'diagonal')
           start.Sigma        = nan(1,yDim,nMixComp);
           for j=1:nMixComp
            start.Sigma(:,:,j)= var(yAll(:,temp2==j),1,2);
           end % for j=1:nMixComp
          elseif strcmp(CovType, 'full')
           start.Sigma        = nan(yDim,yDim,nMixComp);
           for j=1:nMixComp
            start.Sigma(:,:,j)= cov(yAll(:,temp2==j)',1);
           end % for j=1:nMixComp
          end
          if (SharedCov)
            clsz(1,1,:)      	= uhistc(temp2);
            start.Sigma(:,:,1)=...
              sum(bsxfun(@times,start.Sigma,clsz),3)/numel(temp2);
            start.Sigma(:,:,2:end)...
                              = [];
          end % if (SharedCov)
          start.PComponents   = uhistc(temp2)/numel(temp2);
        elseif ~any(~initStatus)
          fprintf('Initializing parameters with input...\n');
          start.mu            = d';

          if strcmp(CovType, 'diagonal')
           start.Sigma        = nan(1,yDim,max(1,nMixComp*~SharedCov));
           for j=1:max(1,nMixComp*~SharedCov)
            start.Sigma(:,:,j)= diag(R(:,:,j));
           end % for j=1:max(1,nMixComp*SharedCov)
          elseif strcmp(CovType, 'full')
            start.Sigma      	= R;
          end
          start.PComponents   = Pi;
        else
          error(['d, R, and Pi must all be specified',...
                 ' if there is no data for training.']);
        end % if all(~initStatus) && ~isempty(seqTrainCut)
      case 'InitComp'
        if isempty(seqTrainCut)
          error('No training data.');
        end % if isempty(seqTrainCut)
        if isempty(mixCompGuess)
          if isempty(clSzThr)
            clSzThr           = yDim;
          end % if isempty(clSzThr)
          start               =...
            kmeansclszthr(yAll',nMixComp,clSzThr,nInitMax);
        else
          start               = cell2mat(mixCompGuess);
        end
      otherwise
        error('Invalid Start specification');
    end % switch(Start)

    if ~isempty(seqTrainCut)
      % =====================
      % Fit model parameters
      % =====================
      fprintf('\nFitting GMM model...\n');

      obj                   	=...
        gmdistribution.fit(yAll', nMixComp,...
                           'Start', start,...
                           'Replicates', Replicates,...
                           'CovType', CovType,...
                           'SharedCov', SharedCov,...
                           'Regularize', Regularize,...
                           'Options', Options);

      % Extract neural trajectories for original, unsegmented trials
      % using learned parameters
      temp                   	= cluster(obj,[seqTrain.y]');
      seqTrain               	=...
        segmentByTrial(seqTrain, temp', 'mixComp');
      [~, LLorig]            	= posterior(obj,[seqTrain.y]');
      LLorig                 	= -LLorig;
    else % if isempty(seqTrainCut)
      obj                     =...
        gmdistribution(start.mu, start.Sigma, start.PComponents);
    end
  else % if (~learning)
    obj                       = loadvars(fname, 'obj');
    if (inference)
      if ~isempty(seqTrain)
        % Extract neural trajectories for original, unsegmented trials
        % using learned parameters
        temp                 	= cluster(obj,[seqTrain.y]');
        seqTrain             	=...
          segmentByTrial(seqTrain, temp', 'mixComp');
        [~, LLorig]          	= posterior(obj,[seqTrain.y]');
        LLorig               	= -LLorig;
      end % if ~isempty(seqTrain)
    end % if (inference)
  end

  if ~isempty(seqTest) % check if there are any test trials
    % =========================================
    % Leave-channel-out prediction on test data
    % =========================================
    if (prediction)
      seqTest                 = predict_gmm(seqTest, obj, varargin{:});
    end % if (prediction)

    % ===================================================
    % Neural trajectories and loglikelihood for test data
    % ===================================================
    if (inference)
      temp                    = cluster(obj,[seqTest.y]');
      seqTest                 =...
        segmentByTrial(seqTest, temp', 'mixComp');
      [~, LLtest]             = posterior(obj,[seqTest.y]');
      LLtest                  = -LLtest;
    end % if (inference)
  end % if ~isempty(seqTest)

  % =============
  % Save results
  % =============
  if (learning)
    vars                      = who;
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
