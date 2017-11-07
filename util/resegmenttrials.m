function [seqOut, cutTrialGuess] = resegmenttrials(seqIn, varargin)
%
% seqOut = resegmenttrials(seqIn, ...)  
%
% Extracts trial segments that are all of the same length.  Uses
% overlapping segments if trial length is not integer multiple
% of segment length. Ignores trials with length shorter than 
% one segment length.
%
% INPUTS:
%
% seqIn         - data structure, whose nth entry (corresponding to
%                 the nth experimental trial) has fields
%                   trialId       -- unique trial identifier
%                   trialType     -- trial type index (Optional)
%                   fs            -- sampling frequency of ECoG data
%                   T (1 x 1)     -- number of timesteps
%                   y (yDim x T)  -- neural data
%
% OUTPUTS:
%
% seqOut        - data structure, whose nth entry (corresponding to
%                 the nth segment) has fields
%                   trialId       -- identifier of trial from which 
%                                    segment was taken
%                   segId         -- segment identifier within trial
%                   trialType     -- trial type index (Optional)
%                   fs            -- sampling frequency of ECoG data
%                   T (1 x 1)     -- number of timesteps in segment
%                   y (yDim x T)  -- neural data
% cutTrialGuess	- data structure has fields
%                   state         -- 'stateGuess' for resegmented trials
%                   mixComp       -- 'mixCompGuess' for resegmented trials
%
% OPTIONAL ARGUMENTS:
%
% method        - method for extracting neural trajectories:
%                 'mhmfa', 'mhmm', 'hmfa' (default), 'mfa', 'hmm', 'gmm'
% segLength     - length of segments to extract, in number of timesteps.
%                 If infinite, entire trials are extracted, i.e., no 
%                 resegmenting. (default: Inf)
% stateGuess  	- cell array of initial state guesses for time points
%                 of trials for MHMFA/MHMM/HMFA/HMM
% mixCompGuess 	- initial mixture component guesses for trials of
%                 MHMFA/MHMM (vector) or time points of trials for
%                 MFA/GMM (cell array)
%
% Code adapted from cutTrials.m by Byron Yu.
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  method                          = 'hmfa';
  segLength                       = Inf;
  stateGuess                      = [];
  mixCompGuess                    = [];
  assignopts(who, varargin);

  if isinf(segLength)
    seqOut                        = seqIn;
    switch(method)
     case {'mhmfa', 'mhmm'}
       cutTrialGuess.state        = stateGuess;
       cutTrialGuess.mixComp      = mixCompGuess;
     case {'hmfa', 'hmm'}
       cutTrialGuess.state        = stateGuess;
     case {'mfa', 'gmm'}
       cutTrialGuess.mixComp      = mixCompGuess;
     otherwise
       error('Invalid method specification');
    end % switch(method)
    return
  end % if isinf(segLength)

  seqOut                          = [];
  switch(method)
   case {'mhmfa', 'mhmm'}
     if iscell(stateGuess)
       cutTrialGuess.state        = {};
     elseif isnumeric(stateGuess) % && isempty(stateGuess)
       cutTrialGuess.state        = [];
     end
     cutTrialGuess.mixComp        = [];
   case {'hmfa', 'hmm'}
     if iscell(stateGuess)
       cutTrialGuess.state        = {};
     elseif isnumeric(stateGuess) % && isempty(stateGuess)
       cutTrialGuess.state        = [];
     end
   case {'mfa', 'gmm'}
     if iscell(mixCompGuess)
       cutTrialGuess.mixComp     	= {};
     elseif isnumeric(mixCompGuess) % && isempty(mixCompGuess)
       cutTrialGuess.mixComp      = [];
     end
   otherwise
     error('Invalid method specification');
  end % switch(method)
  for n=1:numel(seqIn)
    T                             = seqIn(n).T;
    
    % Skip trials that are shorter than segLength
    if (T < segLength)
      fprintf(['Warning: trialId %4d shorter than one ',...
               'segLength...skipping\n'], seqIn(n).trialId); 
      continue
    end % if (T < segLength)
    
    numSeg                        = ceil(T/segLength);
    if (numSeg == 1)
      cumOL                       = 0;
    else
      totalOL                     = (segLength*numSeg) - T;
      probs                       = ones(1,numSeg-1)/(numSeg-1);
      % mnrnd is very sensitive to sum(probs) being even slightly
      % away from 1 due to floating point round-off.
      probs(end)                  = 1-sum(probs(1:end-1));
      randOL                      = mnrnd(totalOL, probs);
      cumOL                       = [0 cumsum(randOL)];
    end

    seg                           = rmfield(seqIn(n),{'T','y'});
    seg.T                         = segLength;
    
    for s=1:numSeg
      tStart                      = -cumOL(s) + segLength * (s-1) + 1;
      seg.segId                   = s;
      seg.y                       =...
       seqIn(n).y(:, tStart:(tStart+segLength-1));
      seqOut                      = [seqOut seg];

      switch(method)
       case {'mhmfa', 'mhmm'}
         if ~isempty(stateGuess)
           cutTrialGuess.state    =...
            [cutTrialGuess.state,...
             stateGuess{n}(:,tStart:(tStart+segLength-1))];
         end % if ~isempty(stateGuess)
         if ~isempty(mixCompGuess)
           cutTrialGuess.mixComp	=...
            [cutTrialGuess.mixComp, mixCompGuess(n)];
         end % if ~isempty(mixCompGuess)
       case {'hmfa', 'hmm'}
         if ~isempty(stateGuess)
           cutTrialGuess.state    =...
            [cutTrialGuess.state,...
             stateGuess{n}(:,tStart:(tStart+segLength-1))];
         end % if ~isempty(stateGuess)
       case {'mfa', 'gmm'}
         if ~isempty(mixCompGuess)
           cutTrialGuess.mixComp  =...
            [cutTrialGuess.mixComp,...
             mixCompGuess{n}(:,tStart:(tStart+segLength-1))];
         end % if ~isempty(mixCompGuess)
       otherwise
         error('Invalid method specification');
      end % switch(method)
    end % for s=1:numSeg
  end % for n=1:numel(seqIn)