function seq = getSeqRMS(dat, binWidth)
%
% seq = getSeqRMS(dat, binWidth)
%
% Converts ECoG voltage to RMS power.
%
% INPUTS:
%
% dat         - structure whose n-th entry (corresponding to the n-th
%               experimental trial) has fields
%                 trialId       -- unique trial identifier
%                 trialType     -- trial type index (OPTIONAL)
%                 fs            -- sampling frequency of ECoG data
%                 ECoG          -- matrix of voltage activity across all
%                                  electrodes. Each row corresponds to an
%                                  electrode. Each column corresponds to a
%                                  (1/fs) sec timestep.
% binWidth    - ECoG window width in sec
%
% OUTPUTS:
%
% seq         - data structure, whose n-th entry (corresponding to
%               the n-th experimental trial) has fields
%                 trialId       -- unique trial identifier
%                 trialType     -- trial type index (OPTIONAL)
%                 fs            -- sampling frequency of ECoG data
%                 T (1 x 1)     -- number of timesteps
%                 y (yDim x T)  -- neural data
%
% Code adapted from getSeq.m by Byron Yu.
%
% @ 2015 Akinyinka Omigbodun    aomigbod@ucsd.edu

  rms                   = @(x,d) sqrt(mean(x.^2,d));

  if isfield(dat, 'T')
    seq                 = dat;
    return
  else % if ~isfield(dat, 'T')
    if isfield(dat, 'trialType')
      seq             	= selectfield(dat, {'trialId', 'trialType', 'fs'});
    else % if ~isfield(dat, 'trialType')
      seq             	= selectfield(dat, {'trialId', 'fs'});
    end % if isfield(dat, 'trialType')
  end

  for n=1:numel(dat)
    yDim                = size(dat(n).ECoG, 1);
    nSamples            = size(dat(n).ECoG, 2);
    nSamplesStep        = max(floor(dat(n).fs * binWidth), 1);
    T                   = ceil((nSamples - nSamplesStep + 1)/nSamplesStep);

    seq(n).T            = T;
    
    if (binWidth <= seq(n).fs^-1)
      seq(n).y          = abs(dat(n).ECoG);
    else % if (binWidth > seq(n).fs^-1)
      iStart            = 1;
      iEnd              = nSamplesStep;
      seq(n).y          = nan(yDim, T);
      for t=1:T
        seq(n).y(:,t)   = rms(dat(n).ECoG(:, iStart:iEnd), 2);

        iStart          = iStart + nSamplesStep;
        iEnd            = iEnd + nSamplesStep;
      end % for t=1:T
    end
  end % for n=1:numel(dat)
  
  % Remove trials that are shorter than one bin width
  if ~isempty(seq)
    trialsToKeep        = ([seq.T] > 0);
    seq                 = seq(trialsToKeep);
  end
end