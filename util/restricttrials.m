function seq = restricttrials(seq, trialType)
%RESTRICTTRIALS
%   
  seq	= seq(ismember([seq.trialType],trialType));
end

