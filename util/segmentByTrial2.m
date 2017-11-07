function seq = segmentByTrial2(seq, X, fn)
%
% seq = segmentByTrial2(seq, X, fn)
%
% Segment and store data by trial.
%  
% INPUT:
%
% seq  	- data structure that has field T, the number of timesteps
% X    	- data to be segmented
%                (any dimensionality x total number of timesteps)
% fn    - new field name of seq where segments of X are stored
%
% OUTPUT:
%
% seq   - data structure with new field specified as a string in fn
%
% Code adapted from segmentByTrial.m by Byron Yu.
%
% @ 2015 Akinyinka Omigbodun    aomigbod@ucsd.edu

  if sum([seq.T]) ~= size(X, 2)
    fprintf('Error: size of X incorrect.\n');
  end
  
  ctr           = 0;
  for n=1:numel(seq)
    T           = seq(n).T;
    idx         = (ctr+1) : (ctr+T);
    seq(n).(fn) = X(:, idx, :);
    
    ctr         = ctr + T;
  end
