function [idx,centroid] = kmeansclszthr(X,k,Thr,nRuns)
%KMEANSCLSZTHR K-means clustering while attempting to satisfy a minimum
%   cluster size threshold
%
% INPUTS:
%
% X         - matrix with points for clustering in rows
%	k         - number of clusters
% Thr       - minimum cluster size threshold
% nRuns     - maximum number of attempts to satisfy minimum cluster size
%             threshold
%
% OUTPUTS:
%
% idx       - cluster assignments
% centroid	- matrix with cluster centroids in rows
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  for i=1:nRuns
    init                          	= true;
    fprintf('K-means run %3d of %d\n', i, nRuns);
    
    [temp, temp2]                   = kmeans(X,k);
    temp3                           = histc(temp,1:k);
    if (i == 1) || (min(nSamples) < min(temp3))
      idx                           = temp;
      centroid                      = temp2;
      nSamples                      = temp3;
    end
    for j=1:k
      if (nSamples(j) < Thr)
        init                        = false;
      end % if (nSamples(j) < yDim)
    end % for j=1:k
    if (init)
      break
    else
      continue
    end
  end % for i=1:nRuns
  if ~(init)
    disp('Minimum cluster size threshold not yet reached');
    programcontrol
  end % if ~(init)
end