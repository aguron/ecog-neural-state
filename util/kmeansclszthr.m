function [idx,centroid] = kmeansclszthr(X,k,Thr,nRuns)
%KMEANSCLSZTHR

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