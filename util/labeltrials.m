function [seqInfo, labelInfo] = labeltrials(seq, datInfo)
%LABELTRIALS
%
% INPUTS:
% 
% datInfo   -   matrix with trial labels in rows ordered according to the
%               lexicographic reordering of seq based on trialId field, OR
%               struct with fields
%               field,
%                 trialId       -- unique trial identifier
%                 trialType     -- index from 1 to the number of
%                                  trial types inclusive
%   
% seq       -   struct of trials without trialType field (not necessarily
%               with the same ordering as datInfo)
%
% OUTPUTS:
%
% seqInfo   -   vector of trialType labels for seq trials, OR
%               struct with seq trials and trialType field
%
% labelInfo	-   cell array of unique datInfo labels (only for datInfo
%               matrix format)


  if isnumeric(datInfo)
    if any(~isint(datInfo(:))) || any(datInfo(:) <= 0)
      error('All entries in datInfo must be positive integers');
    end

    nTrials                               = size(datInfo,1);
    if (numel(seq) ~= nTrials)
      error(['Number of output trajectories is not',...
             ' equal to the number of input trials']);
    end % if (numel(seq) ~= nTrials)

    [~, idxList]                          = sort({seq(:).trialId});
    [~, idxList]                          = sort(idxList);

    seqInfo                               = datInfo(idxList,:);

    nAttributes                           = size(datInfo,2);
    nStimuli                              = zeros(1, nAttributes);
    for a=1:nAttributes
      nStimuli(a)                         = length(unique(datInfo(:,a)));
    end % for a=1:nAttributes


    seqInfo                               =...
      1 + sum(bsxfun(@dotprod,...
                     seqInfo-1,...
                     [1 cumprod(nStimuli(1:end-1))]), 2);
    categories                            = unique(seqInfo);
    nCategories                           = length(categories);
    for c=1:nCategories
      seqInfo(seqInfo == categories(c))   = c;
    end % for i=1:nCategories

    if (nargout == 2)
      labelInfo                           =...
        mat2cell(fliplr(unique(fliplr(datInfo),'rows')),...
                 ones(1,nCategories), nAttributes);
    end % if (nargout == 2)
  elseif isstruct(datInfo)
    seqInfo                               = seq;
    trialId                               = {datInfo.trialId};
    for n=1:numel(seqInfo)
      seqInfo(n).trialType                =...
        datInfo(ismember(trialId,seqInfo(n).trialId)).trialType;
    end % for n=1:numel(seqInfo)
    if (nargout == 2)
      labelInfo                          	= [];
    end % if (nargout == 2)
  else
    error('datInfo must be a matrix or a struct');
  end
end