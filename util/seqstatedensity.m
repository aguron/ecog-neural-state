function seqOut = seqstatedensity(seq, params, varargin)
%SEQSTATEDENSITY
% 
  seqSpec                       = 'state';
  nSegments                     = 1;
  assignopts(who, varargin);

  nSeq                          = length(seq);
  seqOut                        = struct('stateDensity', cell(1,nSeq));
  
  switch(seqSpec)
    case {'state', 'cluster'}
      paramsSpec                = 'nStates';
    case 'mixComp'
      paramsSpec                = 'nMixComp';
    otherwise
      error('Invalid seqSpec');
  end % switch(seqSpec)

  nStates                       = params.(paramsSpec);
  if (numel(nStates) > 1)
    states                      = counter(nStates);
  end % if (numel(nStates) > 1)
  nSegEqLen                     = false;
  for i=1:nSeq
    [~, seqLen]                 = size(seq(i).(seqSpec));
    if (numel(nStates) > 1)
      temp                      = zeros(1, seqLen);
      for s=1:numel(states)
        temp(all(bsxfun(@eq,seq(i).(seqSpec),states{s}')))...
                                = s;
      end % for s=1:numel(states)
      seq(i).(seqSpec)          = temp;
    end % if (numel(nStates) > 1)
    
    if isempty(nSegments) || nSegEqLen
      nSegments                 = seqLen;
      nSegEqLen                 = true;
    end % if isempty(nSegments)
    c_seqLen_j                  = ceil(seqLen/nSegments);
    f_seqLen_j                  = floor(seqLen/nSegments);
    nC                          = seqLen - (f_seqLen_j * nSegments);
    seqOut(i).stateDensity      = zeros(1, prod(nStates)*nSegments);
    bIdx1                       = 1;
    eIdx1                       = prod(nStates);
    bIdx2                       = 1;
    if (nC)
      eIdx2                     = c_seqLen_j;
    else
      eIdx2                     = f_seqLen_j;
    end % if (nC)

    for j=1:nSegments
      if (j < nSegments)
        seqOut(i).stateDensity(bIdx1:eIdx1)...
                                = histc(seq(i).(seqSpec)(bIdx2:eIdx2),...
                                        1:prod(nStates));
      else % (j == nSegments)
        seqOut(i).stateDensity(bIdx1:eIdx1)...
                                = histc(seq(i).(seqSpec)(bIdx2:end),...
                                        1:prod(nStates));
      end % if (j < nSegments)
      seqOut(i).stateDensity(bIdx1:eIdx1)...
                                = normalize(seqOut(i).stateDensity(bIdx1:eIdx1));
      if (j < nSegments)
        bIdx1                   = eIdx1 + 1;
        eIdx1                   = eIdx1 + prod(nStates);
        if (j < nC)
          bIdx2                 = eIdx2 + 1;
          eIdx2                 = eIdx2 + c_seqLen_j;
        else % if (nC <= j < nSegments)
          bIdx2                 = eIdx2 + 1;
          eIdx2                 = eIdx2 + f_seqLen_j;
        end % if (j < nC)
      end % if (j < nSegments)
    end % for j=1:nSegments
  end % for i=1:nSeq
end