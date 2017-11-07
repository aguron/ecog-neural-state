function result = processfilename(str)
%PROCESSFILENAME extracts method, nMixComp, nStates, xDim, faType, covType,
%   sharedCov, cvf, and ncvf from a variety of filename formats
%
% INPUTS:
%
% str      - filename string to parse
%
% OUTPUTS:
%
% result   - structure with fields
%              method     -- method for model fitting, inference, and/or
%                            prediction
%              nMixComp   -- number of mixture components
%              nStates    -- number of states
%              xDim       -- state dimensionality
%              faType     -- factor analyzers specification
%              covType    -- covariance type
%              sharedCov  -- covariance tied or untied
%              cvf        -- cross-validation fold
%              ncvf       -- number of cross-validation folds
%
% Code adapted from parseFilename.m by Byron Yu.
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  str                           = strtok(str, '.');
  undi                          = find(str == '_');
  nUnd                          = numel(undi);
  if isempty(undi)
    error('Invalid filename format');
  end

  [result.method, ~, errmsg]    = sscanf(str(1:undi(1)-1), '%s');
  if errmsg
    error('Invalid filename format');
  end

  i                             = 1;
  if ismember(result.method, {'mhmfa', 'mhmm'})
    [result.nMixComp, ~, errmsg]= sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'nMixComp%d');
    if errmsg
      error('Invalid filename format');
    end
    i                           = i + 1;
  end
  switch(result.method) % OUTER
    case {'pca', 'ppca', 'fa', 'gpfa',...
          'mfa', 'hmfa', 'mhmfa'}
      if (nUnd > 1)
        [result.xDim, ~, errmsg]= sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'xDim%d');
        if errmsg
          error('Invalid filename format');
        end
        p                       = i;
        while (true)
          p                     = p + 1;
          if (p < nUnd)
            [sp, ~, errmsg]   	= sscanf(str(undi(p)+1:undi(p+1)-1), '%d');
            if errmsg
              i                	= p - 1;
              break
            end
          else
            i                   = p - 1;
            break
          end
          result.xDim(end+1)    = sp;
        end
      else % if (nUnd == 1)
        [result.xDim, ~, errmsg]= sscanf(str(undi(i)+1:end), 'xDim%d');
        if errmsg
          error('Invalid filename format');
        end
      end

      switch(result.method) % INNER
        case {'pca', 'ppca', 'fa', 'gpfa'}
          % do nothing
        case {'mfa', 'hmfa', 'mhmfa'}
          i                     = i + 1;
          if isequal(result.method, 'mfa')
            [result.nMixComp, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'nMixComp%d');
          elseif ismember(result.method, {'hmfa', 'mhmfa'})
            [result.nStates, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'nStates%d');
          end
          if errmsg
            error('Invalid filename format');
          end

          p                    	= i;
          while (true)
           p                  	= p + 1;
           if (p < nUnd)
             [sp, ~, errmsg]   	= sscanf(str(undi(p)+1:undi(p+1)-1), '%d');
             if errmsg
               i               	= p - 1;
               break
             end
           else
             i                	= p - 1;
             break
           end
           result.nStates(end+1)= sp;
          end

          i                    	= i + 1;
          result.faType         = {str(undi(i)+4:undi(i)+6)};
          faTypeSpec            = 'tu'; % t - tied; u - untied
          result.faType{1}     	= double(result.faType{1}==faTypeSpec(2));

          if (i < nUnd)
            p                  	= i;
            while (true)
              p                	= p + 1;

              result.faType{end+1}...
                                = str(undi(p)+1:undi(p)+3);
              result.faType{end}=double(result.faType{end}==faTypeSpec(2));
              if (p == nUnd)
                if (numel(str(undi(p)+1:end)) > 3)
                  result.faType(end)...
                                = [];
                  i             = p - 1;
                else
                  i             = p;
                end
                break
              end % if (p == nUnd)
            end
          end
        otherwise
          error('Invalid filename format');
      end % switch(result.method) % INNER
    case {'gmm', 'hmm', 'mhmm'}
      if isequal(result.method, 'gmm')
        [result.nMixComp, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'nMixComp%d');
      elseif ismember(result.method, {'hmm', 'mhmm'})
        [result.nStates, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'nStates%d');
      end
      if errmsg
        error('Invalid filename format');
      end
      i                         = i + 1;
      if (i < nUnd)
        result.covType          = str(undi(i)+1:undi(i+1)-1);
        i                       = i + 1;
        if (i < nUnd)
          result.sharedCov      = str(undi(i)+1:undi(i+1)-1);
          if isequal(result.sharedCov,'tied')
            result.sharedCov    = true;
          else
            error('Invalid filename format');
          end
        else
          result.sharedCov      = str(undi(i)+1:end);
          if isequal(result.sharedCov,'tied')
            result.sharedCov    = true;
          else
            result.sharedCov    = false;
            i                   = i - 1;
          end
        end
      else
        result.covType          = str(undi(i)+1:end);
        result.sharedCov        = false;
      end
      if ~ismember(result.covType, {'diagonal','full'})
        error('Invalid filename format');
      end
    otherwise
      error('Invalid filename format');
  end % switch(result.method) % OUTER

  if (i < nUnd)
    i                           = i + 1;
    [A, ~, errmsg]              = sscanf(str(undi(i)+1:end), 'cv%dof%d');

    if errmsg
      error('Invalid filename format');
    end
    
    if (numel(A) == 1)
      result.cvf                = A(1);
      result.ncvf               = NaN;
    elseif (numel(A) == 2)
      result.cvf               	= A(1);
      result.ncvf              	= A(2);      
    end
  end % if (i < nUnd)

  if ~isfield(result, 'nMixComp')
    result.nMixComp           	= nan;
  end
  if ~isfield(result, 'nStates')
    result.nStates              = nan;
  end
  if ~isfield(result, 'xDim')
    result.xDim                 = nan;
  end
  if ~isfield(result, 'faType')
    result.faType               = nan(1,3);
  end
  if ~isfield(result, 'covType')
    result.covType             	= '';
  end
  if ~isfield(result, 'sharedCov')
    result.sharedCov            = nan;
  end
  if ~isfield(result, 'cvf')
    result.cvf                  = 0;
  end
  if ~isfield(result, 'ncvf')
    result.ncvf                	= 0;
  end

  if (numel(result.faType) == 1)
    result.faType               = cell2mat(result.faType);
  end
  
  result                        = orderfields(result);
end