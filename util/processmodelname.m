function result = processmodelname(str)
%PROCESSMODELNAME
%   
  undi                          = find(str == '_');
  nUnd                          = numel(undi);
  if isempty(undi)
    error('Invalid model name format');
  end % if isempty(undi)

  [result.method, ~, errmsg]    = sscanf(str(1:undi(1)-1), '%s');
  if errmsg
    error('Matching failure in format');
  end % if errmsg

  i                             = 1;
  if ismember(result.method, {'mhmfa', 'mhmm'})
    [result.nMixComp, ~, errmsg]= sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'nMixComp%d');
    if errmsg
      error('Matching failure in format');
    end % if errmsg
    i                           = i + 1;
  end % if ismember(result.method, {'mhmfa', 'mhmm'})
  switch(result.method) % OUTER
    case {'pca', 'ppca', 'fa', 'gpfa',...
          'mfa', 'hmfa', 'phmfa', 'chmfa',...
          'mhmfa'}
      if (nUnd > 1)
        [result.xDim, ~, errmsg]= sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'xDim%d');
        if errmsg
          error('Matching failure in format');
        end % if errmsg
        p                       = i;
        while (true)
          p                     = p + 1;
          if (p < nUnd)
            [sp, ~, errmsg]   	= sscanf(str(undi(p)+1:undi(p+1)-1), '%d');
            if errmsg
              i                	= p - 1;
              break
            end % if errmsg
          else
            i                   = p - 1;
            break
          end
          result.xDim(end+1)    = sp;
        end % while (true)
      else % if (nUnd == 1)
        [result.xDim, ~, errmsg]= sscanf(str(undi(i)+1:end), 'xDim%d');
        if errmsg
          error('Matching failure in format');
        end % if errmsg
      end

      switch(result.method) % INNER
        case {'pca', 'ppca', 'fa'}
          i                     = i + 1;
          if (i < nUnd)
            [result.kernSD, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         'kernSD%g');
            if errmsg
              error('Matching failure in format');
            end % if errmsg
          elseif (i == nUnd)
            [result.kernSD, ~, errmsg]...
                                = sscanf(str(undi(i)+1:end), 'kernSD%g');
            if errmsg
              error('Matching failure in format');
            end % if errmsg
          end
        case 'gpfa'
          % do nothing
        case {'mfa', 'hmfa', 'phmfa', 'chmfa', 'mhmfa'}
          i                     = i + 1;
          if isequal(result.method, 'mfa')
            label               = 'nMixComp';
          elseif ismember(result.method,...
                          {'hmfa', 'phmfa', 'chmfa', 'mhmfa'})
            label               = 'nStates';
          end
          [result.nStates, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         [label,'%d']);
          if errmsg
            error('Matching failure in format');
          end % if errmsg

          p                    	= i;
          while (true)
           p                  	= p + 1;
           if (p < nUnd)
             [sp, ~, errmsg]   	= sscanf(str(undi(p)+1:undi(p+1)-1), '%d');
             if errmsg
               i               	= p - 1;
               break
             end % if errmsg
           else
             i                	= p - 1;
             break
           end
           result.nStates(end+1)= sp;
          end % while (true)

          i                    	= i + 1;
          result.faFtrType     	= str(undi(i)+1:undi(i)+3);
          result.faFtrType     	= eq(result.faFtrType, 'mfv');
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
            end % while (true)
          end % if (i < nUnd)
        otherwise
          error('Invalid method specification');
      end % switch(result.method) % INNER
    case {'gmm', 'hmm', 'mhmm'}
      if isequal(result.method, 'gmm')
        label                   = 'nMixComp';
      elseif ismember(result.method, {'hmm', 'mhmm'})
        label                   = 'nStates';
      end
      [result.nStates, ~, errmsg]...
                                = sscanf(str(undi(i)+1:undi(i+1)-1),...
                                         [label,'%d']);
     if errmsg
       error('Matching failure in format');
     end % if errmsg
     i                          = i + 1;
     if (i < nUnd)
       result.covType          	= str(undi(i)+1:undi(i+1)-1);
     else
       result.covType          	= str(undi(i)+1:end);
     end
    otherwise
      error('Invalid method specification');
  end % switch(result.method) % OUTER

  if ~isfield(result, 'nMixComp')
    result.nMixComp           	= nan;
  end % if ~isfield(result, 'nMixComp')
  if ~isfield(result, 'xDim')
    result.xDim                 = nan;
  end % if ~isfield(result, xDim)
  if ~isfield(result, 'nStates')
    result.nStates              = nan;
  end % if ~isfield(result, nStates)
  if ~isfield(result, 'faFtrType')
    result.faFtrType            = nan(1,3);
  end % if ~isfield(result, 'faFtrType')
  if ~isfield(result, 'faType')
    result.faType               = nan(1,3);
  end % if ~isfield(result, faType)
  if ~isfield(result, 'covType')
    result.covType             	= '';
  end % if ~isfield(result, covType)
  if ~isfield(result, 'kernSD')
    result.kernSD               = 0;
  end % if ~isfield(result, 'kernSD')

  if (numel(result.faType) == 1)
    result.faType               = cell2mat(result.faType);
  end % if (numel(result.faType) == 1)
  
  result                        = orderfields(result);
end