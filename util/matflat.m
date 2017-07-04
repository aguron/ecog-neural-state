function B = matflat(A, arg)
%MATFLAT
%   Convert matrices to column vectors

% INPUT - matrix, cell array or struct array of matrices

  if iscell(A)
    B               = cell(size(A));
    for m=1:numel(A)
      B{m}          = A{m}(:);
    end % for m=1:numel(A)
  elseif isstruct(A)
    B               = A;
    for m=1:numel(A)
      B(m).(arg)    = A(m).(arg)(:);
    end % for m=1:numel(A)
  else
    if (nargin == 1)
      B            	= A(:);
    elseif (nargin == 2)
      if (arg(1) == arg(2)) || (arg(1) > ndims(A))
        B           = A;
        return
      end % if (arg(1) == arg(2)) || (arg(1) > ndims(A))
      sizeDims      = size(A);
      idx          	= cell(ndims(A),1);
      for d=[1:arg(1)-1 arg(1)+1:ndims(A)]
        idx{d}      = 1:sizeDims(d);
      end % for d=[1:arg(1)-1 arg(1)+1:ndims(A)]
      
      C             = cell(sizeDims(arg(1)),1);
      for i=1:sizeDims(arg(1))
        idx{arg(1)}	= i;
        C{i}        = A(idx{:});
      end % for i=1:sizeDims(arg(1))

      B             = cat(arg(2),C{:});
    end
  end % if iscell(A)
end

