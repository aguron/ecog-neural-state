function B = ndarrayflatten(A, arg)
%NDARRAYFLATTEN converts multidimensional arrays (including those within
%   cell and struct arrays) to column vectors. It can also 'flatten' a
%   selected dimension of a multidimensional array into another specified
%   dimension of that array.
%
% INPUTS:
%
% A   - multidimensional array, cell array or struct array of
%       multidimensional arrays
% arg - EITHER:
%              specified field            ] when A is a struct array
%           OR:
%              vector of length 2 with the]
%              first entry specifying the ]
%              dimension to be 'flattened'] when A is a multidimensional
%              and the second entry       ] array (with a dimension to be
%              specifying a dimension that] 'flattened')
%              the first is 'flattened'   ]
%              into                       ]
%
% OUTPUTS:
%
% B   - column vector, cell array or struct array of column vectors, or
%       modified multidimensional array
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

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
  elseif isnumeric(A)
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
      end
      
      C             = cell(sizeDims(arg(1)),1);
      for i=1:sizeDims(arg(1))
        idx{arg(1)}	= i;
        C{i}        = A(idx{:});
      end

      B             = cat(arg(2),C{:});
    end
  end
end

