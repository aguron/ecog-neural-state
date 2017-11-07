function modelclasscomparisons(hyptest,beforeBHP,afterBHP,displayInfo,...
                               varargin)
%MODELCLASSCOMPARISONS performs a multiple comparisons analysis between
%   two model classes with the Benjamini-Hochberg procedure
%   (http://www.biostathandbook.com/multiplecomparisons.html)
%
% INPUTS:
%
% hyptest         - structure array (with D dimensions) with
%                   model comparison hypothesis test information
% beforeBHP       - cell array (of length D) specifying model space
%                   selection before the Benjamini-Hochberg procedure
%                   ([] is used to specify that a model space dimension
%                    is to remain untouched)
% afterBHP        - cell array (of length D) specifying model space
%                   selection after the Benjamini-Hochberg procedure
%                   ([] is used to specify that a model space dimension
%                    is to remain untouched)
% displayInfo     - cell array (of length D) with information to be
%                   displayed on the model space selection. Each entry is
%                   itself a cell array in the following format:
%                   {<Model Space Dimension Name String>,
%                    <Model Space Dimension Selected Values Format Specifier>,
%                    <Vector of Model Space Dimension Selected Values>}
%
% OPTIONAL ARGUMENTS:
%
% metric          - 'predError', 'mInf', or 'accuracy'
%                   (default: 'predError')
% FDR             - Benjamini-Hochberg procedure false discovery rate
%                   (default: 0.05)
% modelClassNames	- cell of length 2 with names of model classes in the
%                   multiple comparisons analysis (default: {'A','B'})
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  metric                    = 'mInf';
  FDR                       = 0.05;
  modelClassNames           = {'A','B'};
  assignopts(who, varargin);

  if ~ismember(metric, {'predError', 'mInf', 'accuracy'})
    error('Invalid metric specification');
  end % if ~ismember(metric, {'predError', 'mInf', 'accuracy'})

  nDim                      = ndims(hyptest);
  if ~isequal(nDim, numel(beforeBHP), numel(afterBHP), numel(displayInfo))
    error(['The number of dimensions in hyptest, ',...
           'beforeBHP, and afterBHP must be the same'])
  end

  % Model comparisons
  missing                   = 0;
  inconclusive             	= 0;
  percent                 	= zeros(1,3);

  % Model selection before Benjamini-Hochberg procedure
  for d=1:nDim
   if isempty(beforeBHP{d})
    continue
   end
   dim                      = cell(1,nDim);
   [dim{:}]                 = deal(':,');
   dim{end}                 = ':';
   dim{d}                   = 1:size(hyptest,d);
   dim{d}                   = setdiff(dim{d},beforeBHP{d});
   dim{d}                   = ['[', num2str(dim{d}), ']'];
   if (d < nDim)
    dim{d}                  = [dim{d}, ','];
   end % if (d < nDim)
   eval(['hyptest(', dim{:}, ') = [];'])
  end % for d=1:nDim
  
  dim                       = size(hyptest);
  % Benjamini-Hochberg procedure for multiple comparisons:
  %  http://www.biostathandbook.com/multiplecomparisons.html
  pVal                     	= nan(prod(dim),3);
  i                         = 0;
  for d=counter(dim)'
   i                        = i + 1;
   d                      	= num2cell(d{1});
   if isempty(hyptest(d{:}).(metric))
    pVal(i,:)               = [];
    i                       = i - 1;
    continue
   end % if isempty(hyptest(d{:}).(metric))
   pVal(i,1)                = hyptest(d{:}).(metric).rslt(1).p;
   pVal(i,2)               	= hyptest(d{:}).(metric).rslt(2).p;
   pVal(i,3)                = hyptest(d{:}).(metric).rslt(3).p;
  end % for d=counter(dim)'
  FDRVal                    = nan(size(pVal));
  FDRVal(:,1)               = mafdr(pVal(:,1),'BHFDR',true);
  FDRVal(:,2)               = mafdr(pVal(:,2),'BHFDR',true);
  FDRVal(:,3)               = mafdr(pVal(:,3),'BHFDR',true);
  i                         = 0;
  for d=counter(dim)'
   i                        = i + 1;
   d                      	= num2cell(d{1});
   if isempty(hyptest(d{:}).(metric))
    i                       = i - 1;
    continue
   end % if isempty(hyptest(d{:}).(metric))
   hyptest(d{:}).(metric).rslt(1).p...
                            = FDRVal(i,1);
   hyptest(d{:}).(metric).rslt(2).p...
                            = FDRVal(i,2);
   hyptest(d{:}).(metric).rslt(3).p...
                            = FDRVal(i,3);
  end % for d=counter(dim)'
  
  % Model selection after Benjamini-Hochberg procedure
  for d=1:nDim
   if isempty(afterBHP{d})
    continue
   end
   dim                      = cell(1,nDim);
   [dim{:}]                 = deal(':,');
   dim{end}                 = ':';
   dim{d}                   = 1:size(hyptest,d);
   dim{d}                   = setdiff(dim{d},afterBHP{d});
   dim{d}                   = ['[', num2str(dim{d}), ']'];
   if (d < nDim)
    dim{d}                  = [dim{d}, ','];
   end % if (d < nDim)
   eval(['hyptest(', dim{:}, ') = [];'])
  end % for d=1:nDim

  dim                     	= size(hyptest);
  for d=counter(dim)'
   d                      	= num2cell(d{1});
   if isempty(hyptest(d{:}).(metric))
    missing                	= missing + 1;
    continue
   end % if isempty(hyptest(d{:}).(metric))
   crit(1)                  = hyptest(d{:}).(metric).rslt(1).p < FDR;
   crit(2)                	= hyptest(d{:}).(metric).rslt(2).p < FDR/2;
   crit(3)                	= hyptest(d{:}).(metric).rslt(3).p < FDR/2;
   if (crit(1))
    % statistically significant
    if (crit(2))
     % statistically significantly smaller
     percent(1)           	= percent(1) + 1;
    elseif (crit(3))
     % statistically significantly greater
     percent(3)            	= percent(3) + 1;
    else
     inconclusive          	= inconclusive + 1;
    end
   else
    % statistically insignificant
    percent(2)            	= percent(2) + 1;
   end
  end % for d=counter(dim)'

  
  % Display model class comparison results
  nHypTest                 	= prod(dim);
  percent                 	= 100*percent/(nHypTest - missing);
  fprintf('********************\n');
  switch(metric)
   case 'mInf'
    fprintf('Mutual Information\n');
   case 'predError'
    fprintf('Prediction Error\n');
   case 'accuracy'
    fprintf('Percentage Accuracy\n');
   otherwise
    error('Invalid metric specification');
  end
  fprintf('********************\n');
  
  for d=1:nDim
   fprintf('%s\n', displayInfo{d}{1});
   fprintf([displayInfo{d}{2},'\n'], displayInfo{d}{3});
  end % for d=1:nDim
  fprintf('%s versus %s\n', modelClassNames{:});

  fprintf('Number of comparisons: %d\n', nHypTest - missing);
  fprintf('metric is statistically significantly greater for %s: %s (%s%%)\n',...
          modelClassNames{2},...
          dadp(percent(1)*(nHypTest - missing)/100),...
          dadp(percent(1),2));
  fprintf('statistically insignificant: %s (%s%%)\n',...
          dadp(percent(2)*(nHypTest - missing)/100), dadp(percent(2),2));
  fprintf('metric is statistically significantly greater for %s: %s (%s%%)\n',...
          modelClassNames{1},...
          dadp(percent(3)*(nHypTest - missing)/100),...
          dadp(percent(3),2));
  fprintf('Inconclusive comparisons: %d\n', inconclusive);
  fprintf('Missing comparisons: %d\n', missing);
  fprintf('\n');
end

