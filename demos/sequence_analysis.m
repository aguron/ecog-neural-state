if isdir(resultsDir)
  fprintf('Using existing directory %s...\n', resultsDir);
else % if ~isdir(resultsDir)
  fprintf('Making directory %s...\n', resultsDir);
  mkdir(resultsDir);
end

cd(resultsDir)

% Number of timesteps in a segment
segLength               = Inf;

% Other specifications
emMaxIters              = 250;

% Select number of cross-validation folds
numFolds                = 4;

% Parallel computing
plComputing             = true;
if (plComputing) && ~exist('pl2','var')
  pl                    = parcluster;
  pl.NumWorkers         = 12;
  pl2                   = parpool(pl);
end

patient                 =...
 {'NY400', 'NY441', 'NY451', 'NY453', 'NY455', 'NY468'};
for pa=pa_select
 if (~synthetic)
  p                    	= patient{pa};
 end
 for iB=binWidth_select
  % Specify directory where models will be saved or loaded from
  runIdx                = pa;
  % runDir              = sprintf('%s/mat_results/run%03d/binWidth_%g',
  %                               resultsDir, runIdx, binWidth(iB));
  for method=methodList
   if ismember(method{1}, {'gmm', 'hmm'})
     xDimList         = zeros();
     covTypeList      = covTypeRange;
     faTypeList       = {[]};
     sharedCovList    = sharedCovRange;
   elseif ismember(method{1}, {'mfa', 'hmfa'})
     xDimList         = xDimRange;
     covTypeList      = {''};
     faTypeList       = faTypeRange;
     sharedCovList    = {[]};
   end
   nStatesList        = nStatesRange;
   for xDim=xDimList
    for nStates=nStatesList
     for faType=faTypeList
      for CovType=covTypeList
       for SharedCov=sharedCovList
        if (~synthetic)
         args         = {dat.(p)};
        else % if (synthetic)
         args         = {seq};
        end
        switch(method{1})
         case 'gmm'
          args        = [args, 'nMixComp', nStates,...
                               'Start', 'randSample',...
                               'CovType', CovType{1},...
                               'SharedCov', SharedCov{1},...
                               'Options', statset('Display','iter',...
                                                  'TolFun',tolGMM,...
                                                  'MaxIter',emMaxIters)];
         case 'hmm'
          args        = [args, 'nStates', nStates,...
                               'tolHMM', tolHMM,...
                               'CovType', CovType{1},...
                               'SharedCov', SharedCov{1},...
                               'Options', statset('Display','iter',...
                                                  'TolFun',tolHMM,...
                                                  'MaxIter',emMaxIters)];
         case 'mfa'
          args        = [args, 'nMixComp', nStates,...
                               'tolMFA', tolMFA,...
                               'faType', faType{1},...
                               'xDim', xDim,...
                               'Options', statset('Display','iter',...
                                                  'TolFun',tolMFA,...
                                                  'MaxIter',emMaxIters)];
         case 'hmfa'
          args        = [args, 'nStates', nStates,...
                               'tolHMFA', tolHMFA,...
                               'faType', faType{1},...
                               'xDim', xDim,...
                               'Options',statset('Display','iter',...
                                                 'TolFun',tolHMFA,...
                                                 'MaxIter',emMaxIters)];
         otherwise
          error('Invalid specification of neural state method');
        end % switch(method{1})
        args          = [args, 'method', method{1},...
                               'segLength', segLength,...
                               'binWidth', binWidth(iB),...
                               'emMaxIters', emMaxIters,...
                               'numFolds', numFolds,...
                               'prediction', prediction,...
                               'fracTrainData', frac];

        try
         neuralstate(runIdx, args{:});
        catch err
         % Parallel computing bookkeeping
         if (plComputing)
          delete(gcp('nocreate'))
          clear pl pl2
         end % if (plComputing)
         displayerror(err)
         keyboard
        end % try
        
        % Two-stage spatiotemporal model optimization
        switch(method{1})
         case {'gmm', 'mfa'}
          sModel.dir        = '.';
          sModel.runIdx     = runIdx;
          sModel.binWidth   = binWidth(iB);
          sModel.method     = method{1};
          sModel.nMixComp   = nStates;
          sModel.numFolds   = numFolds;
          switch(method{1})
           case 'gmm'
            sModel.covType 	= CovType{1};
            sModel.sharedCov= SharedCov{1};

            stModel.dir    	= 'gmm-mm';
            stModel.method 	= 'hmm';
           case 'mfa'
            sModel.xDim     = xDim;
            sModel.faType  	= faType{1};
            
            stModel.dir    	= 'mfa-mm';
            stModel.method 	= 'hmfa';
          end % switch(method{1})

          stModelFlag     	= modelconverter(sModel, stModel);
          if (~stModelFlag)
            continue
          end

          if (~synthetic)
           args             = {dat.(p)};
          else % if (synthetic)
           args             = {seq};
          end
          switch(method{1})
           case 'gmm'
            args           	= [args,...
                               'nStates', nStates,...
                               'tolHMM', tolHMM,...
                               'CovType', CovType{1},...
                               'SharedCov', SharedCov{1},...
                               'Options', statset('Display','iter',...
                                                  'TolFun',tolHMM,...
                                                  'MaxIter',emMaxIters)];
           case 'mfa'
            args            = [args,...
                               'nStates', nStates,...
                               'tolHMFA', tolHMFA,...
                               'faType', faType{1},...
                               'xDim', xDim,...
                               'Options',statset('Display','iter',...
                                                 'TolFun',tolHMFA,...
                                                 'MaxIter',emMaxIters)];
          end % switch(method{1})
          args              = [args,...
                               'method', stModel.method,...
                               'segLength', segLength,...
                               'binWidth', binWidth(iB),...
                               'emMaxIters', emMaxIters,...
                               'numFolds', numFolds,...
                               'prediction', prediction,...
                               'predOverride', true,...
                               'inferOverride', true,...
                               'fracTrainData', frac];          
          cd(stModel.dir);
          try
           neuralstate(runIdx, args{:});
          catch err
           % Parallel computing bookkeeping
           if (plComputing)
            delete(gcp('nocreate'))
            clear pl pl2
           end % if (plComputing)
           displayerror(err)
           keyboard
          end % try
        end % switch(method{1})
       end % for SharedCov=sharedCovList
      end % for CovType=covTypeList
     end % for faType=faTypeList
    end % for nStates=nStatesList
   end % for xDim=xDimList
  end % for method=methodList
 end % for iB=binWidth_select
end % for pa=pa_select