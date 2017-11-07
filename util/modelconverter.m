function newModelFlag = modelconverter(oldModel, newModel, varargin)
%MODELCONVERTER creates an HMM/HMFA (spatiotemporal) model file from a
% GMM/MFA (spatial) model file
% 
% INPUTS:
%
% oldModel      - struct with GMM/MFA model file information with fields
%                   dir               -- model directory information
%                   runIdx            -- model directory information
%                                        (please see RUNIDX in NEURALSTATE)
%                   binWidth          -- ECoG window width in seconds
%                   method            -- method for model fitting,
%                                        inference, and/or prediction
%                   numFolds          -- number of cross-validation folds
%                   nMixComp          -- number of mixture components
%  EITHER (for GMM):
%                   covType           -- covariance type: 'full'
%                                        or 'diagonal'
%                   sharedCov         -- covariance tied (true)
%                                        or untied (false)
%      OR (for MFA):
%                   xDim              -- state dimensionality
%                   faType            -- factor analyzers specification
% newModel     	- struct with HMM/HMFA model file information with fields
%                   dir               -- model directory information
%                   method            -- method for model fitting,
%                                        inference, and/or prediction
%
% OUTPUTS:
%
% newModelFlag	- true only if a new HMM/HMFA model file is created
%
% OPTIONAL ARGUMENTS:
%
% overwrite     - true if current HMM/HMFA model file should be
%                 overwritten (default: false)
%
% @ 2017 Akinyinka Omigbodun    aomigbod@ucsd.edu

  overwrite                         = false;
  assignopts(who, varargin);
  
  newModelFlag                     	= true;

  oldDir                            =...
    sprintf('%s/mat_results/run%03d/binWidth_%g',...
            oldModel.dir, oldModel.runIdx, oldModel.binWidth);

  switch(oldModel.method)
    case 'gmm'
      oldfname                      =...
        sprintf('%s/%s_nMixComp%02d_%s',...
                oldDir, oldModel.method,...
                oldModel.nMixComp, oldModel.covType);
      if oldModel.sharedCov
        oldfname                   	= sprintf('%s_tied', oldfname);
      end % if oldModel.sharedCov
    case 'mfa'
      faTypeSpec                    = 'tu'; % t - tied; u - untied
      oldfname                      =...
        sprintf('%s/%s_xDim%02d_nMixComp%02d_MFV%c%c%c',...
                oldDir, oldModel.method, oldModel.xDim,...
                oldModel.nMixComp, faTypeSpec(oldModel.faType+1));
    otherwise
      error(['Invalid specification of neural ',...
             'trajectory extraction method']);
  end % switch(oldModel.method)

  newDir                            =...
    sprintf('%s/mat_results/run%03d/binWidth_%g',...
            newModel.dir, oldModel.runIdx, oldModel.binWidth);
  switch(newModel.method)
    case 'hmm'
      newfname                      =...
        sprintf('%s/%s_nStates%02d_%s',...
                newDir, newModel.method,...
                oldModel.nMixComp, oldModel.covType);
      if oldModel.sharedCov
        newfname                   	= sprintf('%s_tied', newfname);
      end % if oldModel.sharedCov
    case 'hmfa'
      faTypeSpec                    = 'tu'; % t - tied; u - untied
      newfname                      =...
        sprintf('%s/%s_xDim%02d_nStates%02d_MFV%c%c%c',...
                newDir, newModel.method, oldModel.xDim,...
                oldModel.nMixComp, faTypeSpec(oldModel.faType+1));
    otherwise
      error(['Invalid specification of neural ',...
             'trajectory extraction method']);
  end % switch(newModel.method)
  
  if exist([newfname, '.mat'],'file') && (~overwrite)
    newModelFlag                     = false;
    return
  end % if exist([newfname, '.mat'],'file') && (~overwrite)
  

  if isdir(newDir)
    fprintf('Using existing directory %s...\n', newDir);
  else % if ~isdir(newDir)
    fprintf('Making directory %s...\n', newDir);
    mkdir(newDir);
  end

  for cvf=0:oldModel.numFolds
    if (cvf > 0)
      model                        	=...
        load(sprintf('%s_cv%02dof%02d', oldfname, cvf, oldModel.numFolds));
    else
      model                        	= load(oldfname);
    end
    switch(newModel.method)
      case 'hmm'
        model.estParams.nStates     = model.obj.NComponents;

        model.estParams.piPrior    	= ones(1, model.estParams.nStates);
        model.estParams.transPrior	= ones(model.estParams.nStates);
        
        mixComp                     = cell(1,numel(model.seqTrain));
        if ~isempty(mixComp)
          [mixComp{:}]            	= deal(model.seqTrain.mixComp);
          newParams                	=...
            hmmestimate2(mixComp, mixComp,...
                         model.estParams.nStates, model.estParams.nStates,...
                         'Pseudostarts', model.estParams.piPrior,...
                         'Pseudotransitions', model.estParams.transPrior);
          model.estParams.pi        = newParams.ST;
          model.estParams.trans     = newParams.TR;
        else % if isempty(mixComp)
          model.estParams.pi        = [];
          model.estParams.trans     = [];
        end

        model.estParams.d           = model.obj.mu';
        if isequal(model.obj.CovType, 'full')
          model.estParams.R        	= model.obj.Sigma;
        elseif isequal(model.obj.CovType, 'diagonal')
          model.estParams.R        	= nddiag(model.obj.Sigma);
        end
        
        model.estParams.covType     = model.obj.CovType;
        model.estParams.sharedCov   = model.obj.SharedCov;
      case 'hmfa'
        model.estParams.nStates     = model.estParams.nMixComp;
        
        model.estParams.piPrior    	= ones(1, model.estParams.nStates);
        model.estParams.transPrior	= ones(model.estParams.nStates);
        
        mixComp                     = cell(1,numel(model.seqTrain));
        if ~isempty(mixComp)
          [mixComp{:}]            	= deal(model.seqTrain.mixComp);
          newParams                	=...
            hmmestimate2(mixComp, mixComp,...
                         model.estParams.nStates, model.estParams.nStates,...
                         'Pseudostarts', model.estParams.piPrior,...
                         'Pseudotransitions', model.estParams.transPrior);
          model.estParams.pi        = newParams.ST;
          model.estParams.trans     = newParams.TR;
        else % if isempty(mixComp)
          model.estParams.pi        = [];
          model.estParams.trans     = [];
        end
      otherwise
        error(['Invalid specification of neural ',...
               'trajectory extraction method']);
    end % switch(newModel.method)
    if (cvf > 0)
      savefield(sprintf('%s_cv%02dof%02d',...
                        newfname, cvf, oldModel.numFolds),...
                model)
    else
      savefield(newfname, model)
    end
  end % for cvf=0:oldModel.numFolds
end

