%%
% The analyses below are detailed in the following manuscript to be
% submitted for publication:
%
% Omigbodun, A., Doyle, W.K., Devinsky, O., Friedman, D., Melloni, L., 
% Thesen, T., and Gilja, V. (2017) Discrete state based approaches to 
% analyzing electrocorticographic data.
%%
disp('****************************************');
disp('* Start neural state sequence analysis *');
disp('****************************************');
prevpwd             	= pwd;
%%
% Specify whether the actual dataset or synthetic data will be used
% in the analyses
synthetic             = false;
analyses              = true;
%
if (~synthetic)
 %%
 % Real data preparation
 dataDir              = []; % specify data file location
 load([dataDir, '/dat.mat'],'dat');
 
 % Algorithm configuration
 pa_select            = 1:6; % Number of patients

 binWidth             = [0.05 0.1 0.2];
 binWidth_select      = 1:numel(binWidth);

 tolGMM               = 1e-6;
 tolHMM               = 1e-2;
 tolMFA               = 1e-2;
 tolHMFA              = 1e-2;

 methodList           = {'gmm', 'hmm', 'mfa', 'hmfa'};
 xDimRange            = 1:10;
 nStatesRange         = 1:4; % Number of states or mixture components
 faTypeRange          = {[1 1 1], [1 1 0], [1 0 0]};
 covTypeRange         = {'full', 'diagonal'};
 sharedCovRange       = {true, false};

 prediction           = true;
 
 fracTrainData        = [0.1 0.5 1]; % Fraction of training data
 for frac=fracTrainData
   resultsDir        	=...
    [pwd, sprintf('/results/sequence_analysis/real_data/%g', frac)];
   if (analyses)
     sequence_analysis
   end % if (analyses)
 end % for frac=fracTrainData
else % if (synthetic)
 %%
 % Synthetic data preparation
 load('models/models', 'model', 'trialType')
 nSeq                 = numel(trialType);
 seq                  =...
  struct('trialId', cell(1,nSeq),...
         'fs', cell(1,nSeq),...
         'T', cell(1,nSeq),...
         'y', cell(1,nSeq),...
         'state', cell(1,nSeq));
 genModel             = 'hmfa';
 faType               = [1 1 0];
 faTypeSpec           = 'tu'; % t - tied; u - untied
 switch(genModel)
  case {'gmm', 'hmm'}

  case {'mfa', 'hmfa'}
   seq                = addfield(seq, 'x', cell(1,nSeq));
  otherwise
   error('Invalid generating model');
 end % switch(genModel)
 
 trialLength          = 102;
 for i=1:numunique(trialType)
  seq(i==trialType)   =...
   modelSample(struct('type', genModel,...
                      'params', model(i).(genModel).(faTypeSpec(faType+1)),...
                      'fs', 512),...
             	 trialLength, sum(i==trialType));
 end % for i=1:numunique(trialType)
 
 % Algorithm configuration
 pa_select            = 1;
 binWidth             = 0.05;
 binWidth_select      = 1:numel(binWidth);

 tolGMM               = 1e-6;
 tolHMM               = 1e-2;
 tolMFA               = 1e-2;
 tolHMFA              = 1e-2;
 
 methodList           = {'gmm', 'hmm', 'mfa', 'hmfa'};
 xDimRange            = 1:10;
 nStatesRange         = 1:4; % Number of states or mixture components
 faTypeRange          = {[1 1 1], [1 1 0], [1 0 0]};
 covTypeRange         = {'full', 'diagonal'};
 sharedCovRange       = {true, false};
 
 prediction           = true;
 
 fracTrainData        = [0.1 0.5 1]; % Fraction of training data
 for frac=fracTrainData
   resultsDir        	=...
    [pwd, sprintf('/results/sequence_analysis/synthetic_data/%g', frac)];
   if (analyses)
     sequence_analysis
   end % if (analyses)
 end % for frac=fracTrainData
end
%%
figs                  = true;
if (figs)
  sequence_analysis_visualizations
end
%%
stats                 = true;
if (stats)
  sequence_analysis_stats
end
%%
cd(prevpwd)
disp('*******************************************');
disp('* Neural state sequence analysis complete *');
disp('*******************************************');
%%