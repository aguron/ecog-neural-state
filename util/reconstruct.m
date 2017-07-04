function [seq, cvfolds, obj] = reconstruct(seq, estParams, method, varargin)
%RECONSTRUCT the neural signal according to specifications for the
%   latent variable model

% OPTIONAL ARGUMENTS:
%
% numFolds	- number of folds for K-means and GMMs (default: 0)
% nRuns     - number of K-means/GMM runs
% verbose   - logical that specifies whether to display status messages
%             (default: false)

  numFolds                    	= 0;
  nRuns                         = 100;
  verbose                       = false;

  assignopts(who, varargin);
  
  cvfolds                       = [];
  obj                           = [];
  
  nSeq                          = numel(seq);
  switch(method)
    case 'pca'
      for s=1:nSeq
        seq(s).(method)         =...
          bsxfun(@plus, estParams.Corth * seq(s).xorth, estParams.d);
      end % for s=1:nSeq
    case {'ppca', 'fa'}
      for s=1:nSeq
        seq(s).(method)         =...
          bsxfun(@plus, estParams.C * seq(s).xpost, estParams.d);
      end % for s=1:nSeq
    case 'gpfa'
      for s=1:nSeq
        seq(s).(method)         =...
          bsxfun(@plus, estParams.C * seq(s).xsm, estParams.d);
      end % for s=1:nSeq
    case 'hmfa'
      for s=1:nSeq
        % E[x|s] = \mu_s  
        seq(s).mean             = estParams.d(:,seq(s).state);
        
        T                       = seq(s).T;
        seq(s).factors          = nan(size(seq(s).y));
        for t=1:T
          % E[x|s,z] - E[x|s] = \lambda_s * z
          seq(s).factors(:,t)   =...
            estParams.C(:,:,seq(s).state(t)) * seq(s).x(:,t,seq(s).state(t));
        end % for t=1:T
        % E[x|s,z] = \mu_s + \lambda_s * z
        seq(s).mean_factors     = seq(s).mean + seq(s).factors;
      end % for s=1:nSeq
    case 'k-means'
      nClusters                 = estParams;

      if (numFolds)
        fdiv                 	= floor(linspace(1, nSeq+1, numFolds+1));
        cvfolds                 =...
          struct('centroid', cell(1,numFolds),...
                 'cluster', cell(1,numFolds),...
                 'seqTrain', cell(1,numFolds),...
                 'seqTest', cell(1,numFolds),...
                 'sumdTest', cell(1,numFolds));
        % Set cross-validation folds
        for cvf=1:numFolds
          fprintf('\n===== Cross-validation fold %d of %d =====\n',...
                  cvf, numFolds);
          fprintf('Number of K-means clusters: %d\n', nClusters);
          % Set cross-validation masks
          testMask           	= false(1, nSeq);
          testMask(fdiv(cvf):fdiv(cvf+1)-1)...
                               	= true;
          trainMask           	= ~testMask;
          
          % Randomly reorder trials before partitioning
          % into training and test sets
          rng('default');
          tr                  	= randperm(nSeq);
          trainTrialIdx        	= tr(trainMask);
          testTrialIdx         	= tr(testMask);
          seqTrain             	= seq(trainTrialIdx);
          seqTest              	= seq(testTrialIdx);
          
          yTrain               	= [seqTrain.y];
          [cluster, centroid]  	=...
            kmeanscluster(yTrain', nClusters,...
                          'nRuns', nRuns,...
                          'verbose', verbose);

          yTest                	= [seqTest.y];
          [~, ~, sumdTest]     	=...
            kmeanscluster(yTest', nClusters,...
                          'centroid', centroid);

          cvfolds(cvf).centroid	= centroid;
          cvfolds(cvf).cluster 	= cluster;
          cvfolds(cvf).sumdTest	= sumdTest;

          cvfolds(cvf).seqTrain	=...
            rmfield(seqTrain,...
                    setdiff(fieldnames(seqTrain),{'trialId', 'y', 'T'}));

          cvfolds(cvf).seqTest 	=...
            rmfield(seqTest,...
                    setdiff(fieldnames(seqTest),{'trialId', 'y', 'T'}));
        end % for cvf=1:numFolds
      end % if (numFolds)
      
      fprintf('\n===== Training on all data =====\n');
      fprintf('Number of K-means clusters: %d\n', nClusters);
      
      yAll                      = [seq.y];
      [cluster, centroid]       = kmeanscluster(yAll', nClusters,...
                                                'nRuns', nRuns,...
                                                'verbose', verbose);
      centroid                  = centroid';
      seq                       = segmentByTrial(seq, cluster', 'cluster');
      for s=1:nSeq
        seq(s).kmeans           = centroid(:, seq(s).cluster);
      end % for s=1:nSeq
    case 'gmm'
      nMixComp                  = estParams;

      if (numFolds)
        fdiv                 	= floor(linspace(1, nSeq+1, numFolds+1));
        cvfolds                 =...
          struct('mixComp', cell(1,numFolds),...
                 'obj', cell(1,numFolds),...
                 'seqTrain', cell(1,numFolds),...
                 'seqTest', cell(1,numFolds),...
                 'LLtest', cell(1,numFolds));
        % Set cross-validation folds
        for cvf=1:numFolds
          fprintf('\n===== Cross-validation fold %d of %d =====\n',...
                  cvf, numFolds);
              
          fprintf('Number of GMM mixture components: %d\n', nMixComp);
          % Set cross-validation masks
          testMask           	= false(1, nSeq);
          testMask(fdiv(cvf):fdiv(cvf+1)-1)...
                               	= true;
          trainMask           	= ~testMask;
          
          % Randomly reorder trials before partitioning
          % into training and test sets
          rng('default');
          tr                  	= randperm(nSeq);
          trainTrialIdx        	= tr(trainMask);
          testTrialIdx        	= tr(testMask);
          seqTrain             	= seq(trainTrialIdx);
          seqTest             	= seq(testTrialIdx);
          
          yTrain             	= [seqTrain.y];
          [mixComp, obj]       	=...
            gmcluster(yTrain', nMixComp,...
                      'nRuns', nRuns,...
                      'verbose', verbose);

          yTest                	= [seqTest.y];
          [~, LLtest]         	= posterior(obj,yTest');
          LLtest               	= -LLtest;

          cvfolds(cvf).obj     	= obj;
          cvfolds(cvf).mixComp  = mixComp;
          cvfolds(cvf).LLtest   = LLtest;

          cvfolds(cvf).seqTrain         =...
            rmfield(seqTrain,...
                    setdiff(fieldnames(seqTrain),{'trialId', 'y', 'T'}));

          cvfolds(cvf).seqTest          =...
            rmfield(seqTest,...
                    setdiff(fieldnames(seqTest),{'trialId', 'y', 'T'}));
        end % for cvf=1:numFolds
      end % if (numFolds)
      
      fprintf('\n===== Training on all data =====\n');
      fprintf('Number of GMM mixture components: %d\n', nMixComp);
      
      yAll                      = [seq.y];
      [mixComp, obj]            = gmcluster(yAll', nMixComp,...
                                            'nRuns', nRuns,...
                                            'verbose', verbose);
      mu                        = obj.mu';
      seq                       = segmentByTrial(seq, mixComp', 'mixComp');
      for s=1:nSeq
        seq(s).gmmeans          = mu(:,seq(s).mixComp);
      end % for s=1:nSeq
    otherwise
      error('Invalid method specification');
  end % switch(method)
end