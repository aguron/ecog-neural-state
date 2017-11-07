%%
% Mutual Information and Prediction Error statistical tests
metrics                    	= {'mInf', 'predError'};

metrics_result_selection   	= [13 14 15];
markov_result_selection     = [23 24 25 26 55 56 57];
covariance_result_selection	= [31 37];

comparisons                 = [metrics_result_selection,...
                               markov_result_selection,...
                               covariance_result_selection];

if (~synthetic)
 str                        = 'real_data';
else % if (synthetic)
 str                        = 'synthetic_data';
end
for comparison=comparisons
 if exist(sprintf('results/sequence_analysis/%s/ht.mat', str), 'file')
   break
 end
 switch(comparison)
  case 1 % GMM_{diag} versus GMM_{full}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'gmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = false;
  case 2 % MFA_{CVPN} versus GMM_{full}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 3 % MFA_{CIPN} versus GMM_{full}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 4 % MFA_{CVPN} versus GMM_{diag}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 5 % MFA_{CIPN} versus GMM_{diag}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 6 % MFA_{CIPN} versus MFA_{CVPN}
   str1                   	= '';
   method1                	= 'mfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 7 % GMM_{diag}-MM versus GMM_{full}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '/gmm-mm';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = false;
  case 8 % MFA_{CVPN}-MM versus GMM_{full}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 9 % MFA_{CIPN}-MM versus GMM_{full}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 10 % MFA_{CVPN}-MM versus GMM_{diag}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 11 % MFA_{CIPN}-MM versus GMM_{diag}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 12 % MFA_{CIPN}-MM versus MFA_{CVPN}-MM
   str1                   	= '/mfa-mm';
   method1                	= 'hmfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
   
   
   
  case 13 % HMM_{diag} versus HMM_{full}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = false;
  case 14 % HMFA_{CVPN} versus HMM_{full}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 15 % HMFA_{CIPN} versus HMM_{full}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 16 % HMFA_{CVPN} versus HMM_{diag}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 17 % HMFA_{CIPN} versus HMM_{diag}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 18 % HMFA_{CIPN} versus HMFA_{CVPN}
   str1                   	= '';
   method1                	= 'hmfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
   
   
   
  case 19 % GMM_{full}-MM versus GMM_{full}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '/gmm-mm';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = false;
  case 20 % GMM_{diag}-MM versus GMM_{diag}
   str1                   	= '';
   method1                	= 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '/gmm-mm';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = false;
  case 21 % MFA_{CVPN}-MM versus MFA_{CVPN}
   str1                   	= '';
   method1                	= 'mfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 22 % MFA_{CIPN}-MM versus MFA_{CIPN}
   str1                   	= '';
   method1                	= 'mfa';
   faType1                	= [1 1 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;

   
   
   
  case 23 % HMM_{full} versus GMM_{full}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = false;
  case 24 % HMM_{diag} versus GMM_{diag}
   str1                   	= '';
   method1                	= 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = false;
  case 25 % HMFA_{CVPN} versus MFA_{CVPN}
   str1                   	= '';
   method1                	= 'mfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 26 % HMFA_{CIPN} versus MFA_{CIPN}
   str1                   	= '';
   method1                	= 'mfa';
   faType1                	= [1 1 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;
   
   
   
   
  case 27 % HMM_{full} versus GMM_{full}-MM
   str1                   	= '/gmm-mm';
   method1                	= 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = false;
  case 28 % HMM_{diag} versus GMM_{diag}-MM
   str1                   	= '/hmm0';
   method1                	= 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = false;
  case 29 % HMFA_{CVPN} versus MFA_{CVPN}-MM
   str1                   	= '/mfa-mm';
   method1                	= 'hmfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 1];
   covType2               	= '';
   sharedCov2               = nan;
  case 30 % HMFA_{CIPN} versus MFA_{CIPN}-MM
   str1                   	= '/mfa-mm';
   method1                	= 'hmfa';
   faType1                	= [1 1 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 1 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 31 % GMM_{diag, tied} versus GMM_{full, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                   	= '';
   method2                	= 'gmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 32 % MFA_{tied} versus GMM_{full, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 33 % MFA_{tied} versus GMM_{diag, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = true;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 34 % GMM_{diag, tied}-MM versus GMM_{full, tied}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                     = '/gmm-mm';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 35 % MFA_{tied}-MM versus GMM_{full, tied}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 36 % MFA_{tied}-MM versus GMM_{diag, tied}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = true;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 37 % HMM_{diag, tied} versus HMM_{full, tied}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 38 % HMFA_{tied} versus HMM_{full, tied}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 39 % HMFA_{tied} versus HMM_{diag, tied}
   str1                     = '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = true;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 40 % GMM_{full, tied} versus GMM_{full}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'gmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = true;
  case 41 % GMM_{diag, tied} versus GMM_{diag}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                	= 'gmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 42 % MFA_{tied} versus MFA_{CVPN}
   str1                     = '';
   method1                  = 'mfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 43 % MFA_{tied} versus MFA_{CIPN}
   str1                     = '';
   method1                  = 'mfa';
   faType1                	= [1 1 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'mfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 44 % GMM_{full, tied}-MM versus GMM_{full}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                     = '/gmm-mm';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = true;
  case 45 % GMM_{diag, tied}-MM versus GMM_{diag}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                     = '/gmm-mm';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 46 % MFA_{tied}-MM versus MFA_{CVPN}-MM
   str1                   	= '/mfa-mm';
   method1                	= 'hmfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 47 % MFA_{tied}-MM versus MFA_{CIPN}-MM
   str1                   	= '/mfa-mm';
   method1                	= 'hmfa';
   faType1                	= [1 1 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 48 % HMM_{full, tied} versus HMM_{full}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = false;

   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = true;
  case 49 % HMM_{diag, tied} versus HMM_{diag}
   str1                   	= '';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = false;

   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 50 % HMFA_{tied} versus HMFA_{CVPN}
   str1                   	= '';
   method1                	= 'hmfa';
   faType1                	= [1 1 1];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  case 51 % HMFA_{tied} versus HMFA_{CIPN}
   str1                   	= '';
   method1                	= 'hmfa';
   faType1                	= [1 1 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 52 % GMM_{full, tied}-MM versus GMM_{full, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                     = '/gmm-mm';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = true;
  case 53 % GMM_{diag, tied}-MM versus GMM_{diag, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = true;

   str2                     = '/gmm-mm';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 54 % MFA_{tied}-MM versus MFA_{tied}
   str1                     = '';
   method1                  = 'mfa';
   faType1                	= [1 0 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '/mfa-mm';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 55 % HMM_{full, tied} versus GMM_{full, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;
   
   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = true;
  case 56 % HMM_{diag, tied} versus GMM_{diag, tied}
   str1                     = '';
   method1                  = 'gmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = true;
   
   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 57 % HMFA_{tied} versus MFA_{tied}
   str1                     = '';
   method1                  = 'mfa';
   faType1                	= [1 0 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;




  case 58 % HMM_{full, tied} versus GMM_{full, tied}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'full';
   sharedCov1               = true;

   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'full';
   sharedCov2               = true;
  case 59 % HMM_{diag, tied} versus GMM_{diag, tied}-MM
   str1                     = '/gmm-mm';
   method1                  = 'hmm';
   faType1                	= [nan nan nan];
   covType1               	= 'diagonal';
   sharedCov1               = true;
   
   str2                   	= '';
   method2                  = 'hmm';
   faType2                	= [nan nan nan];
   covType2               	= 'diagonal';
   sharedCov2               = true;
  case 60 % HMFA_{tied} versus MFA_{tied}-MM
   str1                   	= '/mfa-mm';
   method1                	= 'hmfa';
   faType1                	= [1 0 0];
   covType1               	= '';
   sharedCov1               = nan;

   str2                   	= '';
   method2                	= 'hmfa';
   faType2                	= [1 0 0];
   covType2               	= '';
   sharedCov2               = nan;
  otherwise
   error('Invalid comparison index');
 end % switch(comparison)
 
 for runIdx=pa_select
  for iB=binWidth_select
   for iF=1:numel(fracTrainData)
    for nStates=nStatesRange
     for xDim=xDimRange
      for metric=metrics
       fprintf('comparison: %d\n', comparison);
       results(1).dir      	= sprintf('%s%s', resultsDir, str1);
       results(1).runIdx  	= runIdx;
       results(1).binWidth 	= binWidth(iB);
       results(1).nFolds   	= 4;
       results(1).method   	= method1;
       if ismember(method1, {'mfa', 'hmfa'})
        results(1).xDim   	= xDim;
       else % if ismember(method1, {'gmm', 'hmm'})
        results(1).xDim    	= nan;
       end
       if ismember(method1, {'gmm', 'mfa'})
        results(1).nMixComp	= nStates;
        results(1).nStates	= nan;
       else % if ismember(method1, {'hmm', 'hmfa'})
        results(1).nMixComp	= nan;
        results(1).nStates	= nStates;
       end
       results(1).faType   	= faType1;
       results(1).kernSD   	= 0;
       results(1).covType  	= covType1;
       results(1).sharedCov	= sharedCov1;

       
       results(2).dir      	= sprintf('%s%s', resultsDir, str2);
       results(2).runIdx  	= runIdx;
       results(2).binWidth 	= binWidth(iB);
       results(2).nFolds   	= 4;
       results(2).method   	= method2;
       if ismember(method2, {'mfa', 'hmfa'})
        results(2).xDim    	= xDim;
       else % if ismember(method1, {'gmm', 'hmm'})
        results(2).xDim    	= nan;
       end
       if ismember(method1, {'gmm', 'mfa'})
        results(2).nMixComp	= nStates;
        results(2).nStates	= nan;
       else % if ismember(method1, {'hmm', 'hmfa'})
        results(2).nMixComp	= nan;
        results(2).nStates	= nStates;
       end
       results(2).faType   	= faType2;
       results(2).kernSD   	= 0;
       results(2).covType  	= covType2;
       results(2).sharedCov	= sharedCov2;

       try
         ht               	=...
           comparemetric(results,...
                         'metric',metric{1},...
                         'stat', struct('method','exact',...
                                        'alpha',[],...
                                        'maxNumMeasFold',[]));
         if any(ismember({method1, method2}, {'mfa', 'hmfa'}))
           hyptest(runIdx,iB,iF,nStates,xDim,comparison).(metric{1})...
                            = ht;
         else
           hyptest(runIdx,iB,iF,nStates,1,comparison).(metric{1})...
                            = ht;
         end
       catch err
        displayerror(err)
        programcontrol
       end % try
      end % for metric=metrics
     end % for xDim=xDimRange
    end % for nStates=nStatesRange
   end % for iF=1:numel(fracTrainData)
  end % for iB=binWidth_select
 end % for runIdx=pa_select
end % for comparison=comparisons
if ~exist(sprintf('results/sequence_analysis/%s/ht.mat', str), 'file')
  save(sprintf('results/sequence_analysis/%s/ht.mat', str), 'hyptest');
end
%%
for comparison_set={metrics_result_selection,...
                    markov_result_selection,...
                    covariance_result_selection}
 for comparison=comparison_set{1}
  for metric=metrics
   if isequal(comparison_set{1}, metrics_result_selection)
    selectF_selection       = 1:3;
   else
    selectF_selection       = 3;
   end
   for selectF=selectF_selection % fraction of training data
    load(sprintf('results/sequence_analysis/%s/ht.mat', str), 'hyptest');
    selectB                	= 1:3; % bin width
    selectS               	= 2:4; % Number of states or mixture components
    switch(comparison)
     case {1, 7, 13, 19,20 23,24, 27,28,...
           31, 34, 37,...
           40,41, 44,45, 48,49,...
           52,53, 55,56, 58,59}
      selectL              	= [];
      beforeBHP            	= {pa_select,selectB,[],selectS,1,comparison};
      afterBHP             	= {[],[],selectF,[],[],[]};
     case {2,3,4,5, 8,9,10,11, 14,15,16,17,...
           32,33, 35,36, 38,39}
      selectL              	= 1;
      beforeBHP            	=...
       {pa_select,selectB,[],selectS,selectL,comparison};
      afterBHP             	= {[],[],selectF,[],[],[]};
     case {6, 12, 18, 21,22, 25,26, 29,30,...
           42,43, 46,47, 50,51,...
           54, 57, 60}
      selectL              	= 1:size(hyptest,5);
      beforeBHP            	=...
       {pa_select,selectB,[],selectS,selectL,comparison};
      afterBHP             	= {[],[],selectF,[],[],[]};
     otherwise
      error('Invalid comparison index');
    end % switch(comparison)

    clear('str2');
    str2{1}               	=...
     strtok(hyptest(1,1,1,1,1,comparison).(metric{1}).model(2).name,'_');
    [str2{3}, str2{2}]     	=...
     rstrtok(hyptest(1,1,1,1,1,comparison).(metric{1}).model(2).name,'_');
    if isequal(str2{2},'_tied')
     [~, str2{3}]          	= rstrtok(str2{3},'_');
    else % if ~isequal(str2{2},'_tied')
     str2{3}               	= '';
    end
    str2{4}                	=...
     strtok(hyptest(1,1,1,1,1,comparison).(metric{1}).model(1).name,'_');
    [str2{6}, str2{5}]     	=...
     rstrtok(hyptest(1,1,1,1,1,comparison).(metric{1}).model(1).name,'_');
    if isequal(str2{5},'_tied')
     [~, str2{6}]          	= rstrtok(str2{6},'_');
    else % if ~isequal(str2{5},'_tied')
     str2{6}                = '';
    end

    modelclasscomparisons(hyptest, beforeBHP, afterBHP,...
                          {{'Comparison:', '%d', comparison},...
                           {'Patient:', '%d', pa_select},...
                           {'Bin width:', '%g', binWidth(selectB)},...
                           {'Number of states:', '%d', selectS},...
                           {'Latent dimensionality:', '%d', selectL},...
                           {'Fraction of training data:', '%g', fracTrainData(selectF)},...
                           },...
                          'metric', metric{1},...
                          'FDR', 0.05,...
                          'modelClassNames', {[str2{1:3}],[str2{4:6}]})
   end % for selectF=selectF_selection
  end % for metric=metrics
 end % for comparison=comparison_set{1}
end