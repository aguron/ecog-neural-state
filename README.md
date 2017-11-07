# ecog-neural-state
Latent variable models for electrocorticographic (ECoG) signal analysis


==========
COPYRIGHT
==========

@ 2017 Akinyinka Omigbodun	aomigbod@ucsd.edu
       
This code may not be distributed in any form without prior consent
of the author.  While it has been carefully checked for bugs, the
code comes with no warranty of any kind.  The author takes no 
responsibility for support or maintenance of the code, although feedback
is appreciated.  

============
INSTALLATION
============

1. Download addpathrecurse.m (https://www.mathworks.com/matlabcentral/fileexchange/21086-recursive-addpath), 
   and add to the MATLAB path.
2. Run install.m
3. Clone and install the GitHub repository:
     misc-functionality (https://github.com/aguron/misc-functionality)
4. Download GPFA (https://users.ece.cmu.edu/~byronyu/software.shtml) and clone the GitHub repository
   ecog-gpfa-functionality (https://github.com/aguron/ecog-gpfa-functionality). Move the files from 
   ecog-gpfa-functionality to gpfa. Please use install.m (from ecog-gpfa-functionality) instead of 
   startup.m (in gpfa) to set up GPFA.
5. Clone the GitHub repository pmtk3 (https://github.com/probml/pmtk3) and
   initialize with initPmtk3.m.
6. Please see the *_setup files in ecog-neural-state/demos for information
   on the analyses in papers and manuscripts. Models for generating
   synthetic data are in ecog-neural-state/demos/models/models.mat