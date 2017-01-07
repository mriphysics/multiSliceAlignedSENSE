# multiSliceAlignedSENSE
Tools for aligned reconstruction of multi-shot multi-slice MR.

This repository provides tools to implement the reconstruction methods and reproduce the experiments included in Figures 3 and 4 of the manuscript ''3D motion corrected SENSE reconstruction for multi-shot multi-slice MRI: application to neonatal brain imaging'', L. Cordero-Grande, E. J. Hughes, J. Hutter, A. N. Price, and J. V. Hajnal. Unpublished.

The code has been developed in MATLAB and has the following structure:

###### ./
contains the scripts for running the experiments included in Figures 4 and 3 of the manuscript, respectively: *multiSliceAlignedSENSE_Exp[1-2].m*.

###### ./data
contains the datasets used for simulations, respectively for Figures 4 and 3: *yT[1-2].mat*. Data generated when running the scripts is also stored in this folder.

###### ./lib
contains a package used to save the resulting reconstructions as nifti files: *NIfTI_20140122*.

###### ./meth
contains the solvers for motion and reconstruction as well as the alternating method described in the manuscript: *solve[X,T]MS2D.m*, *optimizeLevelMS2D.m*.

###### ./meth/pre
contains functions used for preprocessing: *coilArrayCompression.m*.

###### ./meth/pro
contains functions used by the solvers: *fftGPU.m*, *ifftGPU.m*, *sense.m*, *isense.m*, *precomputeFactors3DTransform.m*, *transform3DSinc.m*, *transform3DSinc[Gradient,Hessian].m*, *filtering.m*, *mirroring.m*, *resampling.m*, *computeROI.m*, *extractROI.m*, *extractSlabs.m*, *sliceProfile.m*.

###### ./meth/pos
contains functions for postprocessing: *gibbsRingingFilter.m*, *rotateMPS.m*, *writeNIFTI.m*.

NOTE 1: When running the scripts, set the *gpu* variable to 0 if your system does not support NVIDIA GPU computing under MATLAB.

NOTE 2: Computation times on an 8(16) x Intel(R) Core(TM) i7-5960X CPU @ 3.00GHz 64GB RAM with a GeForce GTX TITAN X have been:
* Figure 4 / Exp 1:
  * Conventional uncorrected SENSE reconstruction: 4".
  * Uncorrected reconstruction with slice profile filter: 17".
  * Corrected reconstruction without outlier rejection: 17'33".
  * Corrected reconstruction without through-plane motion: 12'55".
  * Fully corrected reconstruction: 15'14".
* Figure 3 / Exp 2:
  * Conventional uncorrected SENSE reconstruction: 3".
  * Uncorrected reconstruction with slice profile filter: 21". 
  * Corrected reconstruction without outlier rejection: 
  * Corrected reconstruction without through-plane motion:
  * Fully corrected reconstruction: 

