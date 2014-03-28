Motion_Correction
=================

Code for motion correcting calcium imaging movies.

Uses spatial subsegments of full movie to compute local x-y shift within individual acquisitions, 
computes global shift across acquisitions, then stitches together acquisition and framewise segment shifts
by fitting an affine geometric transform.

Only files in 'Affine Correction' folder are needed to run code. 'RunAffine' script calls all functions in order.
'InitMovieFiles' specifies filenames, locations, etc.

Uses smart memory management to load only individual segments into memory at a time, including for
calculation of the covariance matrix / PCA processing for ICA_PCA algorithm input. Uses parallel computing toolbox 
to speed up computations. Note: should probably not run multiple instances of this code simultaneously.
