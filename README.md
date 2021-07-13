# PLS-CCA
Matlab code and data used for comparison of partial least squares and canonical correlation
MAT files are matlab workspaces containing the X and Y data blocks for each data example

pls_cancor is the workhorse that does both PLS and CCA on two matrices including permutation and bootstrap resampling

split_half_PLSCCA performs split half tests where singular vectors from PLS & CCA on random splits of the data are calculated and compared (cosine) and the distribution of cosines saved. Z-values are created to ascertain the reproducibility of latent variables

split_half_PLSCCA_TrainTest also performs split half but uses the singular vectors from one half to compute the singular values from the other. The "test" singular values are assessed for the strength effect, similar to an R^2 value computed on test sample in regression analysis
