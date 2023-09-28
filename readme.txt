Some simple functions for tensor method in latent class models

Note: The tensor toolbox is required to run the programs.

The file "tensor for LCM" contains the following useful functions in tensor-EM method:
eigen: Obtain all the eigenpairs of orthogonal decomposable tensor T.
eigenpairs: Obtain an eigenpair of an orthogonal decomposable tensor T.
Whiten: Whiten M3 with M2 and obtain the parameters in M2 and M3.
multi: Use tensor method to learn parameters from latent class model.
MPinv: Calculate the M-P pseudo inverse of A with only first r singular values.
EM2: EM algorithm for latent class model with specified initial values.
EM_random: EM algorithm with n random initials.
GIC: Calculate the GIC value.
GIC_fix: calculate the GIC value for fixed-effect LCM.
likelihood: Calculate the log-likelihood.
CEM: Classification EM for fixed-effect latent class model.
CEM_likeli: Calculate the log-likelihood in CEM algorithm (fixed-effect LCM).
GIC_fix: Calculate the GIC value for fixed-effect LCM.
evaluate_error: Evaluate error rates of estimated membbership.
Match: Match estimates and true parameters. This function is required since latent classes are identified up to permutations.
spectral_cluster: Spectral clustering using the normalized random-walk Laplacian.


Some simulation examples are provided in simu_example.m. Simulations under local dependence are in simu_dependence.m. The real_data.m file contains codes for real data analysis. timss.R contains codes to preprocess the TIMSS data in R (multiple imputaions). Five txt files are the complete datasets from multiple imputations. Compare_clustering.m contains codes for comparisons of different clustering algorithms.

The file "LCM_more_simu" contains similar functions that are modified so that EM-random and tensor-EM start from the same initial values. The EM-random algorithm with refined tolerance is also implemented.

