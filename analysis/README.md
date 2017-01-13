ANALYSIS OF RESULTS
===

We provide here all necessary data to analyze and reproduce the figures of
_Periodic I/O scheduling for super-computers_ by Aupy, Gainaru and Le Fèvre.

The results are for the 10 sets of experiments presented in Table 2.

DATA
---
  * 'expes_syseff' and 'expes_dil' are the results for SysEff and Dilation for
    - Experiments:
    	+ Online (syseff/dilation_online_expes)
    	+ Periodic (syseff/dilation_periodic_expes)
		+ No-additional schedule (syseff_cong_expes)
    - Simulated:
    	+ Online (syseff/dilation_online_simu)
    	+ Periodic (syseff/dilation_periodic_simu)
    - Upper Limit for system efficiency (syseff_UL)

  * study_kp:
  results obtained when varying k' from 1 to 20 + k'=100, for epsilon=0.01.

  * results_kp_10 (resp. results_kp_20):
  results obtained at each iteration of the main loop of Algorithm 2, namely for
  all values of T (T increases by a factor of 1+epsilon until Tmax), for
  epsilon=0.01, and Tmax = 10Tmin (resp. 20Tmin)


SCRIPTS
---
  * perf_expes.r is the R script used to generate Figure 4

  * period.r is the R script used to generate Figure 5

  * study_k'.r is the R script used to generate Figure 6


