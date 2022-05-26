The list of files in the folder \Replication scripts\Online Appendix\Section_S_3_8

Matlab scripts:

Section_S_3_8_results.m - a master script for replicating results in the Section S.3.8 of the Online Appendix. The code structure is essentially the same as in the main paper Sections 6.2.1-6.2.3, the only difference here being that the analysis is using the ML estimate from Table II in SGU (2012) as the benchmark model.



Subfolders:

Global identification - contains the scripts to check global identification of parameters in the baseline SGU (2012) model in double and arbitrary precision

Identification of news shocks- contains the scripts to check global identification of the news shocks in the baseline SGU (2012) model (shutting down news altogether, by horizon, and 1-by-1 )

Equivalence with ARMA shocks - contains the scripts to check global identification across models structures: news shocks vs. MA shocks (replacing all shocks by MA couterparts, as well as 1-by-1)