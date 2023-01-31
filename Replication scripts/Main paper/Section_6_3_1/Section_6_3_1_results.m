%This script replicates the main results in Subsection 6.3.1 (global identification
%of the medium scale model, supplementary Appendix Tables 4-5)

clear
%% Global identification of the SGU (2012) model (Subsection 6.3.1)


run_sgul_nsearch % global id check with neighborhood size c=0.1, 0.5, 1

clear

run_sgul_nsearch_mp % verify the results above using quadruple precision

clear

run_sgul_nsearch_1f % global id check with sigma8_mu fixed, c=0.1, 0.5, c=1

clear

run_sgul_nsearch_1f_mp % % verify the results above using quadruple precision

clear

run_sgul_nsearch_7f % global id check with 7 weakly identified parameters fixed, c=1

clear

run_sgul_nsearch_7f_mp % % verify the results above using quadruple precision

