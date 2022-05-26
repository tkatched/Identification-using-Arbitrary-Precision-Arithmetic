List of files in the General\Double folder:

Matlab functions:

psout.m - output function for particle swarm optimization (PSO). Extracts the swarm every 200 iterations as well as the final one. The extracted swarm provides starting points for the multistart algorithm.

cols.m - returns number of columns in a matrix.

rows.m - returns number of rows in a matrix.

gensys_mod2.m - the gensys algorithm modified to allow for indeterminacy. The modification follows Lubik and Schorfheide (2003). 

legendre_rulem.m - Gaussian quadrature code by John Burkardt, distributed under the GNU LGPL license (see license.txt)

svdrr.m - Matlab translation of the Lubik and Schorfheide (2004) Gauss routine svdrr

struct2array.m - auxiliary function that may be missing on newer Matlab versions, part of Dynare

vec.m - auxiliary function that may be missing on newer Matlab versions, part of Dynare

Matlab scripts:

GA_optim.m - a script file that performs Genetic Algorithm optimization, followed by the Multistart using points from the final population as initial values. 

PSO_optim.m - a script file that performs PSO optimization, followed by the Multistart using points from the final swarm as initial values. 




