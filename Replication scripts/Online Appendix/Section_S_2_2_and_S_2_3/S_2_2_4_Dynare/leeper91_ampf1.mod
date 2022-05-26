//Leeper (1991) cashless version model

//Define variables
var pi b psi th thma  psima;

// Define shocks
varexo eth epsi;

//Define parameters
parameters alpha beta gam mat1 mapsi1;

//Parameter calibration
alpha      = 1.5; 
beta       = 0.9804; 
gam        = 1.2; 
sigt       = 1; 
sigpsi     = 1;
mat1 = 0.5;
mapsi1 = 0.5; 


//Linearized model
model(linear); 

-alpha*pi=-pi(+1)+th;
th=thma+mat1*thma(-1);
thma=eth;

b=-pi*(1/beta)+(1/beta - gam*(1/beta-1))*b(-1)+(alpha/beta)*pi(-1)+(1/beta - 1)*psi+(1/beta)*th(-1);
psi=psima+mapsi1*psima(-1);
psima=epsi;


end;

//Shock variances
shocks;
var eth = sigt^2;
var epsi = sigpsi^2;
end;

//Specify observables
varobs pi b;

//Run identification analysis 
//analytic derivative: =-1 if numerical; =0 if analytical using sylvester (default)  
// normalize_jacobians=1 use normalizations (default), =0 no normalization
identification(order=1,parameter_set=calibration,analytic_derivation_mode = -1, normalize_jacobians = 0,no_identification_strength);