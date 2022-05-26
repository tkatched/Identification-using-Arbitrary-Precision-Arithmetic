
function omega  = createcov_lp(para,msel,RC)

% create covariance matrix for the Leeper (1991) model. 

%paranames(msel=1): theta=%[alpha beta gam sigt sigpsi mat1 mapsi1]';

%paranames(msel=2): theta=%[alpha beta gam sigt sigpsi mat1 mat2 mapsi1 mapsi2]';

%paranames(msel=3): theta=%[alpha beta gam sigt sigpsi rhot1 mat1 rhopsi1 mapsi1]';
sigt   = para(4);
sigpsi   = para(5);

%% Set parameters governing indeterminacy

if RC==[1;1] %determinacy
    m_eta_t =0;
    m_eta_p = 0;
    sigeta = 0;
elseif RC(1,1)==1 && RC(2,1)==0 && msel==1 %MA(1) shock spec has 7 parameters under determinacy
    m_eta_t = para(8);
    m_eta_p = para(9);
    sigeta = para(10);
else %MA(2) or ARMA(1,1) shock spec has 9 parameters under determinacy
    m_eta_t = para(10);
    m_eta_p = para(11);
    sigeta = para(12);
end


%% Create covariance matrix using rotations
omega = zeros(3,3);

omega(1,1) = sigt^2;
omega(2,2) = sigpsi^2;
omega(3,3)=sigeta^2;

rotmat = [1,0, m_eta_t;0, 1,m_eta_p;0, 0,1];
rotmat=rotmat';

omega = rotmat*omega*rotmat';

end