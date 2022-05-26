%% ** initialize matrices **
GAM0 = zeros(neq,neq);
GAM1 = zeros(neq,neq);
C = zeros(neq,1);
PSI = zeros(neq,neps);
PPI = zeros(neq,neta);

%%      1. Capital Accumulation      ok                    

    GAM0(eq_inv, v_k1) = -1;
    GAM0(eq_inv, v_u)  = -del1/mukss;
    GAM0(eq_inv, v_muk)= -(1-del0)/mukss;
    GAM0(eq_inv, v_i)  =  i_k;
    GAM0(eq_inv, v_zi) =  i_k;
    GAM1(eq_inv, v_k1) = -(1-del0)/mukss;
 %%  2. Resource Constraint                           

    GAM0(eq_res, v_c)   = c_y;
    GAM0(eq_res, v_i)   = i_y;
    GAM0(eq_res, v_g)   = g_y*xgss;
    GAM0(eq_res, v_xg)  = g_y*xgss;
    GAM0(eq_res, v_y)   = -1;
    
    %% 3. Consumption Bundle (V_t)                         

    GAM0(eq_v, v_v)   = -(1-badj)*(1-ppsi*hss^(thet)*(1/muyss)^((1-gam)/gam));
    GAM0(eq_v, v_c)   =  1;
    GAM0(eq_v, v_muy) =  badj;
    GAM1(eq_v, v_c)   =  badj;
    GAM0(eq_v, v_h)   = -ppsi*hss^thet*(1-badj)*(1/muyss)^((1-gam)/gam)*thet;
    GAM0(eq_v, v_s)   = -ppsi*hss^thet*(1-badj)*(1/muyss)^((1-gam)/gam);


    %% 4. Consumption Decision                          

    GAM0(eq_c, v_lam)  = -(1-betbar)*(1-eta);
    GAM0(eq_c, v_zet)  =  1;
    GAM0(eq_c, v_v)    = -sig;
    GAM0(eq_c, v_c)    =  eta*(1+betbar*badj)/(1-badj);
    GAM1(eq_c, v_c)    =  eta*badj/(1-badj);
    GAM0(eq_c, v_muy)  =  eta*badj/(1-badj);
    GAM0(eq_c, v_p)    = -eta;
    GAM0(eq_c, v_s)    = -eta;    
    GAM0(eq_c, v_Ezet) = -betbar;
    GAM0(eq_c, v_Ev)   =  betbar*sig;
    GAM0(eq_c, v_Ep)   =  betbar*eta;
    GAM0(eq_c, v_Es)   =  betbar*eta;
    GAM0(eq_c, v_Ec)   = -betbar*eta/(1-badj);
    GAM0(eq_c, v_Emuy) =  betbar*((1-eta)*sig-eta*badj/(1-badj));


    %% 5. Jaimovich-Rebelo Part                        

    GAM0(eq_s,v_s)    = -1;
    GAM0(eq_s,v_c)    =  gam/(1-badj);
    GAM1(eq_s,v_c)    =  gam*badj/(1-badj);
    GAM0(eq_s,v_muy)  =  gam*badj/(1-badj)-(1-gam);
    GAM1(eq_s,v_s)    = -(1-gam);


%%  6. Hours Decision                                
    GAM0(eq_h,v_zet)  =  1;
    GAM0(eq_h,v_v)    = -sig;
    GAM0(eq_h,v_h)    =  thet;
    GAM0(eq_h,v_s)    =  1;
    GAM0(eq_h,v_lam)  = -1;
    GAM0(eq_h,v_y)    = -1;
    GAM0(eq_h,v_mu)   =  1; % The shock is on the gross wage mark-up

%% 7. Shadow Price                                 
    GAM0(eq_p,v_p)    = -1;
    GAM0(eq_p,v_zet)  =  xi;
    GAM0(eq_p,v_v)    = -xi*sig;
    GAM0(eq_p,v_h)    =  xi*thet;
    GAM0(eq_p,v_Emuy) =  (1-xi)*(1-sig);
    GAM0(eq_p,v_Ep)   =  1-xi;
    GAM0(eq_p,v_Es)   =  1-xi;
    GAM0(eq_p,v_s)    = -(1-xi);

%% 8. Output                                        
    GAM1(eq_y, v_k1) = -alpk; 
    GAM0(eq_y, v_y)  = -1;
    GAM0(eq_y, v_z)  =  1;
    GAM0(eq_y, v_u)  =  alpk;
    GAM0(eq_y, v_muk)= -alpk;
    GAM0(eq_y, v_h)  =  alph;

%% 9. Capacity Utilization                          
    GAM0(eq_cap, v_q)   = -1;
    GAM0(eq_cap, v_z)   =  1;
    GAM0(eq_cap, v_u)   =  alpk-1-del2/del1;
    GAM1(eq_cap, v_k1)  = -(alpk-1);
    GAM0(eq_cap, v_muk) = -(alpk-1);
    GAM0(eq_cap, v_h)   =  alph;

%% 10. Euler Equation                               
    GAM0(eq_eul, v_lam)  = -1;
    GAM0(eq_eul, v_q)    = -1;
    GAM0(eq_eul, v_Elam) =  1;
    GAM0(eq_eul, v_Emua) =  1;
    GAM0(eq_eul, v_Emuy) = -sig;
    GAM0(eq_eul, v_Ez)   =  bet*muass*muyss^(-sig)*del1;
    GAM0(eq_eul, v_Eu)   =  bet*muass*muyss^(-sig)*del1*(alpk-1);
    GAM0(eq_eul, v_k1)   =  bet*muass*muyss^(-sig)*del1*(alpk-1);
    GAM0(eq_eul, v_Emuk) = -bet*muass*muyss^(-sig)*del1*(alpk-1);
    GAM0(eq_eul, v_Eh)   =  bet*muass*muyss^(-sig)*del1*alph;
    GAM0(eq_eul, v_Eq)   =  bet*muass*muyss^(-sig)*(1-del0);


%% 11. Value of Firm                             
    GAM0(eq_pinv, v_q)    =  1;
    GAM0(eq_pinv, v_zi)   =  1;
    GAM0(eq_pinv, v_i)    = -kap*(bet*muass*muyss^(-sig)*mukss^3 + mukss^2);
    GAM1(eq_pinv, v_i)    = -kap*mukss^2;
    GAM0(eq_pinv, v_muk)  = -kap*mukss^2;
    GAM0(eq_pinv, v_Ei)   =  kap*bet*muass*muyss^(-sig)*mukss^3;
    GAM0(eq_pinv, v_Emuk) =  kap*bet*muass*muyss^(-sig)*mukss^3;

%% 12-14. Stochastic Trends                         
    GAM0(eq_muy, v_muy) = -1;
    GAM0(eq_muy, v_mua) =  alpk/(alpk-1);
    GAM0(eq_muy, v_mux) =  1;

    GAM0(eq_muk, v_muk) = -1;
    GAM0(eq_muk, v_mua) =  1/(alpk-1);
    GAM0(eq_muk, v_mux) =  1;

    GAM0(eq_xg,  v_xg ) =  1;
    GAM1(eq_xg,  v_xg ) =  rhoxg;
    GAM0(eq_xg,  v_muy) =  1;

%% 15+. Shock Processes                             
 if msel==1 %news   
    % MUA%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    GAM0(eq_mua,v_mua)    = 1;
    GAM1(eq_mua,v_mua)    = rhoa;
    GAM1(eq_mua,v_mua8_1) = 1;
    PSI(eq_mua,e_mua0) = 1;

    GAM0(eq_mua8_1,v_mua8_1) = 1;
    GAM1(eq_mua8_1,v_mua8_2) = 1;

    GAM0(eq_mua8_2,v_mua8_2) = 1;
    GAM1(eq_mua8_2,v_mua8_3) = 1;

    GAM0(eq_mua8_3,v_mua8_3) = 1;
    GAM1(eq_mua8_3,v_mua8_4) = 1;

    GAM0(eq_mua8_4,v_mua8_4) = 1;
    GAM1(eq_mua8_4,v_mua8_5) = 1;
    PSI(eq_mua8_4,e_mua4)    = 1;

    GAM0(eq_mua8_5,v_mua8_5) = 1;
    GAM1(eq_mua8_5,v_mua8_6) = 1;

    GAM0(eq_mua8_6,v_mua8_6) = 1;
    GAM1(eq_mua8_6,v_mua8_7) = 1;

    GAM0(eq_mua8_7,v_mua8_7) = 1;
    GAM1(eq_mua8_7,v_mua8_8) = 1;

    GAM0(eq_mua8_8,v_mua8_8) = 1;
    PSI(eq_mua8_8,e_mua8)    = 1;

    % MUX%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    GAM0(eq_mux,v_mux)    = 1;
    %GAM1(eq_mux,v_mux)    = rhox-0.5; %in the original code, prior is put
    %on rhox+0.5
    GAM1(eq_mux,v_mux)    = rhox; %here just use rhox directly
    GAM1(eq_mux,v_mux8_1) = 1;
    PSI(eq_mux,e_mux0) = 1;

    GAM0(eq_mux8_1,v_mux8_1) = 1;	
    GAM1(eq_mux8_1,v_mux8_2) = 1;	

    GAM0(eq_mux8_2,v_mux8_2) = 1;	
    GAM1(eq_mux8_2,v_mux8_3) = 1;	

    GAM0(eq_mux8_3,v_mux8_3) = 1;	
    GAM1(eq_mux8_3,v_mux8_4) = 1;	

    GAM0(eq_mux8_4,v_mux8_4) = 1;	
    GAM1(eq_mux8_4,v_mux8_5) = 1;
    PSI(eq_mux8_4,e_mux4)     = 1;

    GAM0(eq_mux8_5,v_mux8_5) = 1;
    GAM1(eq_mux8_5,v_mux8_6) = 1;

    GAM0(eq_mux8_6,v_mux8_6) = 1;
    GAM1(eq_mux8_6,v_mux8_7) = 1;

    GAM0(eq_mux8_7,v_mux8_7) = 1;
    GAM1(eq_mux8_7,v_mux8_8) = 1;

    GAM0(eq_mux8_8,v_mux8_8) = 1;
    PSI(eq_mux8_8,e_mux8)   = 1;


    % zi%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_zi,v_zi)    = 1;
    GAM1(eq_zi,v_zi)    = rhozi;
    GAM1(eq_zi,v_zi8_1) = 1;
    PSI(eq_zi,e_zi0) = 1;

    GAM0(eq_zi8_1,v_zi8_1) = 1;
    GAM1(eq_zi8_1,v_zi8_2) = 1;

    GAM0(eq_zi8_2,v_zi8_2) = 1;
    GAM1(eq_zi8_2,v_zi8_3) = 1;

    GAM0(eq_zi8_3,v_zi8_3) = 1;
    GAM1(eq_zi8_3,v_zi8_4) = 1;

    GAM0(eq_zi8_4,v_zi8_4) = 1;
    GAM1(eq_zi8_4,v_zi8_5) = 1;
    PSI(eq_zi8_4,e_zi4)     = 1;

    GAM0(eq_zi8_5,v_zi8_5) = 1;
    GAM1(eq_zi8_5,v_zi8_6) = 1;

    GAM0(eq_zi8_6,v_zi8_6) = 1;
    GAM1(eq_zi8_6,v_zi8_7) = 1;

    GAM0(eq_zi8_7,v_zi8_7) = 1;
    GAM1(eq_zi8_7,v_zi8_8) = 1;

    GAM0(eq_zi8_8,v_zi8_8) = 1;
    PSI(eq_zi8_8,e_zi8)   = 1;

    % z%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    GAM0(eq_z,v_z)    = 1;
    GAM1(eq_z,v_z)    = rhoz;
    GAM1(eq_z,v_z8_1) = 1;
    PSI(eq_z,e_z0) = 1;

    GAM0(eq_z8_1,v_z8_1) = 1;
    GAM1(eq_z8_1,v_z8_2) = 1;

    GAM0(eq_z8_2,v_z8_2) = 1;
    GAM1(eq_z8_2,v_z8_3) = 1;

    GAM0(eq_z8_3,v_z8_3) = 1;
    GAM1(eq_z8_3,v_z8_4) = 1;

    GAM0(eq_z8_4,v_z8_4) = 1;
    GAM1(eq_z8_4,v_z8_5) = 1;
    PSI(eq_z8_4,e_z4)     = 1;

    GAM0(eq_z8_5,v_z8_5) = 1;
    GAM1(eq_z8_5,v_z8_6) = 1;

    GAM0(eq_z8_6,v_z8_6) = 1;
    GAM1(eq_z8_6,v_z8_7) = 1;

    GAM0(eq_z8_7,v_z8_7) = 1;
    GAM1(eq_z8_7,v_z8_8) = 1;

    GAM0(eq_z8_8,v_z8_8) = 1;
    PSI(eq_z8_8,e_z8)   = 1;

    %mu%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_mu,v_mu)    = 1;
    GAM1(eq_mu,v_mu)    = rhom;
    GAM1(eq_mu,v_mu8_1) = 1;
    PSI(eq_mu,e_mu0) = 1;

    GAM0(eq_mu8_1,v_mu8_1) = 1;
    GAM1(eq_mu8_1,v_mu8_2) = 1;

    GAM0(eq_mu8_2,v_mu8_2) = 1;
    GAM1(eq_mu8_2,v_mu8_3) = 1;

    GAM0(eq_mu8_3,v_mu8_3) = 1;
    GAM1(eq_mu8_3,v_mu8_4) = 1;

    GAM0(eq_mu8_4,v_mu8_4) = 1;
    GAM1(eq_mu8_4,v_mu8_5) = 1;
    PSI(eq_mu8_4,e_mu4)     = 1;

    GAM0(eq_mu8_5,v_mu8_5) = 1;
    GAM1(eq_mu8_5,v_mu8_6) = 1;

    GAM0(eq_mu8_6,v_mu8_6) = 1;
    GAM1(eq_mu8_6,v_mu8_7) = 1;

    GAM0(eq_mu8_7,v_mu8_7) = 1;
    GAM1(eq_mu8_7,v_mu8_8) = 1;

    GAM0(eq_mu8_8,v_mu8_8) = 1;
    PSI(eq_mu8_8,e_mu8)   = 1;

    %g%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_g,v_g)    = 1;
    GAM1(eq_g,v_g)    = rhog;
    GAM1(eq_g,v_g8_1) = 1;
    PSI(eq_g,e_g0)    = 1;

    GAM0(eq_g8_1,v_g8_1) = 1;
    GAM1(eq_g8_1,v_g8_2) = 1;

    GAM0(eq_g8_2,v_g8_2) = 1;
    GAM1(eq_g8_2,v_g8_3) = 1;

    GAM0(eq_g8_3,v_g8_3) = 1;
    GAM1(eq_g8_3,v_g8_4) = 1;

    GAM0(eq_g8_4,v_g8_4) = 1;
    GAM1(eq_g8_4,v_g8_5) = 1;
    PSI(eq_g8_4,e_g4)     = 1;

    GAM0(eq_g8_5,v_g8_5) = 1;
    GAM1(eq_g8_5,v_g8_6) = 1;

    GAM0(eq_g8_6,v_g8_6) = 1;
    GAM1(eq_g8_6,v_g8_7) = 1;

    GAM0(eq_g8_7,v_g8_7) = 1;
    GAM1(eq_g8_7,v_g8_8) = 1;

    GAM0(eq_g8_8,v_g8_8) = 1;
    PSI(eq_g8_8,e_g8)   = 1;

    %zet%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_zet,v_zet)    = 1;
    GAM1(eq_zet,v_zet)    = rhozet;
    GAM1(eq_zet,v_zet8_1) = 1;
    PSI(eq_zet,e_zet0) = 1;

    GAM0(eq_zet8_1,v_zet8_1) = 1;
    GAM1(eq_zet8_1,v_zet8_2) = 1;

    GAM0(eq_zet8_2,v_zet8_2) = 1;
    GAM1(eq_zet8_2,v_zet8_3) = 1;

    GAM0(eq_zet8_3,v_zet8_3) = 1;
    GAM1(eq_zet8_3,v_zet8_4) = 1;

    GAM0(eq_zet8_4,v_zet8_4) = 1;
    GAM1(eq_zet8_4,v_zet8_5) = 1;
    PSI(eq_zet8_4,e_zet4)     = 1;

    GAM0(eq_zet8_5,v_zet8_5) = 1;
    GAM1(eq_zet8_5,v_zet8_6) = 1;

    GAM0(eq_zet8_6,v_zet8_6) = 1;
    GAM1(eq_zet8_6,v_zet8_7) = 1;

    GAM0(eq_zet8_7,v_zet8_7) = 1;
    GAM1(eq_zet8_7,v_zet8_8) = 1;

    GAM0(eq_zet8_8,v_zet8_8) = 1;
    PSI(eq_zet8_8,e_zet8)   = 1;
    
 elseif msel==2 %MA(1) shocks
            % MUA%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    GAM0(eq_mua,v_mua)    = 1;
    GAM0(eq_mua,v_mua8_1) = -mamua0;
    GAM1(eq_mua,v_mua)    = rhoa;
    GAM1(eq_mua,v_mua8_4)    = mamua4;
    GAM1(eq_mua,v_mua8_8)    = mamua8;
    
    %Create MA(8) auxiliary states:

    GAM0(eq_mua8_1,v_mua8_1) = 1;
    PSI(eq_mua8_1,e_mua0) = 1;

    GAM0(eq_mua8_2,v_mua8_2) = 1;
    GAM1(eq_mua8_2,v_mua8_1) = 1;

    GAM0(eq_mua8_3,v_mua8_3) = 1;
    GAM1(eq_mua8_3,v_mua8_2) = 1;

    GAM0(eq_mua8_4,v_mua8_4) = 1;
    GAM1(eq_mua8_4,v_mua8_3) = 1;


    GAM0(eq_mua8_5,v_mua8_5) = 1;
    GAM1(eq_mua8_5,v_mua8_4) = 1;

    GAM0(eq_mua8_6,v_mua8_6) = 1;
    GAM1(eq_mua8_6,v_mua8_5) = 1;

    GAM0(eq_mua8_7,v_mua8_7) = 1;
    GAM1(eq_mua8_7,v_mua8_6) = 1;

    GAM0(eq_mua8_8,v_mua8_8) = 1;
    GAM1(eq_mua8_8,v_mua8_7) = 1;
   

    % MUX%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    GAM0(eq_mux,v_mux)    = 1;
    GAM0(eq_mux,v_mux8_1) = -mamux0;
    GAM1(eq_mux,v_mux)    = rhox;
    GAM1(eq_mux,v_mux8_4)    = mamux4;
    GAM1(eq_mux,v_mux8_8)    = mamux8;

    
    %Create MA(8) auxiliary states:

    GAM0(eq_mux8_1,v_mux8_1) = 1;
    PSI(eq_mux8_1,e_mux0) = 1;

    GAM0(eq_mux8_2,v_mux8_2) = 1;
    GAM1(eq_mux8_2,v_mux8_1) = 1;

    GAM0(eq_mux8_3,v_mux8_3) = 1;
    GAM1(eq_mux8_3,v_mux8_2) = 1;

    GAM0(eq_mux8_4,v_mux8_4) = 1;
    GAM1(eq_mux8_4,v_mux8_3) = 1;


    GAM0(eq_mux8_5,v_mux8_5) = 1;
    GAM1(eq_mux8_5,v_mux8_4) = 1;

    GAM0(eq_mux8_6,v_mux8_6) = 1;
    GAM1(eq_mux8_6,v_mux8_5) = 1;

    GAM0(eq_mux8_7,v_mux8_7) = 1;
    GAM1(eq_mux8_7,v_mux8_6) = 1;

    GAM0(eq_mux8_8,v_mux8_8) = 1;
    GAM1(eq_mux8_8,v_mux8_7) = 1;


    % zi%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_zi,v_zi)    = 1;
    GAM0(eq_zi,v_zi8_1) = -mazi0;
    GAM1(eq_zi,v_zi)    = rhozi;
    GAM1(eq_zi,v_zi8_4)    = mazi4;
    GAM1(eq_zi,v_zi8_8)    = mazi8;
  
    
    %Create MA(8) auxiliary states:

    GAM0(eq_zi8_1,v_zi8_1) = 1;
    PSI(eq_zi8_1,e_zi0) = 1;

    GAM0(eq_zi8_2,v_zi8_2) = 1;
    GAM1(eq_zi8_2,v_zi8_1) = 1;

    GAM0(eq_zi8_3,v_zi8_3) = 1;
    GAM1(eq_zi8_3,v_zi8_2) = 1;

    GAM0(eq_zi8_4,v_zi8_4) = 1;
    GAM1(eq_zi8_4,v_zi8_3) = 1;


    GAM0(eq_zi8_5,v_zi8_5) = 1;
    GAM1(eq_zi8_5,v_zi8_4) = 1;

    GAM0(eq_zi8_6,v_zi8_6) = 1;
    GAM1(eq_zi8_6,v_zi8_5) = 1;

    GAM0(eq_zi8_7,v_zi8_7) = 1;
    GAM1(eq_zi8_7,v_zi8_6) = 1;

    GAM0(eq_zi8_8,v_zi8_8) = 1;
    GAM1(eq_zi8_8,v_zi8_7) = 1;

    % z%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    GAM0(eq_z,v_z)    = 1;
    GAM0(eq_z,v_z8_1) = -maz0;
    GAM1(eq_z,v_z)    = rhoz;
    GAM1(eq_z,v_z8_4)    = maz4;
    GAM1(eq_z,v_z8_8)    = maz8;
 
    
    %Create MA(8) auxiliary states:

    GAM0(eq_z8_1,v_z8_1) = 1;
    PSI(eq_z8_1,e_z0) = 1;

    GAM0(eq_z8_2,v_z8_2) = 1;
    GAM1(eq_z8_2,v_z8_1) = 1;

    GAM0(eq_z8_3,v_z8_3) = 1;
    GAM1(eq_z8_3,v_z8_2) = 1;

    GAM0(eq_z8_4,v_z8_4) = 1;
    GAM1(eq_z8_4,v_z8_3) = 1;


    GAM0(eq_z8_5,v_z8_5) = 1;
    GAM1(eq_z8_5,v_z8_4) = 1;

    GAM0(eq_z8_6,v_z8_6) = 1;
    GAM1(eq_z8_6,v_z8_5) = 1;

    GAM0(eq_z8_7,v_z8_7) = 1;
    GAM1(eq_z8_7,v_z8_6) = 1;

    GAM0(eq_z8_8,v_z8_8) = 1;
    GAM1(eq_z8_8,v_z8_7) = 1;

    %mu%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_mu,v_mu)    = 1;
    GAM0(eq_mu,v_mu8_1) = -mamu0;
    GAM1(eq_mu,v_mu)    = rhom;
    GAM1(eq_mu,v_mu8_4)    = mamu4;
    GAM1(eq_mu,v_mu8_8)    = mamu8;
    
    %Create MA(8) auxiliary states:

    GAM0(eq_mu8_1,v_mu8_1) = 1;
    PSI(eq_mu8_1,e_mu0) = 1;

    GAM0(eq_mu8_2,v_mu8_2) = 1;
    GAM1(eq_mu8_2,v_mu8_1) = 1;

    GAM0(eq_mu8_3,v_mu8_3) = 1;
    GAM1(eq_mu8_3,v_mu8_2) = 1;

    GAM0(eq_mu8_4,v_mu8_4) = 1;
    GAM1(eq_mu8_4,v_mu8_3) = 1;


    GAM0(eq_mu8_5,v_mu8_5) = 1;
    GAM1(eq_mu8_5,v_mu8_4) = 1;

    GAM0(eq_mu8_6,v_mu8_6) = 1;
    GAM1(eq_mu8_6,v_mu8_5) = 1;

    GAM0(eq_mu8_7,v_mu8_7) = 1;
    GAM1(eq_mu8_7,v_mu8_6) = 1;

    GAM0(eq_mu8_8,v_mu8_8) = 1;
    GAM1(eq_mu8_8,v_mu8_7) = 1;

    %g%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_g,v_g)    = 1;
    GAM0(eq_g,v_g8_1) = -mag0;
    GAM1(eq_g,v_g)    = rhog;
    GAM1(eq_g,v_g8_4)    = mag4;
    GAM1(eq_g,v_g8_8)    = mag8;
   
    
    %Create MA(8) auxiliary states:

    GAM0(eq_g8_1,v_g8_1) = 1;
    PSI(eq_g8_1,e_g0) = 1;

    GAM0(eq_g8_2,v_g8_2) = 1;
    GAM1(eq_g8_2,v_g8_1) = 1;

    GAM0(eq_g8_3,v_g8_3) = 1;
    GAM1(eq_g8_3,v_g8_2) = 1;

    GAM0(eq_g8_4,v_g8_4) = 1;
    GAM1(eq_g8_4,v_g8_3) = 1;


    GAM0(eq_g8_5,v_g8_5) = 1;
    GAM1(eq_g8_5,v_g8_4) = 1;

    GAM0(eq_g8_6,v_g8_6) = 1;
    GAM1(eq_g8_6,v_g8_5) = 1;

    GAM0(eq_g8_7,v_g8_7) = 1;
    GAM1(eq_g8_7,v_g8_6) = 1;

    GAM0(eq_g8_8,v_g8_8) = 1;
    GAM1(eq_g8_8,v_g8_7) = 1;

    %zet%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GAM0(eq_zet,v_zet)    = 1;
    GAM0(eq_zet,v_zet8_1) = -mazet0;
    GAM1(eq_zet,v_zet)    = rhozet;
    GAM1(eq_zet,v_zet8_4)    = mazet4;
    GAM1(eq_zet,v_zet8_8)    = mazet8;
  
    
    %Create MA(8) auxiliary states:

    GAM0(eq_zet8_1,v_zet8_1) = 1;
    PSI(eq_zet8_1,e_zet0) = 1;

    GAM0(eq_zet8_2,v_zet8_2) = 1;
    GAM1(eq_zet8_2,v_zet8_1) = 1;

    GAM0(eq_zet8_3,v_zet8_3) = 1;
    GAM1(eq_zet8_3,v_zet8_2) = 1;

    GAM0(eq_zet8_4,v_zet8_4) = 1;
    GAM1(eq_zet8_4,v_zet8_3) = 1;


    GAM0(eq_zet8_5,v_zet8_5) = 1;
    GAM1(eq_zet8_5,v_zet8_4) = 1;

    GAM0(eq_zet8_6,v_zet8_6) = 1;
    GAM1(eq_zet8_6,v_zet8_5) = 1;

    GAM0(eq_zet8_7,v_zet8_7) = 1;
    GAM1(eq_zet8_7,v_zet8_6) = 1;

    GAM0(eq_zet8_8,v_zet8_8) = 1;
    GAM1(eq_zet8_8,v_zet8_7) = 1;
 end

    %% 7. Expectation terms                               
    GAM0(eq_Elam, v_lam)  = 1;
    GAM1(eq_Elam, v_Elam) = 1;
    PPI(eq_Elam, sh_Elam) = 1;

    GAM0(eq_Emuy, v_muy)  = 1;
    GAM1(eq_Emuy, v_Emuy) = 1;
    PPI(eq_Emuy, sh_Emuy) = 1;

    GAM0(eq_Emuk, v_muk)  = 1;
    GAM1(eq_Emuk, v_Emuk) = 1;
    PPI(eq_Emuk, sh_Emuk) = 1;

    GAM0(eq_Ei, v_i)  = 1;
    GAM1(eq_Ei, v_Ei) = 1;
    PPI(eq_Ei, sh_Ei) = 1;

    GAM0(eq_Eu, v_u)  = 1;
    GAM1(eq_Eu, v_Eu) = 1;
    PPI(eq_Eu, sh_Eu) = 1;

    GAM0(eq_Ez, v_z)  = 1;
    GAM1(eq_Ez, v_Ez) = 1;
    PPI(eq_Ez, sh_Ez) = 1;

    GAM0(eq_Eh, v_h)  = 1;
    GAM1(eq_Eh, v_Eh) = 1;
    PPI(eq_Eh, sh_Eh) = 1;

    GAM0(eq_Emua, v_mua)  = 1;
    GAM1(eq_Emua, v_Emua) = 1;
    PPI(eq_Emua, sh_Emua) = 1;

    GAM0(eq_Eq, v_q)  = 1;
    GAM1(eq_Eq, v_Eq) = 1;
    PPI(eq_Eq, sh_Eq) = 1;

    GAM0(eq_Ep, v_p)  = 1;
    GAM1(eq_Ep, v_Ep) = 1;
    PPI(eq_Ep, sh_Ep) = 1;

    GAM0(eq_Es, v_s)  = 1;
    GAM1(eq_Es, v_Es) = 1;
    PPI(eq_Es, sh_Es) = 1;

    GAM0(eq_Ezet, v_zet)  = 1;
    GAM1(eq_Ezet, v_Ezet) = 1;
    PPI(eq_Ezet, sh_Ezet) = 1;

    GAM0(eq_Ec, v_c)  = 1;
    GAM1(eq_Ec, v_Ec) = 1;
    PPI(eq_Ec, sh_Ec) = 1;

    GAM0(eq_Ev, v_v)  = 1;
    GAM1(eq_Ev, v_Ev) = 1;
    PPI(eq_Ev, sh_Ev) = 1;


    %%   Lags
% Used for constructing observables to simplify code.
% Helpful for estimation using Kalman filter in the time domain.
%     GAM0(eq_ylag, v_ylag) = 1;
%     GAM1(eq_ylag, v_y)    = 1;
% 
%     GAM0(eq_clag, v_clag) = 1;
%     GAM1(eq_clag, v_c)    = 1;
% 
%     GAM0(eq_ilag, v_ilag) = 1;
%     GAM1(eq_ilag, v_i)    = 1;
% 
%     GAM0(eq_glag, v_glag) = 1;
%     GAM1(eq_glag, v_g)    = 1;
% 
%     GAM0(eq_xglag, v_xglag) = 1;
%     GAM1(eq_xglag, v_xg)    = 1;
% 
%     GAM0(eq_hlag, v_hlag) = 1;
%     GAM1(eq_hlag, v_h)    = 1;
% 
%     GAM0(eq_zlag, v_zlag) = 1;
%     GAM1(eq_zlag, v_z)    = 1; 