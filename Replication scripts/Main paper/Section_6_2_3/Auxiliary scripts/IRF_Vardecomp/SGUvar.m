        %% This is the script that computes the variance decompostion for the wage markup shock
        clear      
        
        parsel=[1:34]; 
        
        
        bc=0; %0 for full spectrum; 1 for business cycle frequencies only
        
        msel1=1;
        msel2=6;
        %% Base parameter values:
        thet      = 4.74; %Frisch elasticity of labor supply (when gam=b=0)
        gam	      = 0.0019; %governs magnitude of wealth elasticity of labor supply
        %NB: small values of gamma crucial for importance of news shocks
        kap	      = 9.11; %investment adj. cost parameter
        del2_del1 = 0.34; %ratio of sensitivity of capacity utilization to variations to rental rate of capital (del2)
        %to the steady state capacity utilization
        b	      = 0.91; %habit formation parameter
        rhoxg     = 0.72; %smoothness in trend of gov. spending
        
        rhoa      = 0.48; %nonstationary investment-specific productivity shock AR
        rhox      = 0.38; %neutral productivity shock
        rhozi     = 0.47; %investment-specific productivity shock AR
        rhoz      = 0.92; %productivity shock AR
        rhom      = 0.98; %wage markup shock
        rhog      = 0.96; %gov. spending shock
        rhozet    = 0.17; %preference shock
        
        %standard deviations of the shocks
        sig_mua0 = 0.21; %
        sig_mua4 = 0.16;
        sig_mua8 = 0.16;
        sig_mux0 = 0.38;
        sig_mux4 = 0.08;
        sig_mux8 = 0.1;
        sig_zi0  = 11.72;
        sig_zi4  = 1.93;
        sig_zi8  = 5.50;
        sig_z0   = 0.65;
        sig_z4   = 0.11;
        sig_z8   = 0.09;
        sig_mu0  = 0.5;
        sig_mu4  = 4.79;
        sig_mu8  = 0.51;
        sig_g0   = 0.62;
        sig_g4   = 0.57;
        sig_g8   = 0.37;
        sig_zet0 = 4.03;
        sig_zet4 = 1.89;
        sig_zet8 = 2.21;
        
        thetab=[thet;gam;kap;del2_del1;b;rhoxg;rhoa;rhox;rhozi;rhoz;rhom;rhog;rhozet;sig_mua0;sig_mua4;
            sig_mua8;sig_mux0;sig_mux4;sig_mux8;sig_zi0;sig_zi4;sig_zi8;sig_z0;sig_z4;sig_z8;sig_mu0;sig_mu4;sig_mu8;
            sig_g0; sig_g4;sig_g8;sig_zet0;sig_zet4;sig_zet8]';
        
        thetaa=thetab;
        intsel=[26,27,28];
        thetaa(intsel)=0.1;
      
        % solve the model
        [TT,TC,TEPS,TETA,RC] = sgusolvl4(thetab,msel1);

        neq=size(TT,2); %number of equations

        ny=7;
        
        
        A=zeros(ny,neq); %selection matrix A
        
        %vector of observables: [y; cons;iinv;hrs;g;tfp;100*a]
         RR = TEPS;

        eq_y = 1;
        eq_c = 2;
        eq_i = 3;
        eq_h = 4;
        eq_g = 5;
        eq_z = 6;
        eq_agr = 7;

        %variable indices for selection
        v_i	     = 69;
        v_c	     = 70;
        v_y	     = 71;
        v_h	     = 73;
        v_g	     = 49;
        v_z	     = 31;
        v_mua    = 4;

        %fill selection matrix
        A(eq_y,v_y)=1; 

        A(eq_c,v_c)=1; 

        A(eq_i,v_i)=1; 

        A(eq_h,v_h)=1; 

        A(eq_g,v_g)=1; 

        A(eq_z,v_z)=1; 

        A(eq_agr,v_mua)=100;

        QQ=zeros(21,21);
        for i=14:34
            QQ(i-13,i-13)=(thetab(i))^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %QQ(i-13,i-13)=theta0(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
 
        %% variance decompostion for the original model, wrs to wage markup shocks
        maxorder = 100;
        V=variance_decomp(RR,TT,A,-A,QQ, maxorder);
        DV = diag(V);
        %indv = 1;
        %inde = 4*3+2;
        %inde = 5*3+1;
        %order = 1:1:(maxorder+1);
        QQ(14,14)=0;
        V2=variance_decomp(RR,TT,A,-A,QQ, maxorder);
        DV2 = diag(V2);
        disp ('the percentage of variance explained by the fourth-horizon anticipated wage markup shock')
        disp (1.-(DV2./DV));

        %decompostion for the ma specification

        load('Sgul_ma2_msel_6.mat');
        parama = xestpso2;
        [TT,TC,TEPS,TETA,RC] = sgusolvl4(parama,msel2); %solve model using Sims algorithm

        neq=size(TT,2); %number of equations

        ny=7;
        
        %truespec=zeros(ny*length(x),ny);
        
        %A=zeros(ny,neq); %selection matrix A
        
        RR = TEPS;
        QQ=zeros(21,21);
        for i=14:34
            QQ(i-13,i-13)=(thetab(i))^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
            %QQ(i-13,i-13)=theta0(i)^2; %put variances on the diagonal (dividing by 100 in HS(2014) code)
        end
        %if msel2>1
        intsel = [26,27,28];
        QQ(intsel(1)-13,intsel(1)-13)=1;
        QQ(intsel(2:3)-13,:)=[];
        QQ(:,intsel(2:3)-13)=[];
        %end

        V=variance_decomp(RR,TT,A,-A,QQ, maxorder);
        DV = diag(V);

        QQ(13,13)=0;
        V2=variance_decomp(RR,TT,A,-A,QQ, maxorder);
        DV2 = diag(V2);
        disp ('the percentage of variance explained by the MA wage markup shock')
        disp (1.-(DV2./DV));
