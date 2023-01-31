  %% This the main function that computes impulse responses on levels 
  % for the SGU model and a variant with MA shocks for wage markups
  
        %% Setup
        clear      
        
        parsel=[1:34]; 
        
        
        bc=0; %0 for full spectrum; 1 for business cycle frequencies only
        
        msel1=1; %for sgu model
        msel2=6; %MA shock for wage markup

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
      
        %% impulse responses for the orginal SGU model        
       
        [TT,TC,TEPS,TETA,RC] = sgusolvl(thetab); %,msel1);

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
        end

        %% compute impulse response for the origina sgu model
        maxorder = 20;
        irf=impulse(RR,TT,A,0,maxorder);

        %% impulse response with MA wage mark up shocks
        load('Sgul_ma2_msel_6.mat');
        parama = xestpso2;
        [TT,TC,TEPS,TETA,RC] = sgusolvl4(parama,msel2); %solve model using Sims algorithm

        neq=size(TT,2); %number of equations
        ny=7;
        RR = TEPS;
        irfma=impulse(RR,TT,A,0,maxorder);
        indema = 4*3+1;%inde;%1*3+1;

        %% plotting results
        figure
        tiledlayout(2,2);

        
        
        %output response
        ax1 = nexttile;
        indv = 1; %variable index
        inde = 4*3+2; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*4.79,'red','LineWidth',0.5); %fourth horizon anticipated
        ylim([-1 0.1]);
        title(ax1,'Output')
        hold on;
        inde = 4*3+1; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.5,'blue'); %zero horizon
        hold on;
        inde = 4*3+3; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.51,'magenta'); %eight horizon
        hold on;
        indema = 4*3+1;%inde;%1*3+1;

        %order = 1:1:(maxorder+1);
        pickma = (order-1)*19+indema;
        irfima = irfma(indv,pickma);%/(-4.8512);
        plot (-irfima,'black');
        hold off;
        [irfi; irfima]

        %Consumption response
        ax1 = nexttile;
        indv = 2; %variable index
        inde = 4*3+2; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*4.79,'red','LineWidth',0.5); %fourth horizon anticipated
        ylim([-1 0.1]);
        title(ax1,'Consumption')
        hold on;
        inde = 4*3+1; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.5,'blue'); %zero horizon
        hold on;
        inde = 4*3+3; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.51,'magenta'); %eight horizon
        hold on;
        indema = 4*3+1;%inde;%1*3+1;

        %order = 1:1:(maxorder+1);
        pickma = (order-1)*19+indema;
        irfima = irfma(indv,pickma);%/(-4.8512);
        plot (-irfima,'black');
        hold off;
        [irfi; irfima]

        %Investment response
        ax1 = nexttile;
        indv = 3; %variable index
        inde = 4*3+2; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*4.79,'red','LineWidth',0.5); %fourth horizon anticipated
        ylim([-3 0.1]);
        title(ax1,'Investment')
        hold on;
        inde = 4*3+1; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.5,'blue'); %zero horizon
        hold on;
        inde = 4*3+3; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.51,'magenta'); %eight horizon
        hold on;
        indema = 4*3+1;%inde;%1*3+1;

        %order = 1:1:(maxorder+1);
        pickma = (order-1)*19+indema;
        irfima = irfma(indv,pickma);%/(-4.8512);
        plot (-irfima,'black');
        hold off;
        [irfi; irfima]

        %Hours response
        ax1 = nexttile;
        indv = 4; %variable index
        inde = 4*3+2; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*4.79,'red','LineWidth',0.5); %fourth horizon anticipated
        ylim([-1.5 0.1]);
        title(ax1,'Hours')
        hold on;
        inde = 4*3+1; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.5,'blue'); %zero horizon
        hold on;
        inde = 4*3+3; %shock
        %inde = 5*3+1;
        order = 1:1:(maxorder+1);
        pick = (order-1)*21+inde;
        irfi = irf(indv,pick);
        plot (irfi*0.51,'magenta'); %eight horizon
        hold on;
        indema = 4*3+1;%inde;%1*3+1;

        %order = 1:1:(maxorder+1);
        pickma = (order-1)*19+indema;
        irfima = irfma(indv,pickma);%/(-4.8512);
        plot (-irfima,'black');
        hold off;
        [irfi; irfima]

        saveas(gcf,'impulse.pdf');



