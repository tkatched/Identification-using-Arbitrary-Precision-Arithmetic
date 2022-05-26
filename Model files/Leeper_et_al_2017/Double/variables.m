% variables.m
% Assigns variable indices.

vnames = {  % Endogenous Variables:
            'cs'        % 1 consumption: Savers
            'cn'        % 2 consumption: Non-Savers
            'R'         % 3 nominal interest rate
            'i'         % 4 investment
            'k'         % 5 effective capital
            'v'         % 6 capital utilization rate
            'l'         % 7 labor
            'y'         % 8 output 
            'gc'		% 9 govt consumption
            'c'         %10 aggregate consumtpion
            'q'         %11 Lagragian multiplier for a unit of investment goods				
            'rk'        %12 real return for private k 
            'w'         %13 real wage 
            'pi'        %14 inflation 
            'b'         %15 govt debt
            'sb'        %16 b/y ratio
            'tauk'      %17 tauk
            'taul'      %18 taul
            'tauc'      %19 tauc
            'r'         %20 real interest rate 
            'z'         %21 transfers 
            'mc'        %22 real marginal cost
            'kbar'      %23 private capital
            'lambda'    %24 household Lagragian multiplier from budget constraint
            'Pb'        %25 price of bonds
            'cstar'     %26 consumption in utility function
            'piL'       %27 long-run inflation
            'rL'        %28 long-run interest rate
            'S'         %29 primary surplus
            'rb'        %30 rb defined
            'Tk'        %31 capital tax revenue
            'Tl'        %32 labor tax revenue
            'Tc'        %33 consumption tax revenue
            % Observables
            'cobs'      %34 consumption
            'iobs'      %35 investment
            'gcobs'     %36 govt spending
            'wobs'      %37 wages
            'bobs'      %38 debt
            'Robs'      %39 interest rate
            'Piobs'     %40 inflation
            'Lobs'      %41 hours worked       
            % Expectations:
            'xlambda'   %42 Cs(t+1)
            'xpi'       %43 pi(t+1)
            'xi'        %44 i(t+1)
            'xq'        %45 q(t+1)
            'xrk'       %46 rk(t+1)
            'xtauk'     %47 tauk(t+1)
            'xw'        %48 w(t+1)
            'xPb'       %49 Pb(t+1)
            'xrL'       %50 rL(t+1)
            % Shocks:     
            'ugc'       %51 gc shock 
            'uz'        %52 z shock
            'ua'        %53 general preference shock
            'ub'        %54 investment shock in adjustment costs
            'um'        %55 technology shock
            'ui'        %56 wage markup shock
            'uw'        %57 price markup shock
            'up'        %58 monetary policy shock
          };  
      

      
for j = 1:size(vnames,1)
   jnk(j,:) = strcat({'N'},vnames(j,:));
   eval([char(jnk(j,:)) ' = j;']);      
end



