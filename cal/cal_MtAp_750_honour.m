% Modified based on Keller cal.MORB.m     Update   19 Nov.2025
% Used for calculating the start components and running the 0D1D model

% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 10;
cal.nmem   = 16;
cal.nmsy   = 7;
cal.ncmp   = 7;



% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','P$_2$O$_5$','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','P','H'};
cal.memStr = {'ant','alb','san','mmt','tmt','mgt','chy','fhy','hyp','mau','fau','aug','ilm','apt','qtz','wat'};
cal.msyStr = {'plg','mgt','opx','cpx','ilm','apt','qtz'};
cal.cmpStr = {'ano','cnr','ogb','tan','rhy','mFe','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2    3   4   5   6    7   8   10]; % oxide indices for viscosity, density functions

cal.imsy = [1 2 4 7];   % Basaltic tetrahedron for this composition

% oxide composition of mineral end-members
%                 SiO2       TiO2      Al2O3       FeO      MgO      CaO       Na2O         K2O      P2O5     H2O
cal.mem_oxd = [ 47.3600         0   33.8000         0         0   17.3400    1.5000         0         0         0
                66.5200         0   21.4100         0         0    2.9300   11.0500    0.0900         0         0
                64.9300         0   19.6200         0         0    0.6900    3.0500   11.7100         0         0

                      0    0.2200    4.0000   89.0200    6.7600         0         0         0         0         0
                      0    5.3800    1.7300   89.1600    3.7300         0         0         0         0         0
                      0    1.0900    0.8400   94.9300    3.1400         0         0         0         0         0

                52.9200         0    4.0000   12.2000   28.9300    1.9500         0         0         0         0
                53.3200         0    1.8900   18.0300   25.3200    1.4400         0         0         0         0
                55.0000         0    0.9000   15.3700   28.1400    0.5900         0         0         0         0

                53.4200         0    1.7900    7.4700   16.4000   20.1400    0.7800         0         0         0
                54.5600         0    0.6500    9.3900   15.2300   18.8100    1.3600         0         0         0
                55.4400         0    0.3100    8.3500   14.7400   19.0800    2.0800         0         0         0

                      0   51.0322         0   42.5365    6.4313         0         0         0         0         0

                      0         0         0         0         0   56.0000         0         0   42.0000         2

                100.00          0         0         0         0         0         0         0         0         0
                      0         0         0         0         0         0         0         0         0  100.0000];


cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0       % plagioclase
               0  0  0  1  1  1  0  0  0  0  0  0  0  0  0  0       % mgt 
               0  0  0  0  0  0  1  1  1  0  0  0  0  0  0  0       % opx 
               0  0  0  0  0  0  0  0  0  1  1  1  0  0  0  0       % cpx 
               0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0       % ilm 
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0       % apt
               0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];     % qtz 

% mineral end-member composition of melting model components
%                  ant         alb       san   |    mmt	     tmt	    mgt|      chy       fhy       hyp|      mau       fau      aug  |     ilm	    apt      qtz      wat
cal.cmp_mem = [  94.0000    5.0000         0         1         0         0         0         0         0         0         0         0         0         0         0         0    % anorthosite (ano)
                 32.6000   16.6000         0    8.7000       0.5       5.8   25.5000         0         0    10.400         0         0         0         0         0         0    % cpx-norite (cnr)
                 29.0000   14.4000    0.4000       0.5    6.6000    1.0000    1.0000   17.0000         0    4.5000   25.6000         0         0         0         0         0    % opx-gabbro (ogb)
                 10.2000   55.5000    6.4000         0         0    2.2000         0         0    6.4000         0         0   13.0000    6.3000         0         0         0    % trachy-andesite (tan)
                       0   22.5000   43.8000         0         0         0         0         0    1.4000         0         0         0    1.3000    1.0000   30.0000         0    % rhyolite (rhy)
                       0         0         0         0   16.1675         0         0   16.7050         0         0         0   42.7373    7.5968   15.0440         0    1.7495    % mfe
                       0         0         0         0         0         0         0         0         0         0         0         0         0         0         0   100.0000];
cal.cmp_mem = cal.cmp_mem./sum(cal.cmp_mem,2)*100;

% mineral systems composition of melting model components
cal.cmp_msy = cal.cmp_mem*cal.msy_mem.';

% oxide composition of melting model components
cal.cmp_oxd = cal.cmp_mem*cal.mem_oxd./100;

% oxide composition of mineral systems in melting model components
for i=1:cal.ncmp
    for j=1:cal.nmsy
        cal.cmp_msy_oxd(i,j,:) = cal.cmp_mem(i,cal.msy_mem(j,:)==1)*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(cal.cmp_mem(i,cal.msy_mem(j,:)==1)+1e-32);
    end
end

% set pure component melting points T_m^i at P=0
cal.T0  = [1650   1105   1057  1005   880];   

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [ 1e16  1e16   1e16   1e16   1e16];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [ 1  1   1   1   1];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [ 35  3.4  3.2  5.8  5.5];

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx  = 300;
cal.Dsf  = 400;

% specify melting point dependence on H2O
cal.dTH2O = 1400 * 1200./cal.T0;  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;                                  % solidus shift from water content exponent

cal.T_liquidus   = 1125;     % Liquidus temperature (°C)
cal.T_solidus    = 697;      % Solidus temperature (°C)
cal.MFE_liquidus = 0.07;     % MFE saturation at liquidus (wt fraction)
cal.k            = 0.75;      % Empirical coefficient

% primary and evolved end-member compositions used in calibration
%                ano      cnr       ogb      bta      rhy       mFe     vol
cal.c0      =  [11.5713  15.5442   15.5544  18.4093  35.3908  1.7111  1.8188]/100;
%cal.c1     = [0.001  0.001  0.001  0.001  0.299  0.697  0.024];

cal.c0_oxd  = [57.2396    0.9866   16.1902    6.8329    4.1644    6.9505    3.5368    1.9812    0.2568    1.8610];
%cal.c1_oxd = [76.34  0.16  11.84  2.80   0.71   1.61  4.42  2.12  2.40];

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               ant   alb   san   mmt   tmt   mgt   mhy   fhy   hyp   mau  fau   aug   ilm   apt   qtz   wat
cal.rhox0   = [2680, 2600, 2550, 4500, 4500, 4650, 3410, 3410, 3410, 3470, 3470, 3470, 4700, 3190, 2540, 1000]; % mem ref densities [kg/m3]
cal.rhof0   = 4000;    %1000          %*** fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%               ant   alb   san   mmt   tmt   mgt   mhy   fhy   hyp   mau  fau   aug   ilm   apt   qtz   wat
cal.etax0   = [1e17, 1e17, 1e17, 1e17, 1e17,  1e17, 1e19,1e19,  1e19,1e19, 1e19, 1e19, 1e17, 1e17 ,1e19,1e0];% mem ref viscosities [Pas]
%cal.etaf0   = 1;% 0.1;                %*** fluid viscosity constant [Pas]
cal.Eax     = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA      =[ 0.25, 0.25, 0.35; ...  % permission slopes
               0.25, 0.25, 0.25; ...  % generally numbers between 0 and 1
               0.25, 0.25, 0.25; ];   % increases permission slopes away from step function 

cal.BB      =[ 0.44, 0.18, 0.38; ...  % permission step locations
               0.60, 0.02, 0.38; ...  % 
               0.72, 0.25, 0.03; ];   % 

cal.CC      =[ 0.30, 0.30, 0.30; ... % permission step widths
               0.60, 0.60, 0.12; ... % 
               0.60, 0.12, 0.60; ];  % 

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;