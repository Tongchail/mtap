function etaf = etamfe(TC)

% ETA_MFE  viscosity of Fe-rich melt (mfe composition) in Pa·s
%   T : temperature in deg C 
%   eta : viscosity in Pa·s 

R   = 8.3145;  % universal gas constant


% read in analysed experimental compositions
%         SiO2       TiO2     Al2O3      FeO(t)    MgO       CaO       Na2O      K2O       P2O5
OX   =  [ 70.4500    0.3000   14.5600    2.6400    0.1700    0.3200    3.3500    7.7800         0 ;  % rhyolite (synthetic)
          64.2300    0.8700   17.2500    3.1500    1.0000    5.6100    4.5200    3.1400    0.1300 ;  % andesite (El Laco sample)
          39.5800    5.4600    1.6800   21.5700   12.8800   12.6100    0.6400    0.1800    2.0100 ;  % Fe-cpx   (synthetic)
          32.6154    4.5000    1.5192   36.3654   10.8558    9.7500    0.5481         0    1.9231 ]; % ***mfe  (new)

       
% combine major oxides into three end-members (EMs)
ind1 = [2 4 5 6 9];  % TiO2, FeOtot, MgO, CaO, P2O5 (Fe-rich EM)
ind2 = [3 7 8];      % Al2O3, Na2O, K2O (Al-Na-K-rich EM)
ind3 = [1];          % SiO2 (Silica EM)

EM  =  [sum(OX(:,ind1),2) sum(OX(:,ind2),2) sum(OX(:,ind3),2)];
EM  =  EM./sum(EM,2); % normalise to unit sum


% read in measured data
Tmp_cpx = [1483 1459 1435 1411 1387 1363 1339 1315];
Eta_cpx = [0.62 0.65 0.69 0.72 0.77 0.82 0.89 0.97];

Tmp_and = [ 1483  1459  1435  1411  1387  1363  1339   1315   1291   1267   1244   1220   1196    1172    1148    1124    1100];
Eta_and = [12.00 15.61 20.68 27.57 37.23 51.01 71.14 100.39 142.86 207.66 307.46 464.98 718.76 1137.84 1860.86 3127.40 5379.65];

Tmp_rhy = [ 1627  1603  1579  1555   1531   1507   1483   1459   1435   1411    1387    1363    1339    1315    1291];
Eta_rhy = [40.51 54.24 73.15 99.44 136.20 188.78 266.62 377.73 546.40 801.64 1181.70 1763.62 2677.91 4117.14 6373.63];


% fit T-dependence of viscosity eta to obtain prefactors A0, activation energies Ea
C = [ones(length(Tmp_rhy),1),1./R./(Tmp_rhy.'+273.15)];
Eta_ft = C\log(Eta_rhy).';   A0(1) = exp(Eta_ft(1)); Ea(1) = Eta_ft(2);
C = [ones(length(Tmp_and),1),1./R./(Tmp_and.'+273.15)];
Eta_ft = C\log(Eta_and).';   A0(2) = exp(Eta_ft(1)); Ea(2) = Eta_ft(2);
C = [ones(length(Tmp_cpx),1),1./R./(Tmp_cpx.'+273.15)];
Eta_ft = C\log(Eta_cpx).';   A0(3) = exp(Eta_ft(1)); Ea(3) = Eta_ft(2);


% fit C-dependence of Ea, A0
C = [ones(3,1),EM(1:3,1),EM(1:3,2)];

Ea_ft = C\log(Ea).';
A0_ft = C\log(A0).';

Ea(4) = exp(Ea_ft(1) + Ea_ft(2).*EM(4,1) + Ea_ft(3).*EM(4,2));
A0(4) = exp(A0_ft(1) + A0_ft(2).*EM(4,1) + A0_ft(3).*EM(4,2));

etaf = A0(4) .* exp(Ea(4)./R./(TC+273.15)); %[Pas]
end