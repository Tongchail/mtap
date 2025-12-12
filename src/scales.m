%*****  calculate and print characteristic scales  ************************

% length scales
D0      =  D/20;
dx0     =  dx0;
df0     =  df0;
L0      =  L0;  L0h  = (L0 +h)/2;
l0x     =  l0x; l0xh = (l0x+h)/2;
l0f     =  l0f; l0fh = (l0f+h)/2;
h0      =  h;

c0_oxd  = c0 *cal.cmp_oxd;

if cal.noxd>9
     c0_oxd_all =  c0_oxd(cal.ioxd);
else
     c0_oxd_all(cal.ioxd) =  c0_oxd;
end

var.c      = c0;                          % component fractions [wt]
var.T      = 1200;                        % temperature [C]
var.P      = Ptop/1e9;                    % pressure [GPa]
var.x      = 0;
var.m      = 1;                           % melt fraction [wt]
var.f      = 0;                           % immiscible mfe fraction [wt]
var.H2O    = c0(:,end);                   % water concentration [wt]
var.MFE    = c0(:,end-1);                 % ** mfe concentration [wt]
var.X      = c0_oxd_all;                  % melt oxide fractions [wt %]
cal.H2Osat = fluidsat(var);               % water saturation [wt]
cal.MFEsat = mfesat(var.T,cal);           % ** mfe saturation [wt] T[K] 

if strcmp(init_mode,'liquidus')          
    [var,cal]  = leappartmfe(var,cal,'T');
    var.T      = cal.Tliq-T0;
    T00        = var.T + 273.15;
    cal.MFEsat = mfesat(var.T,cal);           % ** mfe saturation [wt] T[K] 
else
    T00        = T0;
    var.T      = T0;
    cal.MFEsat = mfesat(var.T,cal);           % ** mfe saturation [wt] T[K] 
end

[var,cal]  = leappartmfe(var,cal,'E');
m0  = var.m;
x0  = var.x;
f0  = var.f;

cm0 = var.cm;
cx0 = var.cx;
cf0 = var.cf;

cal = rmfield(cal,'Tsol');
cal = rmfield(cal,'Tliq');

cm0_oxd = cm0*cal.cmp_oxd;
cx0_oxd = cx0*cal.cmp_oxd;
cf0_oxd = cf0*cal.cmp_oxd;

if cal.noxd>9
    cf0_oxd_all = cf0_oxd(cal.ioxd);
    cm0_oxd_all = cm0_oxd(cal.ioxd);
    cx0_oxd_all = cx0_oxd(cal.ioxd);
else
    cf0_oxd_all(cal.ioxd) = cf0_oxd;
    cm0_oxd_all(cal.ioxd) = cm0_oxd;
    cx0_oxd_all(cal.ioxd) = cx0_oxd;
end

cx0_mem = cx0*cal.cmp_mem;

% material parameter scales
rhom0  = DensityX(cm0_oxd_all,Tref,Pref./1e8);
rhox0  = sum(cx0_mem/100./cal.rhox0,2).^-1;
rhof0  = DensityX(cf0_oxd_all,Tref,Pref./1e8);

rhom0  = rhom0 .* (1 - aTm.*(T00-Tref) + bPm.*(Ptop-Pref));
rhox0  = rhox0 .* (1 - aTx.*(T00-Tref) + bPx.*(Ptop-Pref));
rhof0  = rhof0 .* (1 - aTf.*(T00-Tref) + bPf.*(Ptop-Pref));

rho0   = 1./(m0./rhom0 + x0./rhox0 + f0./rhof0);

etam0  = Giordano08(cm0_oxd_all,T00-273.15);
etaf0  =     etamfe(cf0_oxd_all,T00-273.15);

Drho_x0 =  abs(rhox0-rho0);
Drho_f0 =  abs(rhof0-rho0);
Dchi0   =  x0/20 + 0.0005;
chi0    =  x0 + 0.01;
Dphi0   =  f0/20 + 0.0005;
phi0    =  f0 + 0.01;
Drho0   =  Dchi0*Drho_x0 + Dphi0*Drho_f0;
eta0    =  etam0;

% speed scales
W0l     =  Drho0  *g0*D0 ^2/eta0;  % laminar convection speed
w0l_x   =  Drho_x0*g0*dx0^2/eta0;  % laminar settling speed
w0l_f   =  Drho_f0*g0*df0^2/eta0;  % laminar settling speed


ReL0    =  W0l  *L0/(eta0/rho0);        % laminar convective Reynolds No at L0
Rel0_x  =  w0l_x*l0x/(eta0/rho0);        % laminar settling Reynolds No at L0
Rel0_f  =  w0l_f*l0f/(eta0/rho0);        % laminar settling Reynolds No at L0

digits   =  16;
fReL0    =  vpa(1-exp(-ReL0  ));         % Re-dependent ramp factor
fRel0_x  =  vpa(1-exp(-Rel0_x));         % Re-dependent ramp factor
fRel0_f  =  vpa(1-exp(-Rel0_f));         % Re-dependent ramp factor

W0t  =  sqrt(2*Drho0*g0*D0^3/(L0^2*rho0));  % terminal turbulent convective speed
W0i  =  sqrt(D*Drho0*g0/rho0);              % inertially limited convective speed

if open_cnv
    Ri0 = 1;
else
    Ri0 = W0i./W0t;
end

% general convective speed
W0  = double(D0.*(sqrt(2.*Drho0.*g0.*rho0.*Ri0.^-2.*fReL0.*L0.^2.*D0  + eta0.^2) - eta0)./(Ri0.^-2.*fReL0.*L0.^2.*rho0));

% general settling speed
wx0 = double((sqrt(4.*Drho_x0.*g0.*rho0.*fRel0_x.*l0x.*dx0.^2 + eta0.^2) - eta0)./(2.*fRel0_x.*l0x.*rho0));
wf0 = double((sqrt(4.*Drho_f0.*g0.*rho0.*fRel0_f.*l0f.*df0.^2 + eta0.^2) - eta0)./(2.*fRel0_f.*l0f.*rho0));

% diffusivities
eII0    =  W0/D0/2;
ke0     =  eII0*L0^2;
ks_x0   =  wx0*l0x;
ks_f0   =  wf0*l0f;
kx0     =  double(ks_x0 + fReL0.*ke0);
kf0     =  double(ks_f0 + fReL0.*ke0);
kT0     =  (m0*kTm+x0*kTx+f0*kTf)/rho0/(m0*cPm+x0*cPx+f0*cPf);

% times
tW0     =  D./W0;
twx0    =  D./wx0;
twf0    =  D./wf0;
tkx0    =  D.^2./kx0;
tkf0    =  D.^2./kf0;
ti0     =  rho0.*W0/(Drho0.*g0);
t0      =  ti0/2 + min([tW0, twx0, twf0, tkx0, tkf0]);
dt0     =  min([(h0/2)^2./kx0 , (h0/2)^2./kf0, (h0/2)./(W0+wx0), (h0/2)./(W0+wf0)]);

% noise flux amplitudes
xie0    =  double(Xi*sqrt(     fReL0*ke0  /(L0h /W0 )*(L0h /(L0h +h0))^3));
xiex0   =  double(Xi*sqrt(chi0*fReL0*ke0  /(L0h /W0 )*(L0h /(L0h +h0))^3));
xief0   =  double(Xi*sqrt(phi0*fReL0*ke0  /(L0h /W0 )*(L0h /(L0h +h0))^3));
xisx0   =         Xi*sqrt(chi0*     ks_x0/(l0xh/wx0)*(l0xh/(l0xh+h0))^3);
xisf0   =         Xi*sqrt(phi0*     ks_f0/(l0fh/wf0)*(l0fh/(l0fh+h0))^3);
xix0    =  xisx0 + xiex0;
xif0    =  xisf0 + xief0;


% phase change rate
tau0    =  h0/(W0 + max(wx0, wf0)) + dt0;
G0      =  chi0*rho0./tau0;

% viscosities, stress/pressure
etae0   =  double(fReL0  *ke0  *rho0);
etas_x0 =  double(fRel0_x*ks_x0*rho0);
etas_f0 =  double(fRel0_f*ks_f0*rho0);
p0      =  (eta0+etae0)*eII0;

fReL0    =  double(fReL0);
fRel0_x  =  double(fRel0_x);
fRel0_f  =  double(fRel0_f);

% general dimensionless numbers
Da0     =  G0/(rho0/t0);                % Dahmköhler number
Ne0     =  xie0/W0;                     % Eddy noise number
Nsx0    =  xix0/wx0;                   % Particle noise number
Nsf0    =  xif0/wf0;                   % Particle noise number
Rc_x0   =  W0/wx0;                     % Convection number
Rc_f0   =  W0/wf0;                     % Convection number
Ra_x0   =  W0*D0/(kx0+kf0+kT0);             % Rayleigh number
Ra_f0   =  W0*D0/kf0;                       % Rayleigh number
ReD0    =  W0*D0/((eta0+etae0)/rho0);       % Convection Reynolds number
Red_x0  =  wx0*dx0/((eta0+etas_x0)/rho0);   % Particle Reynolds number
Red_f0  =  wf0*df0/((eta0+etas_f0)/rho0);   % Droplet Reynolds number

% print scaling analysis to standard output
fprintf(1,'\n  Scaled domain depth D0    = %1.0e [m]',D0);
fprintf(1,'\n  Part. size       dx0    = %1.0e [m]',dx0);
fprintf(1,'\n  Fluid size         df0    = %1.0e [m]',df0);
fprintf(1,'\n  Eddy  corrl. length L0    = %1.0e [m]',L0);
fprintf(1,'\n  Segr. corrl. length l0_x  = %1.0e [m]',l0x);
fprintf(1,'\n  Segr. corrl. length l0_f  = %1.0e [m]\n',l0f);

fprintf(1,'\n  Density             rho0  = %1.0f  [kg/m3]',rho0);
fprintf(1,'\n  Density dff. xtal   Drho0 = %1.0f   [kg/m3]',Drho_x0);
fprintf(1,'\n  Density dff. fluid  Drho0 = %1.0f   [kg/m3]',Drho_f0);
fprintf(1,'\n  Density dff. bulk   Drho0 = %1.2f  [kg/m3]',Drho0);
fprintf(1,'\n  Xtal frc contrast   Dchi0 = %1.3f [wt]',Dchi0);
fprintf(1,'\n  Fluid frc contrast  Dphi0 = %1.3f [wt]',Dphi0);
fprintf(1,'\n  Viscosity           eta0  = %1.0e [Pas]\n',eta0);

fprintf(1,'\n  Convection        speed  W0  = %1.2e [m/s]',W0);
fprintf(1,'\n  Part. segregation speed  wx0 = %1.2e [m/s]',wx0);
fprintf(1,'\n  Fluid segregation speed  wf0 = %1.2e [m/s]\n',wf0);

fprintf(1,'\n  Thermal     diffusivity  kT0     = %1.1e [m2/s]',kT0);
fprintf(1,'\n  Eddy        diffusivity  ke0     = %1.1e [m2/s]',ke0);
fprintf(1,'\n  Part. segr. diffusivity  ks_x0   = %1.1e [m2/s]',ks_x0);
fprintf(1,'\n  Fluid segr. diffusivity  ks_f0   = %1.1e [m2/s]',ks_f0);
fprintf(1,'\n  Eddy  viscosity          etae    = %1.1e [Pas]',etae0);
fprintf(1,'\n  Segr. viscosity          etas_x0 = %1.1e [Pas]',etas_x0);
fprintf(1,'\n  Segr. viscosity          etas_f0 = %1.1e [Pas]\n',etas_f0);

fprintf(1,'\n  Eddy        noise rate   xie0 = %1.2e [m/s]',xie0);
fprintf(1,'\n  Part. segr. noise rate   xix0 = %1.2e [m/s]',xix0);
fprintf(1,'\n  Fluid segr. noise rate   xif0 = %1.2e [m/s]\n',xif0);

fprintf(1,'\n  Inertial    time    ti0  = %1.2e [s]',ti0);
fprintf(1,'\n  Convection  time    tW0  = %1.2e [s]',tW0);
fprintf(1,'\n  Part. segr. time    twx0 = %1.2e [s]',twx0);
fprintf(1,'\n  Fluid segr. time    twf0 = %1.2e [s]',twf0);
fprintf(1,'\n  Diffusion   time    tkx0 = %1.2e [s]',tkx0);
fprintf(1,'\n  Diffusion   time    tkf0 = %1.2e [s]\n',tkf0);

fprintf(1,'\n  Dahmköhler No         Da0  = %1.2e [1]',Da0);
fprintf(1,'\n  Eddy Noise No         Ne0  = %1.2e [1]',Ne0);
fprintf(1,'\n  Part. Segr. Noise No  Nsx0 = %1.2e [1]',Nsx0);
fprintf(1,'\n  Fluid Segr. Noise No  Nsf0 = %1.2e [1]\n',Nsf0);


fprintf(1,'\n  Convection No       Rc_x0 = %1.2e [1]',Rc_x0);
fprintf(1,'\n  Convection No       Rc_f0 = %1.2e [1]',Rc_f0);
fprintf(1,'\n  Rayleigh No         Ra_x0 = %1.2e [1]',Ra_x0);
fprintf(1,'\n  Rayleigh No         Ra_f0 = %1.2e [1]',Ra_f0);
fprintf(1,'\n  Domain  Reynolds No ReD0  = %1.2e [1]',ReD0);
fprintf(1,'\n  Part. Reynolds No Red_x0  = %1.2e [1]',Red_x0);
fprintf(1,'\n  Part. Reynolds No Red_f0  = %1.2e [1]\n\n\n',Red_f0);


% adjust scales and units for visualisation
if ndm_op
    tsc = t0;  tun = '1';
    Wsc = W0;  Wun = '1';
    wxsc = max(wx0, wf0);  wun = '1';
    wmsc = (wx0*chi0 + wf0*phi0) /(1 - chi0 - phi0);
    whsc = wx0; %?
    psc = p0;  pun = '1';
    kssc = max(ks_x0, ks_f0);  kun = '1';
    kesc = ke0; 
    kxsc = kx0;
    kfsc = kf0;
    xiesc = xie0;  xieun = '1';
    xixsc = xix0;  xixun = '1';
    xifsc = xif0;  xifun = '1';
    esc   = eta0;  eun   = '1';
    eesc  = etae0; 
    essc_x = etas_x0;   % crystal drag 
    essc_f = etas_f0;   % fluid drag 
    rsc   = rho0;  dun   = '1';
    MFSsc = rho0/t0; MFSun = '1';
    xsc   = chi0;  xun = '1';
    Gsc   = rho0/t0;  Gun = '1';
    ssc   = D;  sun = '1';
    Rasc_x  = Ra_x0;
    Rasc_f  = Ra_f0;
    ReDsc = ReD0;
    Redsc_x = Red_x0;  
    Redsc_f = Red_f0;  
    Rcsc_x  = Rc_x0; 
    Rcsc_f  = Rc_f0; 
else
    kssc = 1;  kun = 'm$^2$/s';
    kesc = 1;  
    kxsc = 1; 
    kfsc = 1;
    esc  = 0;  eun   = 'Pas';
    eesc = 1;
    essc_x  = 1; essc_f  = 1;
    rsc   = 1;  dun   = 'kg/m$^3$';
    MFSsc = 1; MFSun = 'kg/m$^3$/s';
    xsc   = 1/100;  xun = 'wt \%';
    Gsc   = 1;  Gun = 'kg/m$^3$/s';
    Rasc_x = 1;
    Rasc_f = 1;
    ReDsc = 1;
    Redsc_x = 1;
    Redsc_f = 1;
    Rcsc_f = 1;
    Rcsc_x = 1;
if t0 < 1e3
    tsc = 1;
    tun = 's';
elseif t0>= 1e3 && t0 < 1e3*hr
    tsc = hr;
    tun = 'hr';
elseif t0 >= 1e3*hr && t0 < 1e2*yr
    tsc = yr;
    tun = 'yr';
elseif t0 >= 1e2*yr
    tsc = 1e3*yr;
    tun = 'kyr';
end
if D < 1e3
    ssc = 1;
    sun = 'm';
elseif D >= 1e3 && D < 1e6
    ssc = 1e3;
    sun = 'km';
else
    ssc = 1e6;
    sun = 'Mm';
end
if W0 < 1000/yr
    Wsc = 1/yr;
    Wun = 'm/yr';
elseif all(W0 >= 1000/yr) && all(W0 < 1000/hr)
    Wsc = 1/hr;
    Wun = 'm/hr';
elseif W0 >= 1000/hr
    Wsc = 1;
    Wun = 'm/s';
end
if wx0 < 1000/yr
    wxsc = 1/yr;
    wun  = 'm/yr';
elseif wx0 >= 1000/yr && wx0 < 1000/hr
    wxsc = 1/hr;
    wun  = 'm/hr';
elseif wx0 >= 1000/hr
    wxsc = 1;
    wun  = 'm/s';
end
wmsc = wxsc;
if p0 < 1e2
    psc = 1;
    pun = 'Pa';
elseif p0 >= 1e3 && p0 < 1e6
    psc = 1e3;
    pun = 'kPa';
elseif p0 >= 1e6 && p0 < 1e9
    psc = 1e6;
    pun = 'MPa';
else
    psc = 1e9;
    pun = 'GPa';
end
xiesc = Wsc;  xieun = Wun;
xixsc = Wsc;  xixun = Wun;
xifsc = Wsc;  xifun = Wun;
whsc  = Wsc;  
end
Xsc = Xc./ssc;
Zsc = Zc./ssc;
Zsf = Zf./ssc;
