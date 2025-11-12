%% 0-D only use the thermalchemistry not use mechanics; balancing the heat: cooling rate and latent heat  
% Update 11 Nov. 2025

% prepare workspace
clear; close all;

% load default parameters
run('./par_MtAp_default.m')

% set run parameters
runID    =  '0D_MtAp_frc10';     % run identifier  
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nrh      =  1;                   % record diagnostic history every 'nrh' time steps (set to 1 for 0D runs)
nop      =  50;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence
colourmap = 'lapaz';             % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

% set model domain parameters
D        =  1;                   % chamber depth [m] 5km
L        =  D;                   % chamber width [m]
N        =  1;                   % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  1e4;                 % number of time steps to take
tend     =  1e4*hr;              % end time for simulation [s]
dt       =  hr/3;                % initial time step [s]
dtmax    =  hr/3;                % maximum time step [s]
 
% set initial thermo-chemical state
T0       =  1125;                % temperature top  layer [deg C]1290
T1       =  T0;                  % temperature base layer [deg C]
c0       =  [15.4159   10.0390   15.4879   19.3207   39.7364  6  2]/100;  %** components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
c1       =  c0;                  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0,0,0,0];
dcg      =  [0,0,0,0,0,0,0,0];

% set thermo-chemical boundary parameters
fractxtl =  1;                   % fractional crystallisation mode for 0-D (Nz=Nx=1)  1 - fraction crystallization 
fractmlt =  0;                   % fractional melting mode for 0-D (Nz=Nx=1)          0 - do not melt
fractres =  0.10;  %0.05         %** residual fraction for fractionation mode           solid res
dPdT     =  0; %3.00e5;          %** decompression rate for 0D models  
bndmode  =  1;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides) 指定“边界同化（assimilation）模式”——也就是哪些墙面参与热量或化学同化（heat or mass assimilation）。
bnd_w    =  1e16;                % boundary layer width [m]
tau_T    =  D^2/1e-6;            %?? wall cooling/assimilation time [s]       timescale  热量通过岩壁传导出去需要多长时间
Twall    =  [300,300,nan];       %?? [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop     =  1.25e8;              %** top pressure [Pa]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  0;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)

% set thermo-chemical material parameters
calID    =  'MtAp_750_new';              % phase diagram calibration

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
alpha    =  0.75;                % iterative step size parameter
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************



