% Update 17 Nov. 2025 by Tong
% prepare workspace
clear; close all;

% load default parameters
run('./par_MtAp_default.m')

% set run parameters
runID     =  '1D_MtAp';           % run identifier
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nrh       =  1e2;                 % record diagnostic history every 'nrh' time steps
nop       =  1e3;                 % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
colourmap = 'lapaz';              % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

% set model domain parameters
D         =  100;                 % chamber depth [m]
N         =  200;                 % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  h;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  1e6;                 % number of time steps to take
tend      =  100*yr;              % end time for simulation [s]
dt        =  1;                   % initial time step [s]

% set initial thermo-chemical state
init_mode =  'liquidus';          % auto calculate liquidus
T0        =  -50;                  % ? initial temperature [deg C]
c0        =  [15.4159   10.0390   15.4879   19.3207   39.7364  6  2]/100;  % *** components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*0e-4;
dr_trc    =  [1,1,1,-1,-1,-1  ]*0e-4; % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;                   % periodic side boundaries
bndmode   =  3;                   % *** boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  (2*h)^2/1e-6;        % wall cooling/assimilation time [s]
Twall     =  [500,500,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop      =  1.25e8;              % *** top pressure [Pa]

% set thermo-chemical material parameters
calID     =  'MtAp_750_new';      % *** phase diagram calibration

% set effective diffusivity parameters
Delta_cnv =  4*h;                 % correlation length for eddy, convection diffusivity (multiple of h, 0.5-1)
Delta_sgr =  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-4;                % outer its relative tolerance
atol      =  1e-9;                % outer its absolute tolerance
maxit     =  20;                  % maximum outer its
itpar.fp.damp = 1;                % fixed-point iterative damping (0-1)
itpar.aa.m    = 4;                % Anderson acceleration depth (2-5)
itpar.aa.damp = 0.5;              % Anderson acceleration damping (0-1)
itpar.aa.reg  = 0.01;             % Anderson acceleration regularisation (0-1)

%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

