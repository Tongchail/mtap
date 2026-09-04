 % prepare workspace
clear; close all;

% load default parameters
run('./par_MtAp_default.m')

% set run parameters
runID     =  '2D_MtAp_D100';       % run identifier
outdir    =  '/mnt/home/2388600c/sharedscratch/mtap/out5'; %'../out';
restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop       =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op   =  1;                   % switch on to live plot results
save_op   =  1;                   % switch on to save output to file
colourmap = 'lapaz';              % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

% set model domain parameters
D         =  100;                  % chamber depth [m]
N         =  120;                  % number of grid points in z-direction
h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L         =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt        =  1e6;                 % number of time steps to take
tend      =  100*yr;              % end time for simulation [s]
dt        =  10;                  % initial time step [s]

% set initial thermo-chemical state
% simple- make the boundaries of the 2D domain represent the boundaries of the magma chamber.
init_mode =  'liquidus';          % init_mode = 'constant', 'liquidus', 'layer','linear', 'chamber'.
T0        =  -100;                 % initial temperature [deg C] 
%            [ ano  cnr   ogb   tan   rhy  mFe vol(H20)]
c0        =  [12.9  18.8  12.3  16.4  39.6  6   2]/100;  % *** components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
dcr       =  [1,1,1,-1,-1,-1,0]*1e-4;
dr_trc    =  [1,1,0,-1,-1,0  ]*1e-4; % trace elements random noise

% set thermo-chemical boundary parameters
periodic  =  1;                   % periodic side boundaries
bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w     =  h;                   % boundary layer width [m]
tau_T     =  h^2/5e-7*(0.25/h);   % wall cooling time [s]
tau_a     =  tau_T/10;            % wall assimilation time [s]
Twall     =  [500,500,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
Ptop      =  1.25e8;              % top pressure [Pa]

% set physical control parameters
dx0       =  0.001;               % xtal size [mm]
df0       =  0.001;               % droplet size [mm]
L0        =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
l0x       =  dx0*10;              % correlation length for xtal  phase fluctuation diffusivity (multiple of d0, 10-20)
l0f       =  df0*10;              % correlation length for fluid phase fluctuation diffusivity (multiple of d0, 10-20)
Xi        =  0.5;                 % relative amplitude of random noise flux

% set thermo-chemical material parameters
calID     =  'MtAp_750_honour';      % phase diagram calibration

% set numerical model parameters
TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL       =  0.75;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol      =  1e-3;                % outer its relative tolerance
atol      =  1e-6;                % outer its absolute tolerance
maxit     =  15;                  % maximum outer its
itpar.fp.damp = 1;                % fixed-point iterative damping (0-1)
itpar.aa.m    = 4;                % Anderson acceleration depth (2-5)
itpar.aa.damp = 0.3;              % Anderson acceleration damping (0-1)
itpar.aa.reg  = 1e-6;             % Anderson acceleration regularisation (0-1)

%*****  RUN MTAP MODEL  *************************************************
run('../src/main')
%**************************************************************************

