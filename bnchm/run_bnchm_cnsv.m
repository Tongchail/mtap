% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_MtAp_default')

% test decreasing time step
ATOL = [1e-3,1e-6,1e-9];

for atol = ATOL

    % set run parameters
    runID     =  'bnchm_cnsv';           % run identifier
    restart   =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
    nop       =  20;                  % output frame plotted/saved every 'nop' time steps
    nrh       =  1;
    plot_op   =  1;                   % switch on to live plot results
    save_op   =  0;                   % switch on to save output to file
    plot_cv   =  0;                   % switch on to live plot iterative convergence
    colourmap = 'lapaz';              % choose colourmap ('ocean','lipari','lajolla','lapaz','navia','batlow(W/K)','glasgow')

    % set model domain parameters
    D         =  100;                 % chamber depth [m]
    N         =  50;                 % number of grid points in z-direction
    h         =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
    L         =  D;                   % chamber width (equal to h for 1-D mode) [m]

    % set model timing parameters
    Nt        =  nop;               % number of time steps to take
    tend      =  100*yr;              % end time for simulation [s]
    dt        =  10;                  % initial time step [s]

    % set initial thermo-chemical state
    % simple- make the boundaries of the 2D domain represent the boundaries of the magma chamber.
    init_mode =  'constant';          % init_mode = 'constant', 'liquidus', 'layer','linear', 'chamber'.
    T0        =  1025;                 % initial temperature [deg C]
    c0        =  [16   11   16   19   38  7  2]/100;  % *** components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
    dcr       =  [1,1,1,-1,-1,-1,0]*1e-3;
    dr_trc    =  [1,1,1,-1,-1,-1]*1e-3; % trace elements random noise

    % set thermo-chemical boundary parameters
    periodic  =  1;                   % periodic side boundaries
    bndmode   =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
    bnd_w     =  h;                   % boundary layer width [m]
    tau_T     =  (h/4)^2/1e-6;        % wall cooling/assimilation time [s]
    tau_a     =  tau_T/10;
    Twall     =  [500,500,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
    Ptop      =  1.25e8;              % top pressure [Pa]

    % set physical control parameters
    dx0       =  0.01;
    df0       =  0.01;
    L0        =  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
    l0x       =  dx0*20;              % correlation length for xtal  phase fluctuation diffusivity (multiple of d0, 10-20)
    l0f       =  df0*20;              % correlation length for fluid phase fluctuation diffusivity (multiple of d0, 10-20)
    Xi        =  0.5;                 % relative amplitude of random noise flux

    % set thermo-chemical material parameters
    calID     =  'MtAp_750_new';      % phase diagram calibration
    tau_r     =  0;

    % set numerical model parameters
    TINT      =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
    ADVN      =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
    CFL       =  0.9;                 % (physical) time stepping courant number (multiplies stable step) [0,1]
    rtol      =  atol/1e6;            % outer its relative tolerance
    maxit     =  100;                 % maximum outer its
    itpar.fp.damp = 1.0;              % fixed-point iterative damping (0-1)
    itpar.aa.m    = 5;                % Anderson acceleration depth (2-5)
    itpar.aa.damp = 0.3;              % Anderson acceleration damping (0-1)
    itpar.aa.reg  = 1e-9;             % Anderson acceleration regularisation (0-1)

    % create output directory
    if ~isfolder([outdir,'/',runID])
        mkdir([outdir,'/',runID]);
    end

    % run code
    run('../src/main')

    % plot convergence
    EB = rms(diff(hist.EB(Nt/2:Nt))./diff(hist.time(Nt/2:Nt)));
    EM = rms(diff(hist.EM(Nt/2:Nt))./diff(hist.time(Nt/2:Nt)));
    EX = rms(diff(hist.EX(Nt/2:Nt))./diff(hist.time(Nt/2:Nt)));
    EF = rms(diff(hist.EF(Nt/2:Nt))./diff(hist.time(Nt/2:Nt)));
    ES = rms(diff(hist.ES(Nt/2:Nt))./diff(hist.time(Nt/2:Nt)));
    EC = rms(diff(hist.EC(Nt/2:Nt,:),1,1)./diff(hist.time(Nt/2:Nt)).','all');

    clist = [colororder;[0 0 0]];

    fh24 = figure(24);
    loglog(atol,EB,'s','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    loglog(atol,EM,'o','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    loglog(atol,EX,'d','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    loglog(atol,EF,'^','Color',clist(5,:),'MarkerSize',10,'LineWidth',2);
    loglog(atol,ES,'v','Color',clist(6,:),'MarkerSize',10,'LineWidth',2);
    loglog(atol,EC,'p','Color',clist(7,:),'MarkerSize',12,'LineWidth',2);
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Abs. residual tolerance [1]','Interpreter','latex','FontSize',15)
    ylabel('Rel. conservation error rate [1/s]','Interpreter','latex','FontSize',15)
    title('Global conservation with nonlinear convergence','Interpreter','latex','FontSize',18)

    if atol == ATOL(1)
        loglog(ATOL,eps.*ones(size(ATOL)),'k:' ,'LineWidth',2);  % plot trend for comparison
    end
    if atol == ATOL(end)
        legend('error $\bar{\rho}$','error $M$','error $X$','error $F$','error $S$','error $C$','machine prec.','Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [outdir,'/',runID,'/',runID,'_',TINT,'_',ADVN];
print(fh24,name,'-dpng','-r300','-vector');