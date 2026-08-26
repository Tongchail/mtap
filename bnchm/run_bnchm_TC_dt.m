% prepare workspace
clear; close all;

% load default parameters
run('../usr/par_MtAp_default')

% set run parameters
runID    =  'bnchm_TC_dt';       % run identifier
opdir    =  '../out/';           % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  10;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot of results
plot_cv  =  0;                   % switch on to live plot iterative convergence
save_op  =  0;

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  nop;                 % number of time steps to take
dt       =  1;                   % set initial time step

% set initial thermo-chemical state
smth     =  15;
init_mode=  'liquidus';          % init_mode = 'constant', 'liquidus', 'layer','linear', 'chamber'.
T0       =  -200;                 % initial temperature [deg C] 
c0       =  [16 11 16 19 38 10 2]/100;  % *** components (maj comp, H2O) top  layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,1,-1,-1,-1,0]*0e-3;  % amplitude of random noise [wt]
dcg      =  [-1,-1,-1,1,1,1,0]*5e-2;  % amplitude of centred gaussian [wt]
dTg      =  10;
dTr      =  0.0;
dr_trc   =  [1,1,1,-1,-1,-1].*0e-3;
dg_trc   =  [-1,-1,-1,1,1,1].*1e-2;

% closed boundaries for gas flux
periodic =  1;
bndmode  =  0;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
Ptop     =  2.0e8;               % top pressure [Pa]
fin      =  0;
fout     =  0;
tau_r    =  1e32;
calID    =  'MtAp_750_new';              % phase diagram calibration

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1;                   % (physical) time stepping courant number (multiplies stable step) [0,1]
atol     =  1e-9;               % outer its absolute tolerance
rtol     =  atol/1e6;            % outer its absolute tolerance
maxit    =  100;                 % maximum outer its
itpar.fp.damp = 1;                % fixed-point iterative damping (0-1)
itpar.aa.m    = 2;                % Anderson acceleration depth (2-5)
itpar.aa.damp = 0.0;              % Anderson acceleration damping (0-1)
itpar.aa.reg  = 0.01;             % Anderson acceleration regularisation (0-1)

% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

cd ../src

DT    = [h/2,h/4,h/8];  % time step sizes to test relative to grid step
nshft = 2;              % number of grid steps target is shifted from initial

for dti = DT
    
    dt    =  dti;
    dtmax =  dti;
    Nt    =  nshft*h/dti;
    L0    =  mean(L0(:));

    Delta_cnv = h;

    % initialise fields
    init;

    % set velocities to constant values for lateral translation with no segregation
    W(:) = 0;  Wm(:) = 0;  Wx(:) = 0;  Wf(:) = 0;  wx(:) = 0;  wm(:) = 0;  wf(:) = 0;  upd_W(:) = 0;
    U(:) = 0;  Um(:) = 1;  Ux(:) = 1;  Uf(:) = 1; upd_U(:) = 0;   
    P(:) = 0;  upd_P(:) = 0; upd_MFS(:) = 0;

    % set diffusion parameters to zero to isolate advection
    kT(:)   = 0;  ks(:) = 0;  kx(:) = 0;  ks_x(:) = 0;  ks_f(:) = 0;  ke(:) = 0;
    diss(:) = 0;
    res_rho = 0.*rho;

    % set parameters for non-dissipative, non-reactive flow
    rhoin = rho; rhoout = circshift(rho,nshft,2);
    Sin   = S;   Sout   = circshift(S  ,nshft,2);
    Fin   = F;   Fout   = circshift(F  ,nshft,2);
    Min   = M;   Mout   = circshift(M  ,nshft,2);
    Xin   = X;   Xout   = circshift(X  ,nshft,2);

    dt    = dti;
    dtmax = dti;
    time  = 0;

    output;

    % physical time stepping loop
    while time <= tend && step <= Nt

        % time step info
        timing;

        % store previous solution
        store;

        % reset residuals and iteration count
        resnorm  = 1;
        resnorm0 = resnorm;
        iter     = 1;

        % non-linear iteration loop
        while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit

            % solve thermo-chemical equations
            thermochem;

            % update non-linear parameters and auxiliary variables
            update;

            kT(:) = 0;  ks(:) = 0;  kx(:) = 0;  ks_x(:) = 0;  ks_f(:) = 0;  ke(:) = 0;  diss(:) = 0;

            % report convergence
            report;

            iter = iter+1;
        end

        % print model diagnostics
        diagnose;

        % plot model results
        if ~mod(step,nop); output; end

        % increment time/step
        time = time+dt;
        step = step+1;
        if frst; frst=0; end

        figure(100); clf;
        plot(XX(ceil(Nz/2),:), Xout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'k',XX(ceil(Nz/4),:), X(ceil(Nz/2),:)./rho(ceil(Nz/2),:),'r','LineWidth',1.5); axis tight; box on;
        set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
        xlabel('Distance [m]','Interpreter','latex','FontSize',16)
        ylabel('Crystallinity [wt]','Interpreter','latex','FontSize',16)
        drawnow;

    end

    figure(200);
    plot(XX(ceil(Nz/2),:),X(ceil(Nz/2),:)./rho(ceil(Nz/2),:)-Xout(ceil(Nz/2),:)./rhoout(ceil(Nz/2),:),'LineWidth',1.5); axis tight; box on; hold on;
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Distance [m]','Interpreter','latex','FontSize',16)
    ylabel('Residual [wt]','Interpreter','latex','FontSize',16)

    % plot convergence
    % EB = norm(rho-rhoout)./norm(rhoout);
    EM = norm(  M-  Mout)./norm(Mout);
    EX = norm(  X-  Xout)./norm(Xout);
    EF = norm(  F-  Fout)./norm(Fout);
    ES = norm(  S-  Sout)./norm(Sout);

    clist = [colororder;[0 0 0]];

    fh24 = figure(24);
    loglog(dt,EM,'o','Color',clist(2,:),'MarkerSize',10,'LineWidth',2); hold on; box on;
    loglog(dt,EX,'d','Color',clist(1,:),'MarkerSize',10,'LineWidth',2);
    loglog(dt,EF,'^','Color',clist(3,:),'MarkerSize',10,'LineWidth',2);
    loglog(dt,ES,'v','Color',clist(4,:),'MarkerSize',10,'LineWidth',2);
    set(gca,'LineWidth',1.5,'TickLabelInterpreter','latex','FontSize',12)
    xlabel('Time step [s]','Interpreter','latex','FontSize',15)
    ylabel('Rel. numerical error [1]','Interpreter','latex','FontSize',15)
    title('Numerical convergence in time','Interpreter','latex','FontSize',18)

    if round(dt,4) == DT(1)
        % loglog(DT,geomean([EM,EX,EF,ES]).*(DT./DT(1)).^1,'k--','LineWidth',2);  % plot trend for comparison
        loglog(DT,geomean([EM,EX,EF,ES]).*(DT./DT(1)).^2,'k-' ,'LineWidth',2);
    end
    if round(dt,4) == DT(end)
        legend({'error $M$','error $X$','error $F$','error $S$','linear','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end
    drawnow;

end

name = [opdir,'/',runID,'/',runID,'_',TINT];
print(fh24,name,'-dpng','-r300','-vector');