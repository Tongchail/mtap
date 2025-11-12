%    The final monte carlo  10% without AMP,BI  with qtz  15.10.2025
%% prepare workspace
clear all; close all;

addpath(genpath('../'));
addpath('../../cal')
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
Fs = {'FontSize',12};
FS = {'FontSize',15};
FL = {'FontSize',18};
Ms = {'MarkerSize',8};
MS = {'MarkerSize',10};
ML = {'MarkerSize',12};
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
LB = {'Location','best'};
LO = {'Location','bestoutside'};
%% *****  prepare for pseudo-component calibration  *********************** （ Crucial！！！）
% 参数拟合的算法来优化这个结果，以使一个端元能更好地代表分异终点，而另一个端元则更好地代表起始点。

% !!!  Run Section to load end-member calibration and prepare for pseudo-component calibration  !!!
load('MtAp_frac_noampbi10_750_0.mat');

cal_MtAp_750;  % read calibration file

% !!!  Edit end-member appearances in pseudo-components  !!!
% - number and sequence of end-members must correspond to list in cal file
% - in sequence of appearance, add one new mineral
%   system to next pseudo-component
% - in sequence of appearance, add one new end-member of each mineral
%   system to next pseudo-component
% - phase out mineral systems and their end-members in accordance with
%   their fading or disappearance in PHS_frc



%                  ant alb san| mmt	tmt	mgt| mhy fhy hyp | mau fau aug |	ilm	|  qtz |  wat
indmem  = logical([ 1	0	0	 0	 0	 0	  0	  0	  0	    0	0	0		 0	   0	  0
                    1   1   0	 1	 0	 0	  1	  0	  0	    0	0	0		 0	   0      0
                    1	1	1	 1	 1	 1	  1	  1	  0	    1	0	0		 0	   0	  0
                    1	1	1	 0	 1	 1	  0	  0	  1	    0	1	0		 1	   0 	  0
                    0	1	1	 0	 0	 0	  0	  0	  1	    0	0	1		 1	   1	  0
                    0	0	0	 0	 0	 0	  0	  0	  0	    0	0	0		 0	   0 	  1]);


 
%diagonally downwards对角线向下排列的方式，左上角是less evolved更原始的物质，右下角是高度分异后才出现的物质highly differentiated or fractionated。

cmp_oxd = 1.0*EMInt + 0.0*EMExt; %*cal.mem_oxd(1:end-1,1:end-1)/100;%EM_Int，之前通过PCA和端元提取算法得到的内部端元
cmp_oxd = [cmp_oxd,zeros(cal.ncmp-1,1)];
cmp_oxd = [cmp_oxd;zeros(1,cal.noxd)];
cmp_oxd(end,end) =100;

cmp_mem = zeros(ncmp-1,cal.nmem);
for ic = 1:ncmp
    imem  = indmem(ic,:);
    A     = cal.mem_oxd(imem,:).';
    b     = squeeze(cmp_oxd(ic,:));
    cmp_mem(ic,imem) = lsqregcmp(A,b,[0.5,5,0,0.002])*100;
end


%             ant          alb       san   |    mmt	     tmt	    mgt|      chy       fhy       hyp|      mau       fau      aug  |     ilm	    qtz      wat
 cmp_mem = [  94.0000    5.0000         0         1         0         0         0         0         0         0         0         0         0         0         0
              32.6000   16.6000         0    8.7000       0.5       5.8   25.5000         0         0    10.400         0         0         0         0         0
              29.0000   14.4000    0.4000       0.5    6.6000    1.0000    1.0000   17.0000         0    4.5000   25.6000         0         0         0         0
              10.2000   55.5000    6.4000         0         0    2.2000         0         0    6.4000         0         0   13.0000    6.3000         0         0
                    0   23.5000   43.8000         0         0         0         0         0    1.4000         0         0         0    1.3000   30.0000         0
                    0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000];





%              ant          alb       san   |    mmt	    tmt	      mgt|      chy       fhy       hyp|      mau       fau      aug  |     ilm	    qtz      wat
% cmp_mem = [ 90.0000    9.0000         0     1.0000         0         0         0         0         0         0         0         0         0         0         0
%             30.1000   15.3000         0    12.0000    0.1000    5.8000   28.1000         0         0    8.5000         0         0         0         0         0
%             29.5000   14.6000    0.1000          0    6.6000    0.7000    0.9000   16.9000         0    4.2000   26.5000         0         0         0         0
%             10.4000   53.8000    7.0000          0    0.1000    2.1000         0    0.1000    6.3000         0         0   13.1000    7.1000         0         0
%                   0   28.0000   39.3000          0         0    0.1000         0         0    1.9000         0         0         0    1.2000   29.4000         0
%                   0         0         0          0         0         0         0         0         0         0         0         0         0         0  100.0000];



% This is now best
%            ant      alb       san      |   mmt	    tmt	     mgt|      chy      fhy       hyp|       mau       fau        aug  |    mil	     ilm	    fil|      wat
% cmp_mem =[ 86.0000   14.0000         0         0         0         0         0         0         0         0         0         0         0         0         0         0
%            74.1000   18.2000    0.2000    5.5000         0         0    3.4000         0         0         4         0         0       0.5         0         0         0
%             3.0000    6.2000    0.1000   11.6000    2.8000   10.4000   34.9000    1.9000         0   29.0000         0         0         0         0         0         0
%            10.0000   41.9000    0.2000    0.1000    4.6000    0.9000         0    1.1000   13.9000    1.6000   22.3000    2.1000         0         0       0.4         0
%                  0   20.2000   49.0000         0         0    3.3000         0    6.9000    0.1000         0    0.2000   13.4000         0    3.4000    0.5000         0
%                  0    4.8000   77.8000         0         0    0.8000         0         0   12.3000         0         0    4.0000         0    0.9000    0.1000         0
%                  0         0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000];


PHS_frc(:,end+1) = 0;

cmp_mem_init = round(cmp_mem,1);
cmp_mem_init = cmp_mem_init./sum(cmp_mem_init,2)*100;
indmem = logical(cmp_mem_init);
cmp_mem_best = cmp_mem_init;


cmp_oxd_init = cmp_mem_init*cal.mem_oxd/100;
cmp_oxd_best = cmp_oxd_init;


% set initial guess for melting point parameters 在拟合过程中用到的组分-端元初始条件表格
T0_init = [ 1650   1105   1057  1005   880];  T0_best = T0_init; % initial: T0_init = [ 1650   1150   1100  1000   850]; 
A_init  = [ 1e16  1e16   1e16   1e16   1e16];  A_best =  A_init;  %1e16
B_init  = [ 1  1   1   1   1];   B_best =  B_init;
r_init  = [ 35  3.4  3.2  5.8  5.5];   r_best =  r_init;%r_init  = [40.00   3.50   3.50   8.90  3.00   3.00];  
dT_init = 1400 * 1200./T0_init;  dT_best = dT_init;


% compose initial parameter guess
m0     = [T0_init.';A_init.';B_init.';r_init.';dT_init.';cmp_mem_init(:).*indmem(:);];

% set function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitNP(model,Tmp,Prs,SYS_oxdp,PHS_frc,cal,[0.5,5,1.75,1e-3]);

% test fit function for initial guess
[~,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit] = dhatFunc(m0);

% evaluate melting points as function of pressure
% if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
% %PP         = linspace(0.001,max(Psl)*3,50).';
% PP         = linspace(0.001,1.25*3,50).';
% var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
% cal.T0     = T0_init;
% cal.A      = A_init;
% cal.B      = B_init;
% cal.r      = r_init; 
% cal.dTH2O  = dT_init;
% var.c      = repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(PP),1);   % component fractions [wt]
% var.P      = PP/10;
% var.T      = 1000+PP*5e-8;
% var.H2O    = zeros(size(PP));
% cal.H2Osat = var.H2O+0.01;
% [~,cal]    = leappart(var,cal,'T');
% Tm         = cal.Tm;

% rescale mineral phase fractions
PHS_frc(:,2:end) = PHS_frc(:,2:end)./(100-PHS_frc(:,1)+eps)*100;

% plot basic information for initial fit
level = 3;
run('../MCMC/PlotFitFrac.m')
% 拟合算法的目的之一就是确定: which mineral end members are represented at which fration
% in each of the pseudo components. 在每个伪组分中，各矿物端元以什么比例出现。

%% *****  calibrate pseudo-components and melting point parameters  *******   蒙特卡洛算法（brun-in & annealing）

% cal_MORB_hi;  % read calibration file
Psl = 1.25;

% uncomment following lines to run MCMC again with previous best fit as initial guess
T0_init = T0_best; A_init = A_best; B_init = B_best; r_init = r_best; dT_init = dT_best; cmp_mem_init = cmp_mem_best;
m0      = [T0_init.';A_init.';B_init.';r_init.';dT_init.';cmp_mem_init(:).*indmem(:)];

% *** !!!  set MCMC parameters then Run Section dbstack   to execute MCMC routine  !!!
Niter           = 1e4;               % 1e6 number of samples to take; how many samples to take during one Monte Carlo run (how many times you do a random shuffle)
anneal.initstep = 1e-3;              % adjust step size to get reasonable acceptance ratio 20-30% (the step size, so like how far out from your current solution guess do you randomly shuffle to do the next test你在当前解的基础上，随机扰动多远来进行下一次测试。)
anneal.levels   = 1;                 % select number of annealing levels
anneal.burnin   = max(1,Niter/1e3);  % set length of initial burn-in sequence 定义了烧入期长度
anneal.refine   = max(1,Niter/1e3);  % set length of final refinement sequence

% !!!  set data uncertainties to weight likelihood function  !!!
MLT_scl   = max(0.01,(MLT_oxdp(:)-min(MLT_oxdp(:)))./(max(MLT_oxdp(:))-min(MLT_oxdp(:))));
SOL_scl   = max(0.01,(SOL_oxdp(:)-min(SOL_oxdp(:)))./(max(SOL_oxdp(:))-min(SOL_oxdp(:))));
MEM_scl   = max(0.01,(SOL_mem (:)-min(SOL_mem (:)))./(max(SOL_mem (:))-min(SOL_mem (:))));
PHS_scl   = max(0.01,(PHS_frc (:)-min(PHS_frc (:)))./(max(PHS_frc (:))-min(PHS_frc (:))));
sigma_MLT =  0.1  * MLT_scl.^0.25;       % uncertainty of melt oxide composition
sigma_SOL =  0.1  * SOL_scl.^0.25;       % uncertainty of melt oxide composition
sigma_MEM =  0.1  * MEM_scl.^0.25;       % uncertainty of solid end-member composition
sigma_PHS =  0.1  * PHS_scl.^0.25;       % uncertainty of phase fractions
% （****）sigma_TSL =  0.1  * ones(size([Tsol(:);Tliq(:)])); % uncertainty of solidus/liquidus Temp
%sigma = [sigma_MLT;sigma_SOL;sigma_MEM;sigma_PHS;sigma_TSL]; % combine all as in data vector
sigma = [sigma_MLT;sigma_SOL;sigma_MEM;sigma_PHS]; % combine all as in data vector

% load calibration data constraints into data vector
%data   = [MLT_oxdp(:);SOL_oxdp(:);SOL_mem(:);PHS_frc(:);Tsol(:);Tliq(:)];
data   = [MLT_oxdp(:);SOL_oxdp(:);SOL_mem(:);PHS_frc(:)];

% construct lower and upper parameter bounds
dm    =[1*max(20.0,0.02*T0_init.'); ...
        0*max(0.25,0.25* A_init.'); ...   
        0*max(0.25,0.25* B_init.'); ...
        1*max(0.50,0.40* r_init.'); ...%0.4
        0*max(10.0,0.25*dT_init.'); ...
        1*max(5.00,min(20.0,cmp_mem_init(:))).*indmem(:)];
m0_lw  = m0 - dm;
m0_up  = m0 + dm;
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(1*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(0.5,                  mbnds(1*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(2*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(1.0,                  mbnds(2*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(3.0,                  mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(100,                  mbnds(4*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/10,min(99.9,mbnds(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:))); %约束等于 5 × 组分数，再减去 1（用来处理 12345 的情况），再乘上 除水以外每个组分的一个参数。组分数 × 端元数（end members），这就是表格里的边（edges）。如果初始值（上面初始化的）有某个落在分布的边界之外（=0），那么 MCMC 就会正好以这种方式崩溃或者卡住。


mbnds(m0==100 ,:) = 100;
% mbnds(m0==7.0 ,1) = 5.0;


mbnds(m0==T0_init(1),:) = T0_init(1);
mbnds(m0==T0_init(end),:) = T0_init(end);
anneal.initstep = anneal.initstep * diff(mbnds,1,2);  % resize step according to bounded bracket

% set parameter names according to info from calibration file
mNames = cell(cal.ncmp*cal.nmem+4*(cal.ncmp-1),1);
k = 1;
for j=1:5
    for i=1:cal.ncmp-1
        if j==1
            mNames{k} = ['T0:',cal.cmpStr{i}];
        elseif j==2
            mNames{k} = ['A:',cal.cmpStr{i}];
        elseif j==3
            mNames{k} = ['B:',cal.cmpStr{i}];
        elseif j==4
            mNames{k} = ['r:',cal.cmpStr{i}];
        elseif j==5
            mNames{k} = ['dT_H:',cal.cmpStr{i}];
        end
        k = k+1;
    end
end
for j=1:cal.nmem
    for i=1:cal.ncmp
        mNames{k} = [cal.memStr{j},':',cal.cmpStr{i}];
        k = k+1;
    end
end

% set function to apply further constraints to a proposed set of parameter values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% set function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% set function to calculate likelihood of forward model
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.1,0.1,max(Psl)*3,model,cal);

bestfit = m0;  % initialise bestfit from initial conditions


%*****  RUN MCMC PARAMETER FITTING ROUTINE  *******************************
% 	•	这行把烧入期/退火信息一起传给了 MCMC 例程：mcmc 会读取 anneal.burnin，在开始的前 Niter/5 次迭代按“烧入期”处理（通常接受率更宽松/步长更大，具体由 mcmc 内部实现）。
tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),m0,mbnds,anneal,Niter);
RunTime(1) = toc;
%**************************************************************************


% uncomment following line to plot likelihood histograms for fitted parameters (slow!)
% plotmcmc(models, prob, [], mbnds, anneal, mNames); 

% extract best fit parameters
T0_best       = bestfit(               (1:cal.ncmp-1)).';
A_best        = bestfit(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
B_best        = A_best; %bestfit(2*(cal.ncmp-1)+(1:cal.ncmp-1)).'; %A_best;
r_best        = bestfit(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';  % a set of fitting parameters which have to do with the curvature of the phase loops that this model is generating.
dT_best       = 1400 * 1200./T0_best;% bestfit(4*(cal.ncmp-1)+(1:cal.ncmp-1)).';  这个 dt 参数与固相线和液相线的含水依赖性有关。它表示在系统中加入一定量的水后，熔点会降低多少。高熔点的伪组分在熔融行为上受水的影响比低熔点的要小。
cmp_mem_best  = reshape(bestfit(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd_best  = cmp_mem_best*cal.mem_oxd/100;

% evaluate forward model for best fit parameters
%[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,Tsolfit,Tliqfit,~] = dhatFunc(bestfit);
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit] = dhatFunc(bestfit);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

% evaluate melting points as function of pressure
% if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
% PP         = linspace(0.001,max(Psl)*2,50).';  
% var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
% cal.T0     = T0_best;
% cal.A      = A_best;
% cal.B      = B_best;
% cal.r      = r_best;
% cal.dTH2O  = dT_best;
% var.c      =repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(PP),1);   % component fractions [wt] 
% var.P      = PP/10;
% var.T      = 1000+PP*5e-8;
% var.H2O    = zeros(size(PP));
% cal.H2Osat = var.H2O+0.001;
% [~,cal]    = leappart(var,cal,'T');
% Tm         = cal.Tm;

% uncomment following line to retrieve probability distributions (slow!)
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% *****  visualise best fit calibration  *********************************

%!!!  adjust desired level of detail to plot then run section  !!!
%     level = 1     only simple line plots
%     level = 2     add Harker diagrams and T-X pseudo-sections
%     level = 3     add mineral systems and display best fit parameters

level = 3;
run('../MCMC/PlotFitFrac.m')

%% save and display calibration
save('MtAp_noamp10_calibration_2');
