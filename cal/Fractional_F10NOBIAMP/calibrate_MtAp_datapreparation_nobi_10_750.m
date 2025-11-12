% Fractional crystalization - model  (without amp and bi  residual fraction = 10%)    06 Oct. 2025
% Temperature range 1100 to 750℃ 

%**************************************************************************
%*****  CALIBRATE MULTI-COMPONENT MELTING MODEL  **************************
%**************************************************************************

% prepare workspace
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


%% *****  load calibration data  ******************************************

% load MAGEMin results (window with table opens, using default selection
% click Import Selection => Import Data

filename = './fractional_10_750_ig.csv';  % residual fraction = 10%
uiopen(filename,1)
% *** we do not select liq=100 or liq < 10 wt% 

%% *****  unpack calibration data  ****************************************

% !!!  update table name on following line, then Run Section  !!!
DAT = fractional10750ig;       % residual fraction = 10%        % table name must correspond to table header above

% load phase names in order of appearance, liq first
phs = unique(string(DAT.phase),'stable');                                  % load phase list
phs(phs=='system') = [];                                                   % discard system
phs(phs=='qfm') = [];                                                      % discard fO2 buffer
phs(phs=='fl') = [];                                                       % discard fluid phase
nphs = length(phs);                                                        % record number of phases
iliq = find(strcmp(phs,'liq'));                                            % ensure liq comes first
iphs = 1:nphs; iphs(iliq) = []; iphs = [iliq,iphs];
phs  = phs(iphs);

liq = 1;  pl = 2; mgt = 3; opx = 4; cpx = 5; ilm = 6;    % set shortcut phase indices
% afs-alkali feldspar   ru-rutile (TiO2)


% set oxide list in preferred sequence   给氧化物重新排序
oxd  = ["SiO2";"TiO2";"Al2O3";"FeO";"MgO";"CaO";"Na2O";"K2O";"H2O"];       % set major oxides
noxd = length(oxd);                                                        % record number of oxides

% extract calculation points
pts  = unique(DAT.point,'stable'); offset = min(pts)-1; pts = pts-offset;  % point numbers
Tmp  = unique(DAT.TC,'stable');                                            % point temperatures
Prs  = unique(DAT.Pkbar,'stable');                                         % point pressures
npts = length(pts);                                                        % number of points

Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; K = 8; H = 9;      % set shortcut oxide indices

% detect which phases are stable on which points
% 用 0/1 标记每个点是否存在某个矿物相。因为有的时候并不是所有的相都在每一点上存在，比如说一开始结晶的是斜长石还有尖晶石，然后是辉石等
hasphs = zeros(npts,nphs);
for iph = 1:nphs
    for ipt = 1:npts
        hasphs(ipt,iph) = any(table2array(DAT(DAT.point==ipt+offset,'phase'))==phs(iph));
    end
end

% extract phase fractions in [wt%]
PHS_frc = zeros(npts,nphs);
for iph = 1:nphs
    PHS_frc(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'modewt'));
end
PHS_frc = PHS_frc./(sum(PHS_frc,2)+eps)*100;

% extract phase densities in [kg/m3]
RHO = zeros(npts,nphs);
for iph = 1:nphs
    RHO(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'densitykgm3'));
end

% extract phase oxide compositions in [wt%]
PHS_oxd  = zeros(npts,nphs,noxd);
for iph = 1:nphs
    PHS_oxd(hasphs(:,iph)==1,iph,:) = table2array(DAT(DAT.phase==phs(iph),{'SiO2wt','TiO2wt','Al2O3wt','FeOwt','MgOwt','CaOwt','Na2Owt','K2Owt','H2Owt'}));
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
end

% % lump in rutile with ilmenite
% PHS_oxd(:,ilm,:) = (PHS_frc(:,ilm).*PHS_oxd(:,ilm,:) + PHS_frc(:,ru).*PHS_oxd(:,ru,:)) ./ (PHS_frc(:,ilm) + PHS_frc(:,ru) + eps);
% PHS_oxd(:,ru,:) = [];
% RHO(:,ilm)       = (PHS_frc(:,ilm)+PHS_frc(:,ru))./(PHS_frc(:,ilm)./(RHO(:,ilm)+eps) + PHS_frc(:,ru)./(RHO(:,ru)+eps) + eps);
% RHO(:,ru)       = [];
% PHS_frc(:,ilm)   =  PHS_frc(:,ilm) + PHS_frc(:,ru); 
% PHS_frc(:,ru)   = [];
% phs(ru)         = [];
% hasphs(:,ilm)    = max(hasphs(:,ilm),hasphs(:,ru));
% hasphs(:,ru)    = [];
% nphs             = nphs-1;
% % lump in afs with pl
% PHS_oxd(:,pl,:) = (PHS_frc(:,pl).*PHS_oxd(:,pl,:) + PHS_frc(:,afs).*PHS_oxd(:,afs,:)) ./ (PHS_frc(:,pl) + PHS_frc(:,afs) + eps);
% PHS_oxd(:,afs,:) = [];
% RHO(:,pl)       = (PHS_frc(:,pl)+PHS_frc(:,afs))./(PHS_frc(:,pl)./(RHO(:,pl)+eps) + PHS_frc(:,afs)./(RHO(:,afs)+eps) + eps);
% RHO(:,afs)       = [];
% PHS_frc(:,pl)   =  PHS_frc(:,pl) + PHS_frc(:,afs); 
% PHS_frc(:,afs)   = [];
% phs(afs)         = [];
% hasphs(:,pl)    = max(hasphs(:,pl),hasphs(:,afs));
% hasphs(:,afs)    = [];
% nphs             = nphs-1;



liq = 1;  pl = 2; mgt = 3; opx = 4; cpx = 5; ilm = 6;            % update phase indices
% *** here pl = fsp = pl + afs

% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

% remove minor oxides from phases (mean<0.60; max<1.0)
for iph = 1:nphs
    ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.60 & max(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),[],1)<1.00);
    hasoxd(iph,ilim) = false;
    PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
end
PHS_oxdp = PHS_oxd;

% extract melt oxide composition
MLT_oxd = squeeze(PHS_oxd(:,1,:));  % melt oxide composition

% extract solid phase oxide composition
SOL_oxd = zeros(npts,noxd);
wt  = zeros(size(SOL_oxd)) + eps;
for iph = 2:nphs
    SOL_oxd = SOL_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SOL_oxd = SOL_oxd./wt;  % solid oxide composition

% extract system oxide composition
SYS_oxd = zeros(npts,noxd);
wt  = zeros(size(SYS_oxd)) + eps;
for iph = 1:nphs
    SYS_oxd = SYS_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxd = SYS_oxd./wt;  % system oxide composition

%%  Rename and Save.mat

save('MtAp_data_10_750_0.mat', 'PHS_frc', 'PHS_oxd', 'PHS_oxdp', 'MLT_oxd','SOL_oxd','SYS_oxd','RHO','phs','hasphs','pts','Tmp','Prs','npts','nphs');


