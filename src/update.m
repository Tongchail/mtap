%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update phase indicators
hasx = x >= eps^0.5;
hasf = f >= eps^0.5;
hasm = m >= eps^0.5;

% update phase oxide compositions
c_oxd  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cf_oxd = reshape(reshape(cf,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% update phase mineral end-member compositions
c_mem  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cm_mem = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cx_mem = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cf_mem = reshape(reshape(cf,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);

% update phase mineral systems composition for solid assemblage
c_msy  = reshape(reshape( c_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cm_msy = reshape(reshape(cm_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cx_msy = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cf_msy = reshape(reshape(cf_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
cx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    cx_msy_oxd(:,:,j,:) = reshape(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,nnz(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,nnz(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
end

cf_oxd_all = zeros(size(c,1),size(c,2),9);
cm_oxd_all = zeros(size(c,1),size(c,2),9);
cx_oxd_all = zeros(size(c,1),size(c,2),9);
 c_oxd_all = zeros(size(c,1),size(c,2),9);
if cal.noxd>9
    cf_oxd_all = cf_oxd(:,:,cal.ioxd);
    cm_oxd_all = cm_oxd(:,:,cal.ioxd);
    cx_oxd_all = cx_oxd(:,:,cal.ioxd);
     c_oxd_all =  c_oxd(:,:,cal.ioxd);
else
    cf_oxd_all(:,:,cal.ioxd) = cf_oxd;
    cm_oxd_all(:,:,cal.ioxd) = cm_oxd;
    cx_oxd_all(:,:,cal.ioxd) = cx_oxd;
     c_oxd_all(:,:,cal.ioxd) = c_oxd;
end

% get trace element phase compositions
Ktrc = zeros(Nz,Nx,cal.ntrc);
trcm = zeros(Nz,Nx,cal.ntrc);
trcx = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* c_mem(:,:,j)./100; end

    trcm(:,:,i)  = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcx(:,:,i)  = trc(:,:,i)./(m./Ktrc(:,:,i) + x);
end

% update phase densities
rhom0  = reshape(DensityX(reshape(cm_oxd_all,Nz*Nx,9),Tref,Pref./1e8)    ,Nz,Nx);
rhox0  = reshape(sum(reshape(cx_mem/100,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx);
rhof0  = reshape(DensityX(reshape(cf_oxd_all,Nz*Nx,9),Tref,Pref./1e8)    ,Nz,Nx);

rhom   = rhom0 .* (1 - aTm.*(T-Tref) + bPm.*(Pt-Pref));
rhox   = rhox0 .* (1 - aTx.*(T-Tref) + bPx.*(Pt-Pref));
rhof   = rhof0 .* (1 - aTf.*(T-Tref) + bPf.*(Pt-Pref));

rho0   = 1./(m./rhom0 + x./rhox0 + f./rhof0);
rho    = 1./(m./rhom  + x./rhox  + f./rhof );

rhomw  = (rhom(icz(1:end-1),:)+rhom(icz(2:end),:))/2;
rhoxw  = (rhox(icz(1:end-1),:)+rhox(icz(2:end),:))/2;
rhofw  = (rhof(icz(1:end-1),:)+rhof(icz(2:end),:))/2;

rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;

rhoref = mean(rhow,2);

Drhom  = rhomw - rhow;
Drhox  = rhoxw - rhow;
Drhof  = rhofw - rhow;
Drho   = rhow  - rhoref;

rhoW   = rhow.*W(:,2:end-1);
rhoU   = rhou.*U(2:end-1,:);

% convert weight to volume fraction, update bulk density
chi    = max(eps,min(1-eps, x.*rho./rhox));
mu     = max(eps,min(1-eps, m.*rho./rhom));
phi    = max(eps,min(1-eps, f.*rho./rhom));
% chi    = max(0,min(1, x.*rho./rhox ));
% phi    = max(0,min(1, f.*rho./rhof ));
% mu     = max(0,min(1, m.*rho./rhom ));

chiw   = (chi(icz(1:end-1),:)+chi(icz(2:end),:))./2;
muw    = ( mu(icz(1:end-1),:)+ mu(icz(2:end),:))./2;
phiw   = (phi(icz(1:end-1),:)+phi(icz(2:end),:))./2;

xw     = (x(icz(1:end-1),:)+x(icz(2:end),:))./2;
mw     = (m(icz(1:end-1),:)+m(icz(2:end),:))./2;
ffw    = (f(icz(1:end-1),:)+f(icz(2:end),:))./2; % because the name 'fw' is already used, we use 'ffw' instead.

xu     = (x(:,icx(1:end-1))+x(:,icx(2:end)))./2;
muu    = (m(:,icx(1:end-1))+m(:,icx(2:end)))./2;
fu     = (f(:,icx(1:end-1))+f(:,icx(2:end)))./2;

Xw     = (X(icz(1:end-1),:)+X(icz(2:end),:))/2;
Mw     = (M(icz(1:end-1),:)+M(icz(2:end),:))/2;
Fw     = (F(icz(1:end-1),:)+F(icz(2:end),:))/2;

chi_mem = reshape(reshape(cx_mem/100.*rhox,Nz*Nx,cal.nmem)./cal.rhox0,Nz,Nx,cal.nmem);
chi_mem = chi_mem./sum(chi_mem,3);

% update thermal parameters
aT = mu.*aTm + chi.*aTx + phi.*aTf;
bP = mu.*bPm + chi.*bPx + phi.*bPf;
kT = mu.*kTm + chi.*kTx + phi.*kTf;
cP = mu.*cPm + chi.*cPx + phi.*cPf;
RhoCp = mu.*rhom.*cPm + chi.*rhox.*cPx + phi.*rhof.*cPf;
Adbt  = mu.*aTm./rhom./cPm + chi.*aTx./rhox./cPx + phi.*aTf./rhof./cPf;

% % update lithostatic pressure
% Pti = Pt;
% if Nz==1; Pt    = max(Ptop/100,Ptop.*ones(size(Pt)) + Pcouple*(Pchmb + P(2:end-1,2:end-1))); else
%     Pl(1,:)     = repmat(rhoref(1).*g0.*h/2,1,Nx) + Ptop;
%     Pl(2:end,:) = Pl(1,:) + repmat(cumsum(rhoref(2:end-1).*g0.*h),1,Nx);
%     Pt          = max(Ptop/100,Pl + Pcouple*(Pchmb + P(2:end-1,2:end-1)));
% end
% Pt = (Pt + Pti)/2;
% update lithostatic pressure
Pl(1,:)     = repmat(rhoref(1).*g0.*h/2,1,Nx) + Ptop;
Pl(2:end,:) = Pl(1,:) + repmat(cumsum(rhoref(2:end-1).*g0.*h),1,Nx);
Pt          = max(Ptop/100,Pl + P(2:end-1,2:end-1));

% update effective constituent sizes
dm = dm0.*(1-mu ).^0.5;
dx = dx0.*(1-chi).^0.5;
df = df0.*(1-phi).^0.5;

% update pure phase viscosities
etam   = reshape(Giordano08(reshape(cm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);
etax0  = reshape(prod(cal.etax0(1:end-1).^reshape(chi_mem(:,:,1:end-1)+eps,Nz*Nx,cal.nmem-1),2),Nz,Nx);
etax   = etax0 .* ones(size(chi)) .* exp(cal.Eax./(8.3145.*T)-cal.Eax./(8.3145.*(Tref+273.15)));
etaf   = reshape(etamfe(reshape(cf_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);

% get coefficient contrasts
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);

% get permission weights
dd = max(eps^0.5,min(1-eps^0.5,permute(cat(3,dx ,dm ,df ),[3,1,2])));
ff = max(eps^0.5,min(1-eps^0.5,permute(cat(3,chi,mu ,phi),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./cal.BB).^(1./cal.CC);  Sf = Sf./sum(Sf,2);
Xf = sum(cal.AA.*Sf,2).*FF + (1-sum(cal.AA.*Sf,2)).*Sf;

% get momentum flux and transfer coefficients
thtv = squeeze(prod(Mv.^Xf,2));
etai = kv.*thtv;                % the phase viscosities
% Kv   = ff.*kv.*thtv; 
% Cv   = Kv.*(1-ff)./dd.^2;

% get effective viscosity
etax   = squeeze(etai(1,:,:));
etam   = squeeze(etai(2,:,:));
etaf   = squeeze(etai(3,:,:));
etamix = squeeze(sum(ff.*etai,1)); if Nx==1; etamix = etamix.'; end

% update velocity divergence
Div_V    = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                      % get velocity divergence
Div_rhoV = ddz(rhoW        ,h) + ddx(rhoU        ,h);                      % get mass flux divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V/3;                                 % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V/3;                                 % z-normal strain rate
exz = (diff(U,1,1)./h+diff(W,1,2)./h)/2;                                   % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% extract potential densities
rhomp = rhom0 .* (1 - aTm.*(Tp-Tref));
rhoxp = rhox0 .* (1 - aTx.*(Tp-Tref));
rhofp = rhof0 .* (1 - aTf.*(Tp-Tref));

rhop  = 1./(m./rhomp + x./rhoxp + f./rhofp);

% update velocity magnitude
if Nx==1 && Nz==1; V = 0;
elseif Nx==1
    idz = (1:Nz)';  % grid indices
    half_steps = max(1,floor(L0 ./ (2 * h)));  % half mixing length in grid steps

    ip = idz + half_steps;  % upper indices
    im = idz - half_steps;  % lower indices

    ip = min(ip, Nz);  % clamp indices to valid range [1, Nz]
    im = max(im, 1 );  % clamp indices to valid range [1, Nz]

    drhoz   = max(0, -(rhop(ip,:)-rhop(im,:)) ) + 1e-9.*rhop; % central density contrast across mixing length
    for i=1:10
        drhoz = drhoz + diffus(drhoz,1/8*ones(size(drhoz)),1,[1,2],BCD);
    end
    V = 2/9*drhoz.*g0.*(L0/2).^2./eta;
else
    V = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
           + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);
end

Vx   = sqrt(((Wx(1:end-1,2:end-1)+Wx(2:end,2:end-1))/2).^2 ...
          + ((Ux(2:end-1,1:end-1)+Ux(2:end-1,2:end))/2).^2);                % crystal convection speed magnitude
Vf   = sqrt(((Wf(1:end-1,2:end-1)+Wf(2:end,2:end-1))/2).^2 ...
          + ((Uf(2:end-1,1:end-1)+Uf(2:end-1,2:end))/2).^2);                % fluid convection speed magnitude

bndtaperx = (1 - (exp((-ZZ)/l0x) + exp(-(D-ZZ)/l0x)).*(1-open_sgr));
bndtaperf = (1 - (exp((-ZZ)/l0f) + exp(-(D-ZZ)/l0f)).*(1-open_sgr));

vx   = dx0^2./etas_x.*(rhox0-rhom0).*g0.*bndtaperx;                            % xtal segregation speed magnitude vx downward direction
vf   = df0^2./etas_f.*(rhof0-rhom0).*g0.*bndtaperf;  
vm   = ( x .* vx + f .* vf ) ./ max(m, eps);

xie   = sqrt(((xiew(1:end-1,2:end-1)+xiew(2:end,2:end-1))/2).^2 ...
           + ((xieu(2:end-1,1:end-1)+xieu(2:end-1,2:end))/2).^2); 
xisx  = sqrt(((xisxw(1:end-1,2:end-1)+xisxw(2:end,2:end-1))/2).^2 ...
           + ((xisxu(2:end-1,1:end-1)+xisxu(2:end-1,2:end))/2).^2);            % settling noise flux magnitude 
xiex  = sqrt(((xixw(1:end-1,2:end-1)+xixw(2:end,2:end-1))/2).^2 ...
           + ((xixu(2:end-1,1:end-1)+xixu(2:end-1,2:end))/2).^2);            % xtal eddy noise flux magnitude 模长速度噪声大小           % global eddy noise？melt? eddy noise flux magnitude
xix   = xisx + xiex;                                                          % total xtl noise
xisf = sqrt(((xisfw(1:end-1,2:end-1)+xisfw(2:end,2:end-1))/2).^2 ...
          + ((xisfu(2:end-1,1:end-1)+xisfu(2:end-1,2:end))/2).^2);          % settling noise flux magnitude 
xief = sqrt(((xifw(1:end-1,2:end-1)+xixw(2:end,2:end-1))/2).^2 ...
          + ((xifu(2:end-1,1:end-1)+xixu(2:end-1,2:end))/2).^2);            % fluid eddy noise flux magnitude
xif  = xisf + xief;                                                         % total fluid noise

% update diffusion parameters
if Nx==1 && Nz==1; ke = 0; fReL = 1; fRel_x = 1; fRel_f = 1;
elseif Nx==1
    ke     = (ke + V.*L0)/2;                                      % convective mixing diffusivity
    fReL   = 1;
    fRel_x = 1;
    fRel_f = 1;
else
    bndtapere = (1 - (exp((-ZZ)/L0) + exp(-(D-ZZ)/L0).*(1-open_cnv)));
    ke        = eII.*L0.^2 .* bndtapere;                                   % turbulent eddy diffusivity
    fReL    =  1-exp(-ReL  );         % Re-dependent ramp factor
    fRel_x  =  1-exp(-Rel_x);         % Re-dependent ramp factor
    fRel_f  =  1-exp(-Rel_f);         % Re-dependent ramp factor
end

ks_x = vx .*l0x;                                                       % xtl segregation diffusivity
ks_f = vf .*l0f;                                                       % flu segregation diffusivity
k_x  = (ks_x + fReL.*ke);                                             % regularised particle diffusivity
k_f  = (ks_f + fReL.*ke);                                             % regularised particle diffusivity
kh   = (ke./Prt.*fReL + kmin).*rho.*cP./T;                                  % regularised heat diffusion
kc   =  ke./Sct.*fReL + kmin;                                               % turbulent eddy diffusivity

% update viscosities
etae = fReL.*ke.*rho;                                                       % eddy viscosity
eta  = (eta + etamix + etae)/2;                                            % effective viscosity

etat_x = fRel_x.*ks_x.*rho;                                                  % turbulent drag viscosity
etat_f = fRel_f.*ks_f.*rho;                                                  % turbulent drag viscosity
etas_x = (etas_x + etax + etat_x)/2;                                       % crystal effective drag viscosity   
etas_f = (etas_f + etaf + etat_f)/2;                                       % fluid effective drag viscosity   

% limit total viscosity contrast
etamax = geomean(eta(:)).*(etacntr/2);
etamin = geomean(eta(:))./(etacntr/2);
eta    = 1./(1./etamax + 1./eta) + etamin;                                 % total magma viscosity

etamax = geomean(etas_x(:)).*(etacntr/2);
etamin = geomean(etas_x(:))./(etacntr/2);
etas_x   = 1./(1./etamax + 1./etas_x) + etamin;                            % effective crystal settling viscosity

etamax = geomean(etas_f(:)).*(etacntr/2);
etamin = geomean(etas_f(:))./(etacntr/2);
etas_f   = 1./(1./etamax + 1./etas_f) + etamin;                            % effective fluid droplet settling viscosity

% interpolate to staggered nodes
etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

etas_xw = (etas_x(icz(1:end-1),:).*etas_x(icz(2:end),:)).^0.5;
etas_fw = (etas_f(icz(1:end-1),:).*etas_f(icz(2:end),:)).^0.5;

% update dimensionless numbers
ReL   = V .*L0 ./(etamix./rho);                                            % Reynolds number on scaled domain length
ReD   = V .*D0 ./(eta   ./rho);                                            % Reynolds number on scaled domain length
Rel_x = vx.*dx0./(etax./rho);                                              % particle Reynolds number-crystal
Rel_f = vf.*df0./(etaf./rho);                                              % particle Reynolds number-fluid
Red_x = vx.*dx0./(etas_x./rho);                                            % particle Reynolds number-crystal
Red_f = vf.*df0./(etas_f./rho);                                            % particle Reynolds number-fluid
Ra    = V .*D0./(k_x+k_f+ke+kT./rho./cP);                                  % Rayleigh number on scale domain length
Rc_x  = V./vx;                                                             % particle settling number (the ratios between the velocities)
Rc_f  = V./vf;                                                             % particle settling number
Ne    = xie./V;                                                              % eddy noise flux number
Ns_x  = xix./vx;                                                           % settling noise flux number
Ns_f  = xif./vf;                                                           % settling noise flux number
Pr     = (eta./rho)./((kT+kh.*T)./rho./cP);
Sc     = (eta./rho)./( kc                );

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx ,grdTz ] = gradient(T(icz,icx),h);
    diss = kT./T.*(grdTz (2:end-1,2:end-1).^2 + grdTx (2:end-1,2:end-1).^2) ...
         + exx.*txx + ezz.*tzz ...
         + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
            .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
         + chi.^2.*etas_x./dx0^2 .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
         + phi.^2.*etas_f./df0^2 .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2;
end

%% update time step
dtk = (h/2)^2/max([kc(:);k_x(:);k_f(:);(kT(:)+kh(:).*T(:))./rho(:)./cP(:)]); % diffusive time step size
dta =  h/2   /max(abs([Um(:);Wm(:);Ux(:);Wx(:);Uf(:);Wf(:)]+eps)); % advective time step size
dt  = (dt + min([1.1*dto,min(CFL*[dtk,dta]),dtmax,tau_T/10]))/2;                         % time step size

UDtime = UDtime + toc;