%% Mtap --- Tong, last modified 21 May 2026
% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter==1 || restart
    % generate white noise
    re  = randn(Nz+0, Nx+0);
    rsx = randn(Nz+0, Nx+0);
    rsf = randn(Nz+0, Nx+0);
end


% noise decorrelation time
taue  = L0 ./(V +eps);
tausx = l0x./(vx+eps);
tausf = l0f./(vf+eps);
Stx   = txi0x./taue;   % Stokes number for crystals
Stf   = txi0f./taue;   % Stokes number for fluid (?)

% noise flux amplitudes
sge  = Xi * sqrt(     fReL.*ke  ./taue );                     % eddy mixture noise speed
sgex = Xi * sqrt(chi.*fReL.*ke  ./taue  .*Stx./(1+Stx.^2));   % eddy crystal noise speed
sgef = Xi * sqrt(phi.*fReL.*ke  ./taue  .*Stf./(1+Stf.^2));   % eddy fluid noise speed
sgsx = Xi * sqrt(chi.*      ks_x./tausx);                     % crystal settling noise speed
sgsf = Xi * sqrt(phi.*      ks_f./tausf);                     % fluid settling noise speed

% Ornstein–Uhlenbeck time update for evolving noise
taue  = (1./taue  + 1./(1e2*dt)).^-1 + dt;
tausx = (1./tausx + 1./(1e2*dt)).^-1 + dt;
tausf = (1./tausf + 1./(1e2*dt)).^-1 + dt;
 
Fte   = exp(-dt./taue );                                          % eddy noise time evolution factor
Ftsx  = exp(-dt./tausx);                                          % crystal settling noise time evolution factor
Ftsf  = exp(-dt./tausf);                                          % fluid settling noise time evolution factor
 
fL    = 2/sqrt(1 - exp(-h^2/(2*L0h ^2)));                         % scaling factor (eddy)
flx   = 2/sqrt(1 - exp(-h^2/(2*l0xh^2)));                         % scaling factor (crystal settling)
flf   = 2/sqrt(1 - exp(-h^2/(2*l0fh^2)));                         % scaling factor (fluid settling)
 
psie  = Fte  .* psieo  + sqrt(1 - Fte .^2) .* sge .*fL .* re;    % eddy mixture noise stream function
psiex = Fte  .* psiexo + sqrt(1 - Fte .^2) .* sgex.*fL .* re;    % eddy crystal noise potential
psief = Fte  .* psiefo + sqrt(1 - Fte .^2) .* sgef.*fL .* re;    % eddy fluid noise potential
psisx = Ftsx .* psisxo + sqrt(1 - Ftsx.^2) .* sgsx.*flx.* rsx;    % crystal settling noise potential
psisf = Ftsf .* psisfo + sqrt(1 - Ftsf.^2) .* sgsf.*flf.* rsf;    % fluid settling noise potential
 

% filter noise stream function and potentials spatially to decorrelation length
% pad top/bot to avoid periodic bleed over on non-periodic boundaries
psie_padL0  = zeros(Nz+padL0 ,Nx);
psix_padL0  = zeros(Nz+padL0 ,Nx);
psif_padL0  = zeros(Nz+padL0 ,Nx);
psisx_padl0 = zeros(Nz+padl0x,Nx);
psisf_padl0 = zeros(Nz+padl0f,Nx);
 
psie_padL0 (1+padL0 /2:end-padL0 /2,:) = psie;
psix_padL0 (1+padL0 /2:end-padL0 /2,:) = psiex;
psif_padL0 (1+padL0 /2:end-padL0 /2,:) = psief;
psisx_padl0(1+padl0x/2:end-padl0x/2,:) = psisx;
psisf_padl0(1+padl0f/2:end-padl0f/2,:) = psisf;
 
psie_padL0  = real(ifft2(Gkpe_padL0  .* fft2(psie_padL0 )));
psix_padL0  = real(ifft2(Gkpe_padL0  .* fft2(psix_padL0 )));
psif_padL0  = real(ifft2(Gkpe_padL0  .* fft2(psif_padL0 )));
psisx_padl0 = real(ifft2(Gkps_padl0x .* fft2(psisx_padl0)));
psisf_padl0 = real(ifft2(Gkps_padl0f .* fft2(psisf_padl0)));
 
psie_flt  = psie_padL0 (1+padL0 /2:end-padL0 /2,:);
psiex_flt = psix_padL0 (1+padL0 /2:end-padL0 /2,:);
psief_flt = psif_padL0 (1+padL0 /2:end-padL0 /2,:);
psisx_flt = psisx_padl0(1+padl0x/2:end-padl0x/2,:);
psisf_flt = psisf_padl0(1+padl0f/2:end-padl0f/2,:);

% scale back to mean and std of raw noise fields
psie_flt  = (psie_flt -mean(psie_flt(:) )).*std(psie(:) )./(std(psie_flt(:) )+eps) + mean(psie(:) );
psiex_flt = (psiex_flt-mean(psiex_flt(:))).*std(psiex(:))./(std(psiex_flt(:))+eps) + mean(psiex(:));
psief_flt = (psief_flt-mean(psief_flt(:))).*std(psief(:))./(std(psief_flt(:))+eps) + mean(psief(:));
psisx_flt = (psisx_flt-mean(psisx_flt(:))).*std(psisx(:))./(std(psisx_flt(:))+eps) + mean(psisx(:));
psisf_flt = (psisf_flt-mean(psisf_flt(:))).*std(psisf(:))./(std(psisf_flt(:))+eps) + mean(psisf(:));

% take derivative of noise stream function and potentials to get flux components
 
% eddy mixture noise components (NO taper: this perturbs the whole fluid)
psiec = (psie_flt(icz(1:end-1),icx(1:end-1))+psie_flt(icz(1:end-1),icx(2:end)) ...
      +  psie_flt(icz(2:end  ),icx(1:end-1))+psie_flt(icz(2:end  ),icx(2:end)))/4;
xieu  =  ddz(psiec,1); xieu = xieu(icz,:);
xiew  = -ddx(psiec,1); xiew = xiew(:,icx);

% Taper function 
% taper noise potential fields towards chi=0 and chi=0.6 (for crystal phase)
xtaperw = (1-exp(-chiw./1e-4)).*(1-exp(-max(0,muw-0.4)./0.05));
xtaperu = (1-exp(-chiu./1e-4)).*(1-exp(-max(0,muu-0.4)./0.05));
 
% taper noise potential fields towards phi=0 and phi=0.6 (for fluid phase)
ftaperw = (1-exp(-phiw./1e-4)).*(1-exp(-max(0,muw-0.4)./0.05));
ftaperu = (1-exp(-phiu./1e-4)).*(1-exp(-max(0,muu-0.4)./0.05));
 
% particle eddy / settling noise components.
% Per init.m, xi*w are [Nz+1 x Nx+2] and xi*u are [Nz+2 x Nx+1] (ghost cols
% included). ddx/ddz on psi(icz,icx) [Nz+2 x Nx+2] give exactly those sizes.
% The tapers, however, come from chiw/muw etc. which are [Nz+1 x Nx] (*w) and
% [Nz+2 x Nx+1] (*u); the *w tapers therefore need x-ghost columns padded
% (xtaperw(:,icx)) to multiply the [Nz+1 x Nx+2] derivative. Without this the
% .* failed in 2-D (Nx>1); 1-D hid it because Nx=1.
xtaperw_g = xtaperw(:,icx);   % [Nz+1 x Nx+2]  ghost-padded crystal *w taper
ftaperw_g = ftaperw(:,icx);   % [Nz+1 x Nx+2]  ghost-padded fluid   *w taper

% particle eddy noise components (crystal)  % times taper（ xtaper, ftaper）
xiexu = -ddx(psiex_flt(icz,icx),1).*xtaperu;
xixew = -ddz(psiex_flt(icz,icx),1).*xtaperw_g;

% particle eddy noise components (fluid)
xiefu = -ddx(psief_flt(icz,icx),1).*ftaperu;
xiefw = -ddz(psief_flt(icz,icx),1).*ftaperw_g;

% particle settling noise components (crystal)
xisxu = -ddx(psisx_flt(icz,icx),1).*xtaperu;
xisxw = -ddz(psisx_flt(icz,icx),1).*xtaperw_g;

% particle settling noise components (fluid)
xisfu = -ddx(psisf_flt(icz,icx),1).*ftaperu;
xisfw = -ddz(psisf_flt(icz,icx),1).*ftaperw_g;
 
% combine particle noise
% xi*w are [Nz+1 x Nx+2], xi*u are [Nz+2 x Nx+1] (per init.m) -- no padding.
xixw = (xisxw + xixew);
xixu = (xisxu + xiexu);
xifw = (xisfw + xiefw);
xifu = (xisfu + xiefu);

% get corresponding melt flux (mass conservation: melt compensates for both crystal and fluid)
ximw = -x_w(:,icx)./m_w(:,icx).*xixw - f_w(:,icx)./m_w(:,icx).*xifw;
ximu = -x_u(icz,:)./m_u(icz,:).*xixu - f_u(icz,:)./m_u(icz,:).*xifu;










% % filter noise amplitudes spatially to decorrelation length
% % pad top/bot to avoid periodic bleed over on non-periodic boundaries
% sge_padL0  = zeros(Nz+padL0 ,Nx); 
% sgex_padL0 = zeros(Nz+padL0 ,Nx);
% sgef_padL0 = zeros(Nz+padL0 ,Nx);
% sgsx_padl0 = zeros(Nz+padl0x,Nx);
% sgsf_padl0 = zeros(Nz+padl0f,Nx);
% 
% sge_padL0(1+padL0/2:end-padL0/2,:) = sge_raw;
% sge_padL0(            1:padL0/2,:) = repmat(sge_raw(1  ,:),padL0/2,1);
% sge_padL0(end-padL0/2+1:end    ,:) = repmat(sge_raw(end,:),padL0/2,1);
% 
% sgex_padL0(1+padL0/2:end-padL0/2,:) = sgex_raw;
% sgex_padL0(            1:padL0/2,:) = repmat(sgex_raw(1  ,:),padL0/2,1);
% sgex_padL0(end-padL0/2+1:end    ,:) = repmat(sgex_raw(end,:),padL0/2,1);
% 
% sgef_padL0(1+padL0/2:end-padL0/2,:) = sgef_raw;
% sgef_padL0(            1:padL0/2,:) = repmat(sgef_raw(1  ,:),padL0/2,1);
% sgef_padL0(end-padL0/2+1:end    ,:) = repmat(sgef_raw(end,:),padL0/2,1);
% 
% sgsx_padl0(1+padl0x/2:end-padl0x/2,:) = sgsx_raw;
% sgsx_padl0(             1:padl0x/2,:) = repmat(sgsx_raw(1  ,:),padl0x/2,1);
% sgsx_padl0(end-padl0x/2+1:end     ,:) = repmat(sgsx_raw(end,:),padl0x/2,1);
% 
% sgsf_padl0(1+padl0f/2:end-padl0f/2,:) = sgsf_raw;
% sgsf_padl0(             1:padl0f/2,:) = repmat(sgsf_raw(1  ,:),padl0f/2,1);
% sgsf_padl0(end-padl0f/2+1:end     ,:) = repmat(sgsf_raw(end,:),padl0f/2,1);
% 
% sge_padL0  = real(ifft2(Gkpe_padL0 .* fft2(sge_padL0)));
% sgex_padL0 = real(ifft2(Gkpe_padL0 .* fft2(sgex_padL0)));
% sgef_padL0 = real(ifft2(Gkpe_padL0 .* fft2(sgef_padL0)));
% sgsx_padl0 = real(ifft2(Gkps_padl0x.* fft2(sgsx_padl0)));
% sgsf_padl0 = real(ifft2(Gkps_padl0f.* fft2(sgsf_padl0)));
% 
% sge  = sge_padL0 (1+padL0 /2:end-padL0 /2,:);
% sgex = sgex_padL0(1+padL0 /2:end-padL0 /2,:);
% sgef = sgef_padL0(1+padL0 /2:end-padL0 /2,:);
% sgsx = sgsx_padl0(1+padl0x/2:end-padl0x/2,:); 
% sgsf = sgsf_padl0(1+padl0f/2:end-padl0f/2,:); 
% 
% % scale back to mean and std of raw noise amplitude fields
% sge  = (sge -mean(sge(:)) ).*std(sge_raw(:) )./(std(sge(:) )+eps) + mean(sge_raw(:) );
% sgex = (sgex-mean(sgex(:))).*std(sgex_raw(:))./(std(sgex(:))+eps) + mean(sgex_raw(:));
% sgef = (sgef-mean(sgef(:))).*std(sgef_raw(:))./(std(sgef(:))+eps) + mean(sgef_raw(:));
% sgsx = (sgsx-mean(sgsx(:))).*std(sgsx_raw(:))./(std(sgsx(:))+eps) + mean(sgsx_raw(:));
% sgsf = (sgsf-mean(sgsf(:))).*std(sgsf_raw(:))./(std(sgsf(:))+eps) + mean(sgsf_raw(:));
% 
% 
% % Ornstein–Uhlenbeck time update for evolving noise
% Fte   = exp(-dt./taue );                                           % eddy noise time evolution factor
% Ftsx  = exp(-dt./tausx);                                           % settling noise time evolution factor
% Ftsf  = exp(-dt./tausf);                                           % settling noise time evolution factor
% 
% fL    = 2/sqrt(1 - exp(-h^2/(2*L0h ^2)));                          % scaling factor for potential field to noise component variance
% flx   = 2/sqrt(1 - exp(-h^2/(2*l0xh^2)));                          % scaling factor for potential field to noise component variance
% flf   = 2/sqrt(1 - exp(-h^2/(2*l0fh^2)));                          % scaling factor for potential field to noise component variance
% 
% psie  = Fte  .* psieo  + sqrt(1 - Fte .^2) .* sge .*fL .* re ;     % update eddy noise stream function
% psix  = Fte  .* psixo  + sqrt(1 - Fte .^2) .* sgex.*fL .* rex;     % update particle eddy noise potential
% psif  = Fte  .* psifo  + sqrt(1 - Fte .^2) .* sgef.*fL .* ref;     % update particle eddy noise potential
% psisx = Ftsx .* psisxo + sqrt(1 - Ftsx.^2) .* sgsx.*flx.* rsx;     % update particle settling noise potential
% psisf = Ftsf .* psisfo + sqrt(1 - Ftsf.^2) .* sgsf.*flf.* rsf;     % update particle settling noise potential


% % take derivative of noise stream function and potentials to get flux components
% 
% % eddy noise components 
% psiec = (psie(icz(1:end-1),icx(1:end-1))+psie(icz(1:end-1),icx(2:end)) ...
%       +  psie(icz(2:end  ),icx(1:end-1))+psie(icz(2:end  ),icx(2:end)))/4;
% xieu  =  ddz(psiec,1); xieu = xieu(icz,:);
% xiew  = -ddx(psiec,1); xiew = xiew(:,icx);

% % particle eddy noise components
% xixu = -ddx(psix(icz,icx),1);
% xixw = -ddz(psix(icz,icx),1);
% xifu = -ddx(psif(icz,icx),1);
% xifw = -ddz(psif(icz,icx),1);
% 
% % particle settling noise components
% xisxu = -ddx(psisx (icz,icx),1);
% xisxw = -ddz(psisx (icz,icx),1);
% xisfu = -ddx(psisf (icz,icx),1);
% xisfw = -ddz(psisf (icz,icx),1);
% 
% % combine particle noise
% xiwx = (xisxw + xixw);
% xiux = (xisxu + xixu);
% xiwf = (xisfw + xifw);
% xiuf = (xisfu + xifu);
% 
% % get corresponding melt flux
% xiwm = -x_w(:,icx)./m_w(:,icx).*xiwx - f_w(:,icx)./m_w(:,icx).*xiwf;
% xium = -x_u(icz,:)./m_u(icz,:).*xiux - f_u(icz,:)./m_u(icz,:).*xiuf;



