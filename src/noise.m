% update stochastic noise representing subgrid fluctuations from particle
% settling and turbulent eddies

% generate smooth random noise (once per timestep)
if iter==1 || restart
    % generate white noise
    re  = randn(Nz+0, Nx+0);
    rex = randn(Nz+0, Nx+0);
    ref = randn(Nz+0, Nx+0);
    rsx = randn(Nz+0, Nx+0);
    rsf = randn(Nz+0, Nx+0);

    % filter white noise spatially to decorrelation length
    re  = real(ifft2(Gkpe .* fft2(re)));
    rex = real(ifft2(Gkpe .* fft2(rex)));
    ref = real(ifft2(Gkpe .* fft2(ref)));
    rsx = real(ifft2(Gkpx .* fft2(rsx)));
    rsf = real(ifft2(Gkpf .* fft2(rsf)));

    % rescale filtered noise to zero mean and unit variance
    re  = (re - mean(re(:))) ./ std(re(:));
    rex = (rex - mean(rex(:))) ./ std(rex(:));
    ref = (ref - mean(ref(:))) ./ std(ref(:));
    rsx = (rsx - mean(rsx(:))) ./ std(rsx(:));
    rsf = (rsf - mean(rsf(:))) ./ std(rsf(:));
end


% noise decorrelation time
taue  = L0h ./(V +eps);
tausx = l0xh./(vx+eps);
tausf = l0fh./(vf+eps);

% noise flux amplitudes
sge_raw  = Xi * sqrt(     fReL.*ke  ./taue  .* (L0h ./(L0h +h)).^3);    % eddy mixture noise speed
sgex_raw = Xi * sqrt(chi.*fReL.*ke  ./taue  .* (L0h ./(L0h +h)).^3);    % eddy crystal noise speed
sgef_raw = Xi * sqrt(phi.*fReL.*ke  ./taue  .* (L0h ./(L0h +h)).^3);    % eddy fluid noise speed
sgsx_raw = Xi * sqrt(chi.*      ks_x./tausx .* (l0xh./(l0xh+h)).^3);    % crystal segregation noise speed
sgsf_raw = Xi * sqrt(phi.*      ks_f./tausf .* (l0fh./(l0fh+h)).^3);    % fluid segregation noise speed


% filter noise amplitudes spatially to decorrelation length
% pad top/bot to avoid periodic bleed over on non-periodic boundaries
sge_padL0  = zeros(Nz+padL0 ,Nx); 
sgex_padL0 = zeros(Nz+padL0 ,Nx);
sgef_padL0 = zeros(Nz+padL0 ,Nx);
sgsx_padl0 = zeros(Nz+padl0x,Nx);
sgsf_padl0 = zeros(Nz+padl0f,Nx);

sge_padL0(1+padL0/2:end-padL0/2,:) = sge_raw;
sge_padL0(            1:padL0/2,:) = repmat(sge_raw(1  ,:),padL0/2,1);
sge_padL0(end-padL0/2+1:end    ,:) = repmat(sge_raw(end,:),padL0/2,1);

sgex_padL0(1+padL0/2:end-padL0/2,:) = sgex_raw;
sgex_padL0(            1:padL0/2,:) = repmat(sgex_raw(1  ,:),padL0/2,1);
sgex_padL0(end-padL0/2+1:end    ,:) = repmat(sgex_raw(end,:),padL0/2,1);

sgef_padL0(1+padL0/2:end-padL0/2,:) = sgef_raw;
sgef_padL0(            1:padL0/2,:) = repmat(sgef_raw(1  ,:),padL0/2,1);
sgef_padL0(end-padL0/2+1:end    ,:) = repmat(sgef_raw(end,:),padL0/2,1);

sgsx_padl0(1+padl0x/2:end-padl0x/2,:) = sgsx_raw;
sgsx_padl0(             1:padl0x/2,:) = repmat(sgsx_raw(1  ,:),padl0x/2,1);
sgsx_padl0(end-padl0x/2+1:end     ,:) = repmat(sgsx_raw(end,:),padl0x/2,1);

sgsf_padl0(1+padl0f/2:end-padl0f/2,:) = sgsf_raw;
sgsf_padl0(             1:padl0f/2,:) = repmat(sgsf_raw(1  ,:),padl0f/2,1);
sgsf_padl0(end-padl0f/2+1:end     ,:) = repmat(sgsf_raw(end,:),padl0f/2,1);

sge_padL0  = real(ifft2(Gkpe_padL0 .* fft2(sge_padL0)));
sgex_padL0 = real(ifft2(Gkpe_padL0 .* fft2(sgex_padL0)));
sgef_padL0 = real(ifft2(Gkpe_padL0 .* fft2(sgef_padL0)));
sgsx_padl0 = real(ifft2(Gkps_padl0x.* fft2(sgsx_padl0)));
sgsf_padl0 = real(ifft2(Gkps_padl0f.* fft2(sgsf_padl0)));

sge  = sge_padL0 (1+padL0 /2:end-padL0 /2,:);
sgex = sgex_padL0(1+padL0 /2:end-padL0 /2,:);
sgef = sgef_padL0(1+padL0 /2:end-padL0 /2,:);
sgsx = sgsx_padl0(1+padl0x/2:end-padl0x/2,:); 
sgsf = sgsf_padl0(1+padl0f/2:end-padl0f/2,:); 

% scale back to mean and std of raw noise amplitude fields
sge  = (sge -mean(sge(:)) ).*std(sge_raw(:) )./(std(sge(:) )+eps) + mean(sge_raw(:) );
sgex = (sgex-mean(sgex(:))).*std(sgex_raw(:))./(std(sgex(:))+eps) + mean(sgex_raw(:));
sgef = (sgef-mean(sgef(:))).*std(sgef_raw(:))./(std(sgef(:))+eps) + mean(sgef_raw(:));
sgsx = (sgsx-mean(sgsx(:))).*std(sgsx_raw(:))./(std(sgsx(:))+eps) + mean(sgsx_raw(:));
sgsf = (sgsf-mean(sgsf(:))).*std(sgsf_raw(:))./(std(sgsf(:))+eps) + mean(sgsf_raw(:));


% Ornstein–Uhlenbeck time update for evolving noise
Fte   = exp(-dt./taue );                                           % eddy noise time evolution factor
Ftsx  = exp(-dt./tausx);                                           % settling noise time evolution factor
Ftsf  = exp(-dt./tausf);                                           % settling noise time evolution factor

fL    = 2/sqrt(1 - exp(-h^2/(2*L0h ^2)));                          % scaling factor for potential field to noise component variance
flx   = 2/sqrt(1 - exp(-h^2/(2*l0xh^2)));                          % scaling factor for potential field to noise component variance
flf   = 2/sqrt(1 - exp(-h^2/(2*l0fh^2)));                          % scaling factor for potential field to noise component variance

psie  = Fte  .* psieo  + sqrt(1 - Fte .^2) .* sge .*fL .* re ;     % update eddy noise stream function
psix  = Fte  .* psixo  + sqrt(1 - Fte .^2) .* sgex.*fL .* rex;     % update particle eddy noise potential
psif  = Fte  .* psifo  + sqrt(1 - Fte .^2) .* sgef.*fL .* ref;     % update particle eddy noise potential
psisx = Ftsx .* psisxo + sqrt(1 - Ftsx.^2) .* sgsx.*flx.* rsx;     % update particle settling noise potential
psisf = Ftsf .* psisfo + sqrt(1 - Ftsf.^2) .* sgsf.*flf.* rsf;     % update particle settling noise potential


% take derivative of noise stream function and potentials to get flux components

% eddy noise components
psiec = (psie(icz(1:end-1),icx(1:end-1))+psie(icz(1:end-1),icx(2:end)) ...
      +  psie(icz(2:end  ),icx(1:end-1))+psie(icz(2:end  ),icx(2:end)))/4;
xieu  =  ddz(psiec,1); xieu = xieu(icz,:);
xiew  = -ddx(psiec,1); xiew = xiew(:,icx);

% particle eddy noise components
xixu = -ddx(psix(icz,icx),1);
xixw = -ddz(psix(icz,icx),1);
xifu = -ddx(psif(icz,icx),1);
xifw = -ddz(psif(icz,icx),1);

% particle settling noise components
xisxu = -ddx(psisx (icz,icx),1);
xisxw = -ddz(psisx (icz,icx),1);
xisfu = -ddx(psisf (icz,icx),1);
xisfw = -ddz(psisf (icz,icx),1);

% combine particle noise
xiwx = (xisxw + xixw);
xiux = (xisxu + xixu);
xiwf = (xisfw + xifw);
xiuf = (xisfu + xifu);

% get corresponding melt flux
xiwm = -xw(:,icx)./mw (:,icx).*xiwx - ffw(:,icx)./mw (:,icx).*xiwf;
xium = -xu(icz,:)./muu(icz,:).*xiux - fu(icz,:)./muu(icz,:).*xiuf;



