% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

Ktrc    = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc  
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* cx_mem(:,:,j)./100; end
end  

% update trace element phase compositions
trcmq = trc./(m + x.*Ktrc);
trcxq = trc./(m./Ktrc + x);

% get trace element advection
adv_TRC = - advect(M.*trcm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
          - advect(X.*trcx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

% major component diffusion (regularisation)
[diff_TRCm,qz_dffn_TRCm,qx_dffn_TRCm] = diffus(m.*trcm,rho.*kc,h,[1,2],BCD);
[diff_TRCx,qz_dffn_TRCx,qx_dffn_TRCx] = diffus(x.*trcx,rho.*kc,h,[1,2],BCD);

diff_TRC    = diff_TRCm + diff_TRCx;

qz_dffn_TRC = qz_dffn_TRCm + qz_dffn_TRCx;
qx_dffn_TRC = qx_dffn_TRCm + qx_dffn_TRCx;

% get trace element assimilation
bnd_TRC = zeros(size(TRC));
if ~isnan(trcwall(1)); bnd_TRC = bnd_TRC + (permute(repmat(trcwall(1,:).',1,Nz,Nx),[2,3,1]).*rho-TRC)./(tau_a+dt) .* topshape; end
if ~isnan(trcwall(2)); bnd_TRC = bnd_TRC + (permute(repmat(trcwall(2,:).',1,Nz,Nx),[2,3,1]).*rho-TRC)./(tau_a+dt) .* botshape; end
if ~isnan(trcwall(3)); bnd_TRC = bnd_TRC + (permute(repmat(trcwall(3,:).',1,Nz,Nx),[2,3,1]).*rho-TRC)./(tau_a+dt) .* sdsshape; end

% get total rate of change
dTRCdt = adv_TRC + dff_TRC + bnd_TRC;

% residual of trace element evolution
res_TRC = (a1*TRC-a2*TRCo-a3*TRCoo)/dt - (b1*dTRCdt + b2*dTRCdto + b3*dTRCdtoo);

% semi-implicit update of trace element density
[TRC,GHST.TRC,FHST.TRC,specrad.TRC] = iterate(TRC,res_TRC*dt/a1,specrad.TRC,GHST.TRC,FHST.TRC,itpar,iter);

% convert from densites to concentrations
trc = TRC./sum(PHS,3);
