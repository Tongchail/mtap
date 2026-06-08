% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TRC = zeros(Nz,Nx,cal.ntrc);
Ktrc    = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* cx_mem(:,:,j)./100; end

    % update trace element phase compositions
    % guard denominators against division by zero: Ktrc can be 0 if the only mineral phases present have zero partitioning for element i, and m/x can vanish in pure-solid/pure-melt limits (esp. 0-D runs).
    Ktrc_safe   = max(Ktrc(:,:,i),eps);
    trcm(:,:,i) = trc(:,:,i)./max(m + x.*Ktrc_safe, eps);
    trcx(:,:,i) = trc(:,:,i)./max(m./Ktrc_safe + x, eps);

    % get trace element advection
    adv_TRC(:,:,i) = - advect(M.*trcm(:,:,i),Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
                     - advect(X.*trcx(:,:,i),Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

    % get trace element diffusion (regularisation)
    dff_TRC(:,:,i) = diffus(trcm(:,:,i),M.*kc,h,[1,2],BCD) + diffus(trcx(:,:,i),X.*kc,h,[1,2],BCD);

    % get trace element assimilation
    if ~isnan(trcwall(1,i)); bnd_TRC(:,:,i) = bnd_TRC(:,:,i) + (rho.*trcwall(1,i)-TRC(:,:,i)).*mu./tau_a .* topshape; end
    if ~isnan(trcwall(2,i)); bnd_TRC(:,:,i) = bnd_TRC(:,:,i) + (rho.*trcwall(2,i)-TRC(:,:,i)).*mu./tau_a .* botshape; end
    if ~isnan(trcwall(3,i)); bnd_TRC(:,:,i) = bnd_TRC(:,:,i) + (rho.*trcwall(3,i)-TRC(:,:,i)).*mu./tau_a .* sdsshape; end
end

% get total rate of change
dTRCdt = adv_TRC + dff_TRC + bnd_TRC;

% residual of trace element evolution
res_TRC = (a1*TRC-a2*TRCo-a3*TRCoo)/dt - (b1*dTRCdt + b2*dTRCdto + b3*dTRCdtoo);

% semi-implicit update of trace element density
[TRC,GHST.TRC,FHST.TRC,specrad.TRC] = iterate(TRC,res_TRC*dt/a1,specrad.TRC,GHST.TRC,FHST.TRC,itpar,iter*~frst);

% convert from densites to concentrations
for i = 1:cal.ntrc; trc(:,:,i) = TRC(:,:,i)./rho; end
