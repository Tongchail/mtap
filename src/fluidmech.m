tic;

if ~bnchm && step>0 && ~restart

%***  update mixture mass density
drhodt  = advn_rho;             % advection term

% residual of mixture mass evolution
res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);

% volume source and background velocity passed to fluid-mechanics solver
[MFS,GHST.MFS,FHST.MFS,specrad.MFS] = iterate(MFS,res_rho./b1,specrad.MFS,GHST.MFS,FHST.MFS,itpar,iter*~frst);

MFSmean  = mean(MFS,'all');

MFBG     = MFSmean .* ZZw;

end

%% 0-D run does not require fluidmech solve
if Nz==1 && Nx==1  
    W  = WBG; Wm = W;  Wx = W;  Wf = W;
    U  = UBG; Um = U;  Ux = U;  Uf = U;
    P  = zeros(Nz+2,Nx+2);
    resnorm_VP = 0;
else

if Nx==1
    % update 1D velocity
    W(:,2) = -flipud(cumsum(flipud([MFS*h;-WBG(end)])));

    % update 1D pressure
    Div_tz = ddz(tzz(icz,:),h);           % z-stress divergence
    PrsSrc = - Div_tz;
    P(:,2) =  flipud(cumsum(flipud([PrsSrc*h;0])));
    P(:,2) = P(:,2) - P(Nz/2,2);

    W(:,1) = W(:,2); W(:,end) = W(:,2);
    P(:,1) = P(:,2); P(:,end) = P(:,2);
else


%% assemble coefficients for matrix velocity diagonal and right-hand side

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R


% assemble coefficients of z-stress divergence
    
% left boundary
if periodic
    ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,end-1);
else
    ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);    % keep at the moment
end
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
if periodic
    ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,2);
else
    ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
end
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% top boundary
ii  = MapW(1,2:end-1); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii)) + WBG(1,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary (*xcore replace)   the operator for divergence of velocity → the divergence for the mass flux
ii  = MapW(end,2:end-1); jj1 = ii; jj2 = MapW(end-1,2:end-1); jj3 = MapU(end-1,2:end); jj4 = MapU(end-1,1:end-1);
rho1 = rhow(end  ,:    );
rho2 = rhow(end-1,:    );
rho3 = rhou(end,2:end  );
rho4 = rhou(end,1:end-1);
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;+         rho1(:)/h];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-open_cnv*rho2(:)/h];
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+open_cnv*rho3(:)/h];
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-open_cnv*rho4(:)/h];
aa  = open_cnv.*MFS(end,:) + (1-open_cnv).*MFBG(end, 2:end-1)/h;
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
% ii  = MapW(end,2:end-1); jj = ii;
% aa  = zeros(size(ii));
% IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
% aa  = zeros(size(ii)) + WBG(end,2:end-1);
% IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 =  etaco(2:end-1,1:end-1);   EtaC2 =  etaco(2:end-1,2:end);
EtaP1 =  eta  (1:end-1,:      );   EtaP2 =  eta  (2:end,:      );

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa  = a1.*rhow(2:end-1,:)./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = ddz(rho,h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right


% z-RHS vector
advn_mz = advect(rhow(2:end-1,:).*W(2:end-1,2:end-1),(U(2:end-2,:)+U(3:end-1,:))/2,(W(1:end-1,2:end-1)+W(2:end,2:end-1))/2,h,{ADVN,''},[1,2],BCA);
rr  = + Drho(2:end-1,:) .* g0 ...
      + (a2.*rhoWo(2:end-1,:)+a3.*rhoWoo(2:end-1,:))/dt ...
      - advn_mz;
if bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+top_cnv];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+bot_cnv];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

if ~periodic
    % left boundary
    ii  = MapU(2:end-1,1); jj = ii;
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
    aa  = zeros(size(ii));
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

    % right boundary
    ii  = MapU(2:end-1,end); jj = ii;
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
    aa  = zeros(size(ii));
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
end

% internal points
if periodic
    ii    = MapU(2:end-1,:);
    EtaC1 = etaco(1:end-1,:    );  EtaC2 = etaco(2:end,:    );
    EtaP1 = eta  (:,icx(1:end-1));  EtaP2 = eta  (:,icx(2:end));
else
    ii    = MapU(2:end-1,2:end-1);
    EtaC1 = etaco(1:end-1,2:end-1);  EtaC2 = etaco(2:end,2:end-1);
    EtaP1 = eta  (:      ,1:end-1);  EtaP2 = eta  (:      ,2:end);
end

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
if periodic
    jj1 = MapU(2:end-1,ifx(1:end-2)); jj2 = MapU(2:end-1,ifx(3:end)); jj3 = MapU(1:end-2,ifx(2:end-1)); jj4 = MapU(3:end,ifx(2:end-1));
else
    jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);
end
if periodic
    aa  = (a1+gamma).*rhou./dt;
else
    aa  = a1.*rhou(:,2:end-1)./dt;
end
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % U one below

% what shall we do with the drunken sailor...
if ~bnchm
    if periodic
        aa  = ddx(rho(:,icx),h).*g0.*dt;
    else
        aa  = ddx(rho,h).*g0.*dt;
    end
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
if periodic
    jj1 = MapW(1:end-1,1:end-1); jj2 = MapW(1:end-1,2:end); jj3 = MapW(2:end,1:end-1); jj4 = MapW(2:end,2:end);
else
    jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);
end
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
if periodic
    advn_mx = advect(rhou.*U(2:end-1,:),(U(2:end-1,ifx(1:end-1))+U(2:end-1,ifx(2:end)))/2,(W(:,1:end-1)+W(:,2:end))/2,h,{ADVN,''},[1,2],BCA);
    advn_mx(:,[1 end]) = repmat((advn_mx(:,1)+advn_mx(:,end))/2,1,2);
    rr  = + (a2.*rhoUo+a3.*rhoUoo)/dt ...
          - advn_mx;
else
    advn_mx = advect(rhou(:,2:end-1).*U(2:end-1,2:end-1),(U(2:end-1,1:end-1)+U(2:end-1,2:end))/2,(W(:,2:end-2)+W(:,3:end-1))/2,h,{ADVN,''},[1,2],BCA);
    rr  = + (a2.*rhoUo(:,2:end-1)+a3.*rhoUoo(:,2:end-1))/dt ...
          - advn_mx;
end
if bnchm
    if periodic
        rr = rr + src_U_mms(2:end-1,:); 
    else
        rr = rr + src_U_mms(2:end-1,2:end-1); 
    end
end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV  = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RV  = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for gradient operator

if ~exist('GG','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    % coefficients for z-gradient
    ii  = MapW(2:end-1,2:end-1);
    
    %         top              ||          bottom
    jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
    
    % coefficients for x-gradient
    if periodic
        ii  = MapU(2:end-1,:);
    else
        ii  = MapU(2:end-1,2:end-1);
    end
    
    %         left             ||           right
    if periodic
        jj1 = MapP(2:end-1,1:end-1); jj2 = MapP(2:end-1,2:end);
    else
        jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
    end
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
    
    % assemble coefficient matrix
    GG  = sparse(IIL,JJL,AAL,NW+NU,NP);
end


%% assemble coefficients for divergence of matrix mass flux (DM)

IIL = [];       % equation indeces into A
JJL = [];       % variable indeces into A
AAL = [];       % coefficients for A

%internal points
ii  = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
rho1 = rhou(:,1:end-1);
rho2 = rhou(:,2:end  );
rho3 = rhow(1:end-1,:);
rho4 = rhow(2:end  ,:);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; -rho1(:)/h];  % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; +rho2(:)/h];  % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; -rho3(:)/h];  % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; +rho4(:)/h];  % W one below

% Assemble coefficient matrix
DM  = sparse(IIL,JJL,AAL,NP,NW+NU);

%% assemble coefficients for matrix pressure diagonal and right-hand side

if ~exist('KP','var') || bnchm

    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    % boundary points
    ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
    jj1 = ii;
    jj2 = [MapP(2,:).'; MapP(end-1,:).'];
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    ii  = [MapP(2:end-1,1    ); MapP(2:end-1,end)]; % left & right
    jj1 = ii;
    if periodic
        jj2 = [MapP(2:end-1,end-1); MapP(2:end-1,2    )];
    else
        jj2 = [MapP(2:end-1,2    ); MapP(2:end-1,end-1)];
    end
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
      
% if ~exist('KP','var') || bnchm || lambda1+lambda2>0
% 
%     % internal points
%     ii  = MapP(2:end-1,2:end-1);
%     jj1 = MapP(1:end-2,2:end-1);
%     jj2 = MapP(3:end-0,2:end-1);
%     jj3 = MapP(2:end-1,1:end-2);
%     jj4 = MapP(2:end-1,3:end-0);
% 
%     % coefficients multiplying matrix pressure P
%     aa  = zeros(size(ii)) + lambda1*eps*h^2./eta;
%     IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre
% 
%     kP  = lambda2*h^2./eta;
%     kP1 = (kP(icz(1:end-2),:).*kP(icz(2:end-1),:)).^0.5;   kP2 = (kP(icz(2:end-1),:).*kP(icz(3:end-0),:)).^0.5;
%     kP3 = (kP(:,icx(1:end-2)).*kP(:,icx(2:end-1))).^0.5;   kP4 = (kP(:,icx(2:end-1)).*kP(:,icx(3:end-0))).^0.5;
% 
%     aa  = (kP1+kP2+kP3+kP4)/h^2;
%     IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL;-aa(:)     ];      % P on stencil centre
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; kP1(:)/h^2];      % P one above
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; kP2(:)/h^2];      % P one below
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; kP3(:)/h^2];      % P one above
%     IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; kP4(:)/h^2];      % P one below
% 
% end

% assemble coefficient matrix
KP  = sparse(IIL,JJL,AAL,NP,NP);

end

% RHS
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

ii  = MapP(2:end-1,2:end-1);

rr  = MFS;       % add volume source term
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RP  = sparse(IIR,ones(size(IIR)),AAR,NP,1);

% set P = 0 in fixed point
% nzp = 2; %round((Nz+2)/2);
% nxp = round((Nx+2)/2);
% np0 = MapP(nzp,nxp);
% % nxp = 2:Nx-1;
% % np0 = MapP(nzp,nxp);
% KP(np0,:  ) = 0;
% KP(np0,np0) = 1;%KP(np0,np0) + 1e-6.*h^2./geomean(eta(:));
% DD(np0,:  ) = 0;
% RP(np0    ) = 0;
% 
% if bnchm; RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp); end

if bnchm
    % fix P = P_mms in middle of domain
    nzp = round((Nz+2)*3/8);
    nxp = round((Nx+2)/2);
    np0 = MapP(nzp,nxp);
    KP(np0,:  ) = 0;
    KP(np0,np0) = 1;
    DD(np0,:  ) = 0;
    RP(np0    ) = P_mms(nzp,nxp);
    
    % fix U = U_mms in middle of domain
    nzu = round((Nz+2)/2);
    nxu = round((Nx+2)/2);
    nu0 = MapU(nzu,nxu);
    KV(nu0,:  ) = 0;
    KV(nu0,nu0) = 1;
    GG(nu0,:  ) = 0;
    RV(nu0    ) = U_mms(nzu,nxu);
else
    % set P = 0 in fixed point
    if open_cnv
        nzp = Nz+1;
        nxp = 1:Nx+2;
    else
        nzp = round(Nz/2);
        nxp = round(Nx/2);
    end
    np0 = MapP(nzp,nxp);
    KP(np0,:  ) = 0;
    KP(np0,np0) = speye(length(np0));
end

% if bnchm
%     % fix P = P_mms in middle of domain
%     nzp = 1;%round((Nz+2)/2);
%     nxp = round((Nx+2)/2);
%     np0 = MapP(nzp,nxp);
%     DD(MapP(nzp,nxp),:) = 0;
%     KP(MapP(nzp,nxp),:) = 0;
%     KP(MapP(nzp,nxp),MapP(nzp,nxp)) = 1;
%     RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp);
% 
%     % fix U = U_mms in middle of domain
%     nzu = round((Nz+2)/2);
%     nxu = round((Nx+2)/2);
%     KV(MapU(nzu,nxu),:) = 0;
%     GG(MapU(nzu,nxu),:) = 0;
%     KV(MapU(nzu,nxu),MapU(nzu  ,nxu)) = 1;
%     RV(MapU(nzp,nxp),:) = U_mms(nzu,nxu);
% end

%% assemble and scale global coefficient matrix and right-hand side vector
% differences in how the matrix is scaled (because the mass flux divergence)
LL  = [KV GG  ; ...
       DM KP ];

RR  = [RV; RP;];

% etagh = ones(size(P));  etagh(2:end-1,2:end-1) = eta;
scl = ones(size(P));  scl(2:end-1,2:end-1) = rho./eta;
SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse( 1./(SCL + sqrt([zeros(NU+NW,1); 1./scl(:)])) ));

% FF  = LL*[W(:);U(:);P(:)] - RR;
% FF  = SCL*FF;

FF  = SCL*(LL*SOL - RR);
LL  = SCL*LL*SCL;
% RR  = SCL*RR;


%% Solve linear system of equations for W, U, P

if ~exist('pcol','var') || bnchm; pcol = colamd(LL); end % get column permutation for sparsity pattern once per run
dLL         = decomposition(LL(:,pcol), 'lu');  % get LU-decomposition for consistent performance of LL \ RR
UPD(pcol,1) = dLL \ FF;                         % solve permuted decomposed system
UPD         = SCL*UPD;

% map solution vector to 2D arrays
upd_W = -full(reshape(UPD(MapW(:))        ,Nz+1,Nx+2));  % matrix z-velocity
upd_U = -full(reshape(UPD(MapU(:))        ,Nz+2,Nx+1));  % matrix x-velocity
upd_P = -full(reshape(UPD(MapP(:)+(NW+NU)),Nz+2,Nx+2));  % matrix dynamic pressure

% update solution
W = W + upd_W;
U = U + upd_U;
P = P + upd_P;
SOL = [W(:);U(:);P(:)];


end

%% Update phase segregation speeds
if ~bnchm && step>=1

    % taper towards boundaries if closed segregation top/bot boundaries
    bndtaperwx = (1 - (exp((-ZZw)/l0xh) + exp(-(D-ZZw)/l0xh)).*(1-open_sgr));
    bndtaperwf = (1 - (exp((-ZZw)/l0fh) + exp(-(D-ZZw)/l0fh)).*(1-open_sgr));

    % terminal xtal segregation speed  
    wx(:,2:end-1) = dx0^2./etas_xw.*Drhox.*g0;  %*take the segregation coefficient apart to address the viscosity
    wx = wx .* bndtaperwx;
    wx(:,[1 end]) = wx(:,[end-1 2]);

    % fluid segregation speed 
    wf(:,2:end-1) = df0^2./etas_fw.*Drhof.*g0;
    wf = wf .* bndtaperwf;
    wf(:,[1 end]) = wf(:,[end-1 2]);
   
    % melt segregation speed   (mass fraction, xw​*wx ​+ fw*​wf​ + mw*​wm​ = 0  mass flux)
    wm  = -xw(:,icx)./mw(:,icx).*wx - ffw(:,icx)./mw(:,icx).*wf;

    % phase diffusion fluxes and speeds
    [~,wqx,uqx] = diffus(chi,k_x,h,[1,2],BCD);
    [~,wqf,uqf] = diffus(phi,k_f,h,[1,2],BCD);
    wqx  = wqx .* bndtaperwx;
    wqf  = wqf .* bndtaperwf;
    wqm  = -xw(:,icx)./mw (:,icx).*wqx - ffw(:,icx)./mw (:,icx).*wqf;
    uqm  = -xu(icz,:)./muu(icz,:).*uqx - fu (icz,:)./muu(icz,:).*uqf;

    % wqm = -wqx-wqf;
    % uqm = -uqx-uqf;

    % update stochastic noise speeds
    noise;  

    % update phase velocities
    Wx = W + wx + wqx + xiwx + xiew; % xtl z-velocity
    Ux = U + 0. + uqx + xiux + xieu; % xtl x-velocity
    
    Wf = W + wf + wqf + xiwf + xiew; % mfe z-velocity 
    Uf = U + 0. + uqf + xiuf + xieu; % mfe x-velocity

    Wm = W + wm + wqm + xiwm + xiew; % mlt z-velocity
    Um = U + 0. + uqm + xium + xieu; % mlt x-velocity

end

end

FMtime = FMtime + toc;