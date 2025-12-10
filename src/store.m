[Soo    , So    ] = deal(So    , S    );
[Coo    , Co    ] = deal(Co    , C    );
[Moo    , Mo    ] = deal(Mo    , M    );
[Xoo    , Xo    ] = deal(Xo    , X    );
[Foo    , Fo    ] = deal(Fo    , F    );
[rhooo  , rhoo  ] = deal(rhoo  , rho  );
[TRCoo  , TRCo  ] = deal(TRCo  , TRC  );
[Too    , To    ] = deal(To    , T    );
[Tpoo   , Tpo   ] = deal(Tpo   , Tp   );
%[Pchmboo, Pchmbo] = deal(Pchmbo, Pchmb);
[rhoWoo , rhoWo ] = deal(rhoWo , rhoW );
[rhoUoo , rhoUo ] = deal(rhoUo , rhoU );

[dSdtoo    , dSdto    ] = deal(dSdto    , dSdt    );
[dCdtoo    , dCdto    ] = deal(dCdto    , dCdt    );
[dMdtoo    , dMdto    ] = deal(dMdto    , dMdt    );
[dXdtoo    , dXdto    ] = deal(dXdto    , dXdt    );
[dFdtoo    , dFdto    ] = deal(dFdto    , dFdt    );
[drhodtoo  , drhodto  ] = deal(drhodto  , drhodt  );
[dTRCdtoo  , dTRCdto  ] = deal(dTRCdto  , dTRCdt  );
[dTpdtoo   , dTpdto   ] = deal(dTpdto   , dTpdt   );
[dTdtoo    , dTdto    ] = deal(dTdto    , dTdt    );
%[dPchmbdtoo, dPchmbdto] = deal(dPchmbdto, dPchmbdt);

[Gmo,Gxo,Gfo] = deal(Gm,Gx,Gf);

psieo = psie;
psixo = psix;
psifo = psif;
psisxo = psisx;
psisfo = psisf;

dto = dt;
Pto = Pt;

% reset update history
specrad.S.est   = specrad.S.mean + 0.*specrad.S.est;
specrad.C.est   = specrad.S.mean + 0.*specrad.C.est;
specrad.PHS.est = specrad.S.mean + 0.*specrad.PHS.est;
specrad.TRC.est = specrad.S.mean + 0.*specrad.TRC.est;
specrad.MFS.est = specrad.S.mean + 0.*specrad.MFS.est;
FHST.S          = 0.*FHST.S;
FHST.C          = 0.*FHST.C;
FHST.PHS        = 0.*FHST.PHS;
FHST.TRC        = 0.*FHST.TRC;
FHST.MFS        = 0.*FHST.MFS;