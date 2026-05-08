%% Copyright (c) Adam Ouzeri 2025

Cxy = ArealDeformation;
Cz = 1./Cxy./Cxy;
kxysquared = Cxy;

% Gsharp for calculations
Gxy = Y(:,1);
Gxylat = Y(:,2);
Gz = Y(:,3);

% Apical
Cmixed_xy = Cxy.*Gxy;
Iap = Cmixed_xy + Cmixed_xy;
Jap = sqrt(Cmixed_xy.*Cmixed_xy);

% lateral parallel side
Cmixed_11lat = Cxy.*Gxylat;
Cmixed_22lat = Cz.*Gz;
Ilat = Cmixed_11lat + Cmixed_22lat;
Jlat = sqrt(Cmixed_11lat.*Cmixed_22lat);

%% implementing the  barrier
kxyThresholdHigh = 10.0;
kxyThresholdLow = 0.1;
K_high = 0.005;
K_low = 0.005;

gammaMembraneApHigh = 3*K_high*heaviside(kxysquared - kxyThresholdHigh).*(kxysquared - kxyThresholdHigh).^2;
gammaMembraneLatHigh = 3*K_high*heaviside(1./sqrt(kxysquared) - kxyThresholdHigh).*(1./sqrt(kxysquared) - kxyThresholdHigh).^2;

gammaMembraneApLow =   -3*K_low*heaviside(kxyThresholdLow - kxysquared).*(kxysquared - kxyThresholdLow).^2;
gammaMembraneLatLow =  -3*K_low*heaviside(kxyThresholdLow - 1./sqrt(kxysquared)).*(1./sqrt(kxysquared) - kxyThresholdLow).^2;


gammaMembraneAp = gammaMembraneApLow + gammaMembraneApHigh;
gammaMembraneLat = gammaMembraneLatLow + gammaMembraneLatHigh;

%% Mass balance
Aap = Aap0.*kxysquared;
Alat = ah*V0./(sqrt(Aap));
dAlatdepsdot = -ah*V0./(2.*Aap.^1.5).*Aap0;
dAapdepsdot = Aap0;


%% Rayleighian
% Free energy

Fdotap = rhoAp.*Aap.*(2.*lambda.*Jap.*(Jap-1)./kxysquared);
Fdotlat = rhoLat.*Alat.*(mu./Jlat.*(Gxylat - 2.*Gz./kxysquared.^3 + Ilat/2./kxysquared) - lambda.*Jlat.*(Jlat - 1)./kxysquared);
ApicoBasalViscoElasticTension = 2*Fdotap/Aap0;
LateralViscoElasticTension = Fdotlat/Aap0;
TotalViscoElasticTension = ApicoBasalViscoElasticTension + LateralViscoElasticTension;

% active power (Eulerian)
ApicoBasalActiveTension = 2*rhoAp*xiap.*dAapdepsdot/Aap0;

LateralCorticalTension = rhoLat.*xilat;
LateralActiveTension = LateralCorticalTension.*dAlatdepsdot/Aap0;

TotalActiveTension = ApicoBasalActiveTension + LateralActiveTension;
Barrier = 2*gammaMembraneAp.*dAapdepsdot/Aap0 + gammaMembraneLat.*dAlatdepsdot/Aap0;

TotalTension = TotalActiveTension + TotalViscoElasticTension + Barrier;

% Gflat
Gxy = 1./Y(:,1);
Gxylat = 1./Y(:,2);
Gz = 1./Y(:,3);