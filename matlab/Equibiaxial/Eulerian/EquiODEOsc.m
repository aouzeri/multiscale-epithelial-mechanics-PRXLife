function [res] = EquiODEOsc(t,y,yp,params, stimulus)
%EQUIODE Function that encodes the system of DAE for uniaxial creep and
%stress relaxation
% Cs are all flat
% Gs are all sharp
%   y(1) = Gxy sharp apical
%   y(2) = Gxylat sharp lateral 
%   y(3) = Gz sharp lateral 
%   y(4) = strain 
%   y(5) = tauxy
%   y(6) = rho apical
%   y(7) = rho lateral

kp = params.kp; kd = params.kd; V0 = params.V0; xilat = params.xilat;
Aap0 = params.Aap0; mu = params.mu; lambda = params.lambda; eta = params.eta; ah = params.ah;
xiap = params.xiap; 

%% Solving
epsilon = y(4);
epsilondot = yp(4);

tauxy = stimulus(t); % uncomment for creep

% densities
rhoAp = y(6);
rhoLat = y(7);

%% for the remodelling dissipation
% deformation
Cxy = (epsilon + 1);
Cz = 1/(epsilon + 1)/(epsilon + 1);

% Apical
Gxy = y(1);
Cmixed_xy = Cxy*Gxy;
Iap = Cmixed_xy + Cmixed_xy;
Jap = sqrt(Cmixed_xy*Cmixed_xy);


% lateral parallel side
Gxylat = y(2);
Gz = y(3);

Cmixed_11lat = Cxy*Gxylat;
Cmixed_22lat = Cz*Gz;

Ilat = Cmixed_11lat + Cmixed_22lat;
Jlat = sqrt(Cmixed_11lat*Cmixed_22lat);

%% Calculating all Gdotsharp
dy1dt = -Gxy*lambda/eta*Jap*(Jap - 1);
dy2dt = -mu/eta*Gxylat/Jlat*(Gxylat*Cxy - Ilat/2) - lambda/eta*Jlat*(Jlat -1)*Gxylat;
dy3dt = -mu/eta*Gz/Jlat*(Gz*Cz - Ilat/2) - lambda/eta*Jlat*(Jlat -1)*Gz;

%% implementing the  barrier

J_thresholdHigh = 10;
J_thresholdLow = 0.1;

K_high = 0.005;
K_low = 0.005;

gammaMembraneApHigh = 3*K_high*heaviside((epsilon + 1) - J_thresholdHigh)*((epsilon + 1) - J_thresholdHigh)^2;
gammaMembraneLatHigh = 3*K_high*heaviside(1./sqrt(epsilon + 1) - J_thresholdHigh)*(1./sqrt(epsilon + 1) - J_thresholdHigh)^2;
 
gammaMembraneApLow =   -3*K_low*heaviside(J_thresholdLow - (epsilon + 1))*((epsilon + 1) - J_thresholdLow)^2;
gammaMembraneLatLow =  -3*K_low*heaviside(J_thresholdLow - 1./sqrt(epsilon + 1))*(1./sqrt(epsilon + 1) - J_thresholdLow)^2;

gammaMembraneAp = gammaMembraneApLow + gammaMembraneApHigh;
gammaMembraneLat = gammaMembraneLatLow + gammaMembraneLatHigh;

%% Mass balance

Aap = Aap0*(epsilon + 1);
Alat = ah*V0/(sqrt(Aap));
Aapdot = Aap0*epsilondot;
Alatdot = -ah*V0/(2*Aap^1.5)*epsilondot*Aap0;
dAlatdepsdot = -ah*V0/(2*Aap^1.5)*Aap0;
dAapdepsdot = Aap0;

rhoApdot= kp - rhoAp*(Aapdot/Aap + kd);
rhoLatdot = kp - rhoLat*(Alatdot/Alat + kd);

%% Rayleighian

Fdotap = rhoAp*Aap*(2*lambda*Jap*(Jap-1)/(epsilon + 1));
Fdotlat = rhoLat*Alat*(mu/Jlat*(Gxylat - 2*Gz/(epsilon + 1)^3 + Ilat/2/(epsilon + 1)) - lambda*Jlat*(Jlat - 1)/(epsilon + 1));
dFdotdepsdot = 2*Fdotap + Fdotlat;

ActivetermAp = 2*rhoAp*xiap*dAapdepsdot;
ActivetermLat = rhoLat*xilat*dAlatdepsdot;
Barrier = 2*gammaMembraneAp*dAapdepsdot + gammaMembraneLat*dAlatdepsdot;
dPadepsdot = ActivetermAp + ActivetermLat + Barrier;

% External power
dPextdepsdot = Aap0*tauxy;

%% System of equations

res1 = yp(1) - dy1dt;
res2 = yp(2) - dy2dt;
res3 = yp(3) - dy3dt;
res4 = dFdotdepsdot - dPextdepsdot + dPadepsdot;
res5 = yp(5);
res6 = yp(6) - rhoApdot;
res7 = yp(7) - rhoLatdot;

%% grouping everything together
res = [res1;res2;res3;res4;res5;res6;res7];

end

