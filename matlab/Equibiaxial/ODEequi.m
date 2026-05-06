function [dydt] = EquiODE45(t,y,params)
%EQUIODE45 Function that encodes the system of DAE for uniaxial creep and
%stress relaxation
% Cs are all flat
% Gs are all sharp
%   y(1) = Gxy sharp apical
%   y(2) = Gxylat sharp lateral 
%   y(3) = Gz sharp lateral  
%   y(4) = rho apical
%   y(5) = rho lateral

kp = params.kp; kd = params.kd; V0 = params.V0;
Aap0 = params.Aap0; mu = params.mu; lambda = params.lambda; eta = params.eta; ah = params.ah;
xiLat_active_Eshelby = params.xiLat_active_Eshelby; xiap_active_Eshelby = params.xiap_active_Eshelby;

%% Ramping the strain
% if t < 0.02
%     stimulus = 5; 
% elseif t > 5 && t < 5.02
%      stimulus = -5;
% else
%     stimulus = 0;
% end

%% Ramping the strain
if t > -1.0
    stimulus = 0.0001; 
end

if t > 80000
     stimulus = -stimulus;
end

if  y(6) < 5.4801
     stimulus = 0.0001; 
else
    stimulus = 0.0;
end


%% Solving
epsilon = y(6);
epsilondot = stimulus;

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

% densities
rhoAp = y(4);
rhoLat = y(5);

%% calculating all Gdotsharp
dG1dt = -Gxy*lambda/eta*Jap*(Jap - 1);
dG2dt = -mu/eta*Gxylat/Jlat*(Gxylat*Cxy - Ilat/2) - lambda/eta*Jlat*(Jlat -1)*Gxylat;
dG3dt = -mu/eta*Gz/Jlat*(Gz*Cz - Ilat/2) - lambda/eta*Jlat*(Jlat -1)*Gz;

% Adding the active Eshelby stress contribution
Active_Eshelby_ap = 1/2/eta * xiap_active_Eshelby * Gxy ;
Active_Eshelby_lat = 1/2/eta * xiLat_active_Eshelby * Gxylat ;
Active_Eshelby_z = 1/2/eta * xiLat_active_Eshelby * Gz ;


%% Mass balance
Aap = Aap0*(epsilon + 1);
Alat = ah*V0/(sqrt(Aap));
Aapdot = Aap0*epsilondot;
Alatdot = -ah*V0/(2*Aap^1.5)*epsilondot*Aap0;

rhoApdot= kp - rhoAp*(Aapdot/Aap + kd);
rhoLatdot = kp - rhoLat*(Alatdot/Alat + kd);

dydt = [dG1dt+ Active_Eshelby_ap; dG2dt + Active_Eshelby_lat; dG3dt  + Active_Eshelby_z; rhoApdot; rhoLatdot; epsilondot];

end

