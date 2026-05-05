%% Copyright (c) Adam Ouzeri 2025

function [res] = ODEuniaxial(t,y,yp,params)
%ODEUNIAXIAL Function that encodes the system of DAE for uniaxial creep and
%stress relaxation
%   y(1) = G11 sharp apical
%   y(2) = G22 sharp apical
%   y(3) = G11 sharp lateral parallel
%   y(4) = G22 sharp lateral parallel
%   y(5) = G11 sharp lateral transverse
%   y(6) = G22 sharp lateral transverse
%   y(7) = kappa
%   y(8) = tau
%   y(9) = rho apical
%   y(10) = rho lateral parallel
%   y(11) = rho lateral transverse

res = zeros(11,1);

kp = params.kp; kd = params.kd; Aap0 = params.Aap0; A_lat_p0 = params.A_lat_p0; A_lat_t0 = params.A_lat_t0;
mu = params.mu; lambda = params.lambda; eta = params.eta; C0 = params.C0;
xiap = params.xiap; xilat = params.xilat; appliedTension = params.appliedTension;

% deformation coefficients
kappa = y(7);
kappadot = yp(7);

% applied tension
tau = appliedTension*(1 - exp(-t/0.0001)); 

% densities
rho_ap = y(9);
rho_lat_p = y(10);
rho_lat_t = y(11);

gamma_ap = rho_ap*xiap;
gamma_lat_p = rho_lat_p*xilat;
gamma_lat_t = rho_lat_t*xilat;

%% metric components
Gsharp_ap = [y(1), 0; 0, y(2)];
Gsharp_lat_p = [y(3), 0; 0, y(4)];
Gsharp_lat_t = [y(5), 0; 0, y(6)];

Gflat_ap = Gsharp_ap\eye(2);
Gflat_lat_p = Gsharp_lat_p\eye(2);
Gflat_lat_t = Gsharp_lat_t\eye(2);

Gsharpdot_ap = [yp(1), 0; 0, yp(2)];
Gsharpdot_lat_p = [yp(3), 0; 0, yp(4)];
Gsharpdot_lat_t = [yp(5), 0; 0, yp(6)];

%% kinematics
J_ap = kappa;
J_lat_p = 1;
J_lat_t = sqrt(kappa^2 + 3)/2/kappa;

Jdot_ap = kappadot;
Jdot_lat_p = 0;
Jdot_lat_t = -3*kappadot/(2*kappa^2*sqrt(kappa^2 + 3));

Cflat_ap = [kappa^2, 0; 0, 1];
Cflat_lat_p = [kappa^2, 0; 0, 1/kappa^2];
Cflat_lat_t = [(kappa^2 + 3)/4, 0; 0, 1/kappa^2];

C_inverse_sharp_ap = Cflat_ap\eye(2);
C_inverse_sharp_lat_p = Cflat_lat_p\eye(2);
C_inverse_sharp_lat_t = Cflat_lat_t\eye(2);

Cmixed_ap = Cflat_ap*Gsharp_ap;
Cmixed_lat_p  = Cflat_lat_p *Gsharp_lat_p;
Cmixed_lat_t  = Cflat_lat_t *Gsharp_lat_t;

I1_ap = trace(Cmixed_ap);
I3_ap = sqrt(det(Cmixed_ap));
I1_lat_p = trace(Cmixed_lat_p);
I3_lat_p = sqrt(det(Cmixed_lat_p));
I1_lat_t = trace(Cmixed_lat_t);
I3_lat_t = sqrt(det(Cmixed_lat_t));

dCflat_ap_dkappa = [2*kappa,0;0,0];
dCflat_lat_p_dkappa = [2*kappa,0;0,-2/kappa^3];
dCflat_lat_t_dkappa = [kappa/2,0;0,-2/kappa^3];

%% variation of free energy with respect to Gsharpdot
dPsi_ap_dI1_ap = mu/I3_ap;
dPsi_ap_dI3_ap = 2*lambda*(I3_ap - 1) - mu*I1_ap/I3_ap^2;
dPsi_lat_p_dI1_lat_p = mu/I3_lat_p;
dPsi_lat_p_dI3_lat_p = 2*lambda*(I3_lat_p - 1) - mu*I1_lat_p/I3_lat_p^2;
dPsi_lat_t_dI1_lat_t = mu/I3_lat_t;
dPsi_lat_t_dI3_lat_t = 2*lambda*(I3_lat_t - 1) - mu*I1_lat_t/I3_lat_t^2;

dI1_ap_dGsharp_ap = Cflat_ap;
dI1_lat_p_dGsharp_lat_p = Cflat_lat_p;
dI1_lat_t_dGsharp_lat_t = Cflat_lat_t;

dI3_ap_dGsharp_ap = I3_ap*Gflat_ap/2;
dI3_lat_p_dGsharp_lat_p = I3_lat_p*Gflat_lat_p/2;
dI3_lat_t_dGsharp_lat_t = I3_lat_t*Gflat_lat_t/2;

dPsi_ap_dGsharp_ap = dPsi_ap_dI1_ap*dI1_ap_dGsharp_ap + dPsi_ap_dI3_ap*dI3_ap_dGsharp_ap;
dPsi_lat_p_dGsharp_lat_p = dPsi_lat_p_dI1_lat_p*dI1_lat_p_dGsharp_lat_p + dPsi_lat_p_dI3_lat_p*dI3_lat_p_dGsharp_lat_p;
dPsi_lat_t_dGsharp_lat_t = dPsi_lat_t_dI1_lat_t*dI1_lat_t_dGsharp_lat_t + dPsi_lat_t_dI3_lat_t*dI3_lat_t_dGsharp_lat_t;

%% variation of dissipation potential with respect to Gsharpdot
dD_ap_dGsharpdot_ap = eta*Gflat_ap*Gsharpdot_ap*Gflat_ap;
dD_lat_p_dGsharpdot_lat_p = eta*Gflat_lat_p*Gsharpdot_lat_p*Gflat_lat_p;
dD_lat_t_dGsharpdot_lat_t = eta*Gflat_lat_t*Gsharpdot_lat_t*Gflat_lat_t;

%% minimisaiton with respect to Gsharpdot
res(1:2) = diag(dPsi_ap_dGsharp_ap + dD_ap_dGsharpdot_ap);

res(3:4) = diag(dPsi_lat_p_dGsharp_lat_p + dD_lat_p_dGsharpdot_lat_p);

res(5:6) = diag(dPsi_lat_t_dGsharp_lat_t + dD_lat_t_dGsharpdot_lat_t);

%% variation of free energy with respect to kappadot
dI1_ap_dCflat_ap = Gsharp_ap;
dI1_lat_p_dCflat_lat_p = Gsharp_lat_p;
dI1_lat_t_dCflat_lat_t = Gsharp_lat_t;

dI3_ap_dCflat_ap = I3_ap*C_inverse_sharp_ap/2;
dI3_lat_p_dCflat_lat_p = I3_lat_p*C_inverse_sharp_lat_p/2;
dI3_lat_t_dCflat_lat_t = I3_lat_t*C_inverse_sharp_lat_t/2;

dPsi_ap_dCflat_ap = dPsi_ap_dI1_ap*dI1_ap_dCflat_ap + dPsi_ap_dI3_ap*dI3_ap_dCflat_ap;
dPsi_lat_p_dCflat_lat_p = dPsi_lat_p_dI1_lat_p*dI1_lat_p_dCflat_lat_p + dPsi_lat_p_dI3_lat_p*dI3_lat_p_dCflat_lat_p;
dPsi_lat_t_dCflat_lat_t = dPsi_lat_t_dI1_lat_t*dI1_lat_t_dCflat_lat_t + dPsi_lat_t_dI3_lat_t*dI3_lat_t_dCflat_lat_t;

E_ap = rho_ap*trace(dPsi_ap_dCflat_ap'*dCflat_ap_dkappa)*J_ap*Aap0;
E_lat_p = rho_lat_p*trace(dPsi_lat_p_dCflat_lat_p'*dCflat_lat_p_dkappa)*J_lat_p*A_lat_p0;
E_lat_t = rho_lat_t*trace(dPsi_lat_t_dCflat_lat_t'*dCflat_lat_t_dkappa)*J_lat_t*A_lat_t0;

dPsi_dkappadot = 2*E_ap + E_lat_p + E_lat_t;

%% variation of active power input with respect to kappadot
dJdot_ap_dkappadot = 1;
dJdot_lat_p_dkappadot = 0;
dJdot_lat_t_dkappadot = -3/(2*kappa^2*sqrt(kappa^2 + 3));

dPa_ap_dkappadot = gamma_ap*dJdot_ap_dkappadot*Aap0;
dPa_lat_p_dkappadot = gamma_lat_p*dJdot_lat_p_dkappadot*A_lat_p0;
dPa_lat_t_dkappadot = gamma_lat_t*dJdot_lat_t_dkappadot*A_lat_t0;

dPatot_dkappadot = 2*dPa_ap_dkappadot + dPa_lat_p_dkappadot + dPa_lat_t_dkappadot;

%% variation of external power input with respect to kappadot
dPextdkappadot = Aap0*tau;
 

%% minimisation with respect to kappadot
res(7) = dPsi_dkappadot + dPatot_dkappadot - dPextdkappadot;

%% Mass balance

rhodot_ap = kp*C0 - rho_ap*(Jdot_ap/J_ap + kd);
rhodot_lat_p = kp*C0 - rho_lat_p*(Jdot_lat_p/J_lat_p + kd);
rhodot_lat_t = kp*C0 - rho_lat_t*(Jdot_lat_t/J_lat_t + kd);

res(8:10) = [yp(9) - rhodot_ap; yp(10) - rhodot_lat_p; yp(11) - rhodot_lat_t];

%% Applied stiumulus
res(11) = yp(8);

%% fractional element approximation
dDdkappadot = (params.passedF)*2*params.materialConstant/gamma(1-params.alphaPower)/kappa*params.deltat*params.V0;

res(7) = res(7) + dDdkappadot;

end