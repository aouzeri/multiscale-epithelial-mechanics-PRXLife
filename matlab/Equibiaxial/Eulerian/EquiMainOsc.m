%% IagVM
%% Copyright (c) Adam Ouzeri 2025

tic
% clc
clear all;
close all; 

%% normalising parameters
kd = 1;                 % depolymerization rate
AssumedTaut = 1;       % for comparison purposes
s0 = 1;                 % cell dimension
xiap = 1.0;               % apical stress in reference configuraion

%% kinetic parameters
kp = 0.03;           % polymerisation rate; 
rho0 = kp;                       % equilibrium cortical thickness
lambda = 1.0;
mu = lambda;               % shear modulus 
eta = 0.3;     % viscous remodeling;

%% Geometric parameters
h0 = 2;
Aap0 = 3*sqrt(3)/2*s0^2;            % reference cell area, [micro.m^2]
V0 = h0*Aap0;
Peri0 = 6*s0;
Apar0 = s0*h0;
Atra0 = s0*h0;
Alat0 = 2*Apar0 + 4*Atra0;
Acell_0 = 2*Aap0 + Alat0;			% reference total cell area, [micro.m^2]
ah = sqrt(72/3/sqrt(3));
K = ah*h0/sqrt(Aap0)/2;

%% Tensional parameters
xilat = 0.5;

%% Initial conditions
CreepTension = 1.93*rho0; % 1.8; 1.93; 1.96

%% Calculating the fixed point for initial conditions (tissue pulsation)
options = optimoptions('fsolve','Display','iter', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 10000, ...
    'OptimalityTolerance', 1e-12);
FP = fsolve(@(x) myDAEfun(x, CreepTension, eta, lambda, mu, kp, kd, s0, xiap, xilat, h0, ...
    rho0, Aap0, V0, Peri0,Apar0, Atra0 , Alat0, Acell_0, ah),[    0.9086    0.9086    1.2114    0.0858    0.3431 0.1006], options); 
y0 = [FP(1); FP(2); FP(3); FP(6); CreepTension; FP(4) /(Aap0*(FP(6) +1)); FP(5)/(V0*ah/sqrt(Aap0*(FP(6) + 1)))];  
y0(4) = y0(4)*1.001; % perturbing strain to see oscillation around fixed point
params.SSstrain = FP(6);

yp0 = zeros(size(y0));
tol = 2e-6; 
options = odeset('AbsTol', tol, 'RelTol', tol);

%% Mechanical sitmuli
teq = 0;
tfinal = 100;
tStep = 0.0001;
tInterval = [0 tfinal]; 
stimulus = @(x) CreepTension;

%% creating parameter structure for easier implementatikxyon and play
params.xiap = xiap;  params.V0 = V0; params.xilat = xilat;
params.Aap0 = Aap0; params.kp = kp; params.kd = kd; params.rho0 = rho0;
params.eta = eta; params.Apar0 = Apar0; params.Atra0 = Atra0;
params.mu = mu; params.lambda = lambda; params.rho0 = rho0; params.s0 = s0;
params.Alat0 = Alat0; params.h0 = h0; params.Peri0 = Peri0;
params.ah = ah;

%% Solving the DAE
[T,Y] = ode15i(@(t,g,gp) EquiODEOsc(t,g,gp,params, stimulus), tInterval ,y0, yp0, options);

%% Post-processing
n1 = length(Y);
Gxy = 1./Y(:,1);
Gxylat = 1./Y(:,2);
Gz = 1./Y(:,3);
strain = Y(:,4);
TotalTension = Y(:,5);
rhoAp = Y(:,6); % apical and basal density
rhoLat = Y(:,7); % lateral density

ArealDeformation = strain + 1;
Aap = Aap0 * ArealDeformation;
Alat = Alat0./sqrt(ArealDeformation);
h = V0./(Aap0*ArealDeformation);

%% Plotting as a function of time
Time = T;
NormalisingTau = rho0; 

figure(1)
plot(rhoAp/rho0,strain,'linewidth',2)
xlabel('apical density')
ylabel('Strain')

toc

function Fun = myDAEfun( x, CreepTension, eta, lambda, mu, kp, kd, s0, xiap, xilat, h0, ...
    rho0, Aap0, V0, Peri0,Apar0, Atra0 , Alat0, Acell_0, ah)

%% calculating the number of fixed points
sigma0 = CreepTension;
taurmu = eta/mu;
taurlambda = eta/lambda;

%% Calculating the matrix
G1ap = x(1);
G1lat = x(2);
G2lat = x(3);
Nap = x(4);
Nlat = x(5);
epsilon = x(6);
G2ap = G1ap;

Gap = [G1ap, 0; 0, G2ap];
Cap = [epsilon + 1, 0 ; 0 , epsilon + 1];

Glat = [G1lat, 0; 0, G2lat];
Clat = [epsilon + 1, 0 ; 0 , 1/(epsilon  + 1)^2];

Iap = sum(sum(Cap.*Gap));
Jap = sqrt(det(Cap)*det(Gap));

Ilat = sum(sum(Clat.*Glat));
Jlat = sqrt(det(Clat)*det(Glat));

Fun(1) = -G1ap*G1ap/taurlambda*(Jap*(Jap - 1)/G1ap);
Fun(2) = -G1lat*G1lat*((epsilon + 1/Jlat - Ilat/(2*Jlat*G1lat))/taurmu + Jlat*(Jlat-1)/(taurlambda*G1lat));
Fun(3) = -G2lat*G2lat*((1/(epsilon  + 1)^2/Jlat - Ilat/(2*Jlat*G2lat))/taurmu + Jlat*(Jlat-1)/(taurlambda*G2lat));
Fun(4) = kp*Aap0*(epsilon + 1) - kd*Nap;
Fun(5) = kp*V0*ah/sqrt(Aap0*(epsilon + 1)) - kd*Nlat;

Fun(6) = (Nap/Aap0*( 2 + (2*lambda*(Jap -1))*Jap ) - Nlat/Aap0*( xilat/2 + (-mu*Ilat/Jlat^2 + 2*lambda*(Jlat - 1))*Jlat/2 ) )/(epsilon + 1) ...
    - (Nlat/Aap0*mu/Jlat*2*G2lat)/(epsilon + 1)^3 ...
    + Nlat/Aap0*mu/Jlat*G1lat - sigma0;

Fun(1) = Fun(1);
Fun(2) = Fun(2);
Fun(3) = Fun(3);

end

