%% Idealised VertexModel
%% Copyright (c) Adam Ouzeri 2026

tic
% clc
clear all;
close all;

%% 
load('Y0_active_eshelby_0p5_5p4801_strain.mat')

%% normalising parameters
kd = 1;                 % depolymerization rate
s0 = 1;                 % cell dimension
xiap_active_Eshelby = 1.0; % apical active Eshelby parameter

%% kinetic parameters
kp = 0.03;           % polymerisation rate; 
rho0 = kp;           % equilibrium cortical thickness
lambda = 1.0;
mu = lambda;         % shear modulus
eta = 0.3;           % viscous remodeling; 

%% Geometric parameters
h0 = 2;
Aap0 = 3*sqrt(3)/2*s0^2;           
V0 = h0*Aap0;
Peri0 = 6*s0;
Apar0 = s0*h0;
Atra0 = s0*h0;
Alat0 = 2*Apar0 + 4*Atra0;
Acell_0 = 2*Aap0 + Alat0;			
ah = sqrt(72/3/sqrt(3));
K = ah*h0/sqrt(Aap0)/2;

%% Tensional parameters
xiLat_active_Eshelby = 0.5; 

%% initial conditions for lagrangian activity with lambda = 1.0, eta = 0.3,
CreepTension = 1.93*rho0;
y0 = [Yend0p5_5p4801_strain(1:3)'; 5.4856; CreepTension; rho0; rho0]; % strain 5.4801 (perturbed = 5.4856) <-> tensio 1.93

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
params.V0 = V0;
params.Aap0 = Aap0; params.kp = kp; params.kd = kd; params.rho0 = rho0;
params.eta = eta; params.Apar0 = Apar0; params.Atra0 = Atra0;
params.mu = mu; params.lambda = lambda; params.rho0 = rho0; params.s0 = s0;
params.Alat0 = Alat0; params.h0 = h0; params.Peri0 = Peri0;
params.ah = ah;
params.tfinal = tfinal; params.teq = teq; params;
params.xiLat_active_Eshelby = xiLat_active_Eshelby; 
params.xiap_active_Eshelby = xiap_active_Eshelby;

%% Solving the DAE
[T,Y] = ode15i(@(t,g,gp) EquiODEOsc(t,g,gp,params, stimulus), tInterval ,y0, yp0, options);

%% Post-processing
n1 = length(Y);
Gxy = 1./Y(:,1);
Gxylat = 1./Y(:,2);
Gz = 1./Y(:,3);
strain = Y(:,4);
TotalTension = Y(:,5);
rhoAp = Y(:,6); 
rhoLat = Y(:,7);

ArealDeformation = strain + 1;
Aap = Aap0 * ArealDeformation;
Alat = Alat0./sqrt(ArealDeformation);
h = V0./(Aap0*ArealDeformation);

%% Plots
Time = T;
NormalisingTau = rho0;

figure(1)
plot(rhoAp/rho0,strain,'linewidth',2)
xlabel('apical density')
ylabel('Strain')
