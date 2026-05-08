%% IagVM
% Copyright (c) Adam Ouzeri 2025

tic
% clc
clear all;
close all;

%% 
load('../../PaperData/casares_data/casares_data_stretch.csv')
load('../../PaperData/casares_data/casares_data_unstretch.csv')
load('Y0_active_eshelby_0p55.mat')

%% normalising parameters
kd = 1;                 % depolymerization rate
AssumedTaut = 1;       % for comparison purposes
s0 = 1;                 % cell dimension 
xiap_active_Eshelby = 1.0; % apical active Eshelby parameter
%% kinetic parameters
kp = 0.03;           % polymerisation rate; 
rho0 = kp;                       % equilibrium cortical thickness
lambda = 0.4;  
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
xiLat_active_Eshelby = 0.55; % lateral active eshelby paramter

%% Initial conditions
y0 = Yend';

yp0 = zeros(size(y0));
tol = 1e-13;             % tolerance value for numerical solver
options = odeset('AbsTol', tol, 'RelTol', tol);

%% Mechanical sitmuli
teq = 0;
tfinal = 10;
tInterval = [0 tfinal]; 

%% creating parameter structure for easier implementatikxyon and play
params.V0 = V0; 
params.Aap0 = Aap0; params.kp = kp; params.kd = kd; params.rho0 = rho0;
params.eta = eta; params.Apar0 = Apar0; params.Atra0 = Atra0;
params.mu = mu; params.lambda = lambda; params.rho0 = rho0; params.s0 = s0;
params.Alat0 = Alat0; params.h0 = h0; params.Peri0 = Peri0;
params.ah = ah; params.xiLat_active_Eshelby = xiLat_active_Eshelby; 
params.xiap_active_Eshelby = xiap_active_Eshelby;
params.tfinal = tfinal; params.teq = teq; params;


%% Solving
sol = ode89(@(t,g) EquiODESR(t,g,params), tInterval ,y0, options);
T = sol.x';
Y = sol.y';

%% Post-processing
nTimeSteps = length(Y);
Gxy = 1./Y(:,1);
Gxylat = 1./Y(:,2);
Gz = 1./Y(:,3);
rhoAp = Y(:,4);
rhoLat = Y(:,5);
strain = Y(:,6);

ArealDeformation = strain + 1;
Aap = Aap0 * ArealDeformation;
Alat = Alat0./sqrt(ArealDeformation);
h = V0./(Aap0*ArealDeformation);


%% computing the tension
calculateTensions;

%% Plotting as a function of time
Time = T;
NormalizingTension = TotalTension(length(T));

figure(1)
hold on
plot(Time, TotalTension/NormalizingTension, 'DisplayName', 'tension', 'LineWidth', 2)
plot(casares_data_stretch(:,1), casares_data_stretch(:,2),'d','markerfacecolor','b', 'Markers', 20)
plot(casares_data_unstretch(:,1)-4.5, casares_data_unstretch(:,2),'d','markerfacecolor','b', 'Markers', 20)
title("Tension vs time")

toc

