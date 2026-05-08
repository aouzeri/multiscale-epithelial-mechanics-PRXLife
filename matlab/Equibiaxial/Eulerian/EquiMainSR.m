%% IagVM
% Copyright (c) Adam Ouzeri 2025

tic
% clc
clear all;
close all;

%% 
load('../../PaperData/casares_data/casares_data_stretch.csv')
load('../../PaperData/casares_data/casares_data_unstretch.csv')

%% normalising parameters
kd = 1;                 % depolymerization rate
AssumedTaut = 1;       % for comparison purposes
s0 = 1;                 % cell dimension
xiap = 1.0;                % apical activity
%% kinetic parameters
kp = 0.03;           % polymerisation rate; 
rho0 = kp;                       % equilibrium cortical thickness
lambda = 0.4; 
mu = lambda;               % shear modulus
eta = 0.3;     % viscous remodeling; 

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
xilat = 0.55;

%% Initial conditions
epsilon0 = 0;
G0 = 1;
y0 = [G0; G0; G0; rho0; rho0; epsilon0];

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
params.ah = ah; 
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
rhoAp = Y(:,4); % apical and basal density
rhoLat = Y(:,5); % lateral density
strain = Y(:,6);

ArealDeformation = strain + 1;
Aap = Aap0 * ArealDeformation;
Alat = Alat0./sqrt(ArealDeformation);
h = V0./(Aap0*ArealDeformation);

%% computing the tension
calculateTensions;

%% Plotting as a function of time
Time = T;
% T = Time*AssumedTaut/60;

NormalizingTension = TotalTension(length(T));

figure(2)
hold on
plot(Time, TotalActiveTension/NormalizingTension, 'DisplayName', 'Active tension', 'LineWidth', 2)
plot(Time, TotalViscoElasticTension/NormalizingTension, 'DisplayName', 'VE tension', 'LineWidth', 2)
plot(Time, TotalTension/NormalizingTension, 'DisplayName', 'tension', 'LineWidth', 2)
plot(casares_data_stretch(:,1), casares_data_stretch(:,2),'d','markerfacecolor','b', 'Markers', 20)
plot(casares_data_unstretch(:,1)-4.5, casares_data_unstretch(:,2),'d','markerfacecolor','b', 'Markers', 20)
title("Tension vs time")

figure(4)
hold on
plot(Time, strain + 1, 'DisplayName', 'Strain + 1', 'LineWidth', 2)
plot(Time, Gxy, 'DisplayName', 'Gxy', 'LineWidth', 2)
plot(Time, Gxylat, 'DisplayName', 'Gxylat', 'LineWidth', 2)
plot(Time, Gz, 'DisplayName', 'Gz', 'LineWidth', 2)
title("Strain vs time")

figure(3)
hold on
plot(Time, rhoAp/rho0, 'DisplayName', 'rhoAp', 'LineWidth', 2)
plot(Time, rhoLat/rho0, 'DisplayName', 'rhoLat', 'LineWidth', 2)
title("Density vs strain")



%% Quasistatic
strain = -0.95:0.001:10;
s0 = 1;
Aap0 = 3*sqrt(3)/2*s0^2;            % reference cell area
Aap = Aap0 * (strain+1);
h0 = 2;
V0 = h0*Aap0;
xilat = 0.44;
ah = sqrt(72/3/sqrt(3));
K = ah*h0/sqrt(Aap0)/2;
tauQS = (2 - xilat*K./sqrt((strain + 1).^3));
% figure(1)
% hold on
% plot(strain,(tauQS - tauQS(1))*0.5)
% plot(strain,tauQS, 'DisplayName', 'Quasistatic analytical', 'LineWidth', 2)

toc

