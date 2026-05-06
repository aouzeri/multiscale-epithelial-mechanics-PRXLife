%% Stress relaxation under uniaxial deformaion
%% Copyright (c) Adam Ouzeri 2025

tic
clc
clear all;
close all;

%% normalising parameters
kd = 1;                 % depolymerization rate [1/s]
s0 = 1;                 % cell dimension [micro.m]
xiap = 1;               % apical stress in reference configuraion, [N/micro.m2]
AssumedTaut = 30;       % for comparison purposes

%% independent parameters
kp = 0.03*kd/s0;           % polymerisation rate; kp/kd = rho0 = 0.03; need to keep it two orders of magnitude smaller than cell dimension
rho0 = kp;
C0 = 1;
lambda = 0.6;           % bulk modulus; mu/l3mbda = 1;
mu = lambda;               % shear modulus
eta = 0.2;     % viscous remodeling; need to keed remodelling timescale one order of magnitude smaller to kd

%% dependent parameters
xilat = 0.4;
h0 = 2;
% h0 = sqrt(3)/xiLatCortex;
Aap0 = 3*sqrt(3)/2*s0^2;            % reference cell area, [micro.m^2]
V0 = h0*Aap0;
Peri0 = 6*s0;
A_lat_p0 = 2*s0*h0;
A_lat_t0 = 4*s0*h0;			% reference total cell area, [micro.m^2]
ah = sqrt(72/3/sqrt(3));
K = ah*h0/sqrt(Aap0)/2;
V0 = h0*Aap0;

%% Mechanical sitmuli
InitialTension =  2*rho0*xiap - K*rho0*xilat;
params.InitialTension = InitialTension;

% SR
params.deltat = 5e-3/AssumedTaut;
cyclePeriod = 20/AssumedTaut;
tfinal = 120/AssumedTaut;

%% Fractional elements
params.alphaPower = 0.4; 
params.materialConstant = 0.1;

%% Initial conditions
tau0 = 0;
kappa0 = 1;
G0 = 1;
rho0 = kp/kd;                       % reference thickness
y0 = [G0; G0; G0; G0; G0; G0; kappa0; tau0; rho0; rho0; rho0];
yp0 = zeros(size(y0));

%%
tol = 1e-08;             % tolerance value for numerical solver
options = odeset('AbsTol', tol, 'RelTol', tol);

% creating parameter structure for easier implementation
params.xiap = xiap; params.Aap0 = Aap0; params.kp = kp; params.kd = kd; params.rho0 = rho0;
params.eta = eta; params.A_lat_p0 = A_lat_p0; params.A_lat_t0 = A_lat_t0;
params.mu = mu; params.lambda = lambda; params.rho0 = rho0; params.s0 = s0;
params.xilat = xilat; params.C0 = C0; params.V0 = V0;

%% Simulation
% keyboard;

tn = 0;
tnp1 = tn + params.deltat;
time = tn;

kappa = kappa0;
kappadot = 0;
Tau = tau0;
rho_ap = rho0;
rho_lat_p = rho0;
rho_lat_t = rho0;

MetricTensorSharp = ones(1,6);
nIter = 0;
n = 0;
tstart = 0;
passedF = [];

while(tn < tfinal)
            
    % fractional element approximation
    params.passedF = (tnp1 - time').^(-params.alphaPower) * (kappadot./kappa);    

    % solving
    params.tnp1 = tnp1;
    sol = ode15i(@(t,g,gp) ODEuniaxialSR(t,g,gp,params), [tn tnp1],y0,yp0);%,options);
    T = sol.x;    
    
    if length(T) > 2
                
        [y0,yp0] = deval(sol,sol.x(end)); % updating y0 and yp0 with last

        MetricTensorSharp = [MetricTensorSharp;y0(1:6)'];
        kappa = [kappa; y0(7)];
        Tau = [Tau;y0(8)];
        %     Taux = [Taux;params.appliedTension];
        rho_ap = [rho_ap; y0(9)];
        rho_lat_p = [rho_lat_p; y0(10)];
        rho_lat_t = [rho_lat_t; y0(11)];
        
        kappadot = [kappadot;yp0(7)];
        
        tn = tnp1; 
        tnp1 = tnp1 + params.deltat;  
        time = [time;tn];
        disp(['Reached time : ', num2str(tn*AssumedTaut)]);
        tol = 1e-08;  % tol = 1e-05;            % tolerance value for numerical solver
        options = odeset('AbsTol', tol, 'RelTol', tol);
        nIter = nIter +1;
        
    else
        tnp1 = tnp1 + params.deltat; 
        [tn tnp1]
    end

end
    
% Retrieving the different tensions
calculateTensions;

%% Plotting

% Data from simulations
experimentalDataRaw = load('../PaperData/khalilgharibi_data/Khalilgharibi2019_FigS3.csv');
% 

EstimatedHeight = 1e-5; % from Khaligharibi2019_sup h = ~10 um
experimentalData = experimentalDataRaw(:,2)*EstimatedHeight * 1000; % 1000 mN/N
experimentalTime = experimentalDataRaw(:,1) - experimentalDataRaw(1,1);

% Data from experiments
B = min(experimentalData);
A = 1100 - B;
paperSimulation = A*experimentalTime.^(-0.3).*exp(-experimentalTime/14.9) + B;


% close all;
T = time*AssumedTaut;
% T = T - 0.6;
DimensionalTension = 5e-6 * 6.5e3 * 1000; % s0 * xi_ap [m * N/m^2] * 1000 [mN/n]
Tension = Tau;
stretch = kappa;

% Setting the 0 at the end of the step
startIdx = find(kappa >= 1.3, 1) - 2;
T = (time - time(startIdx))*AssumedTaut;

alphasign = 1;
NormalisingTau = rho0*xiap; % apical tension

figure(1)
hold on
loglog(T, stretch, 'linewidth',2, 'DisplayName', 'k^2 Uniaxial')
xlabel('Time [s]'); ylabel('Strain [-]');

idxEnd = find(T > 60, 1);
NoCytoplasmicTension = ActiveTension + ViscoelasticTension;
NoViscoElasticTension = ActiveTension + SpringpotTension;

figure(2)
hold on

loglog(experimentalTime, experimentalData,'d','markerfacecolor','b', 'Markers', 20)
loglog(T, (Tension)*DimensionalTension, 'linewidth', 4)
loglog(T, (Tension - Tension(idxEnd))*DimensionalTension, 'linewidth', 4)
plot(T, SpringpotTension*DimensionalTension, 'linewidth', 4)
plot(T, ViscoelasticTension*DimensionalTension, 'linewidth', 4)
plot(T, ActiveTension*DimensionalTension, 'linewidth', 4)
plot(T, NoCytoplasmicTension*DimensionalTension, 'linewidth', 4)
plot(T, NoViscoElasticTension*DimensionalTension, 'linewidth', 4)


interval1 = T > 0.1 & T< 2;
interval2 = T > 10 & T< 45;
funLog = (Tension(interval1) - Tension(idxEnd))*DimensionalTension; 
funExp = (Tension(interval2) - Tension(idxEnd))*DimensionalTension;
fitCoefflog = polyfit(log(T(interval1)),log(funLog),1);
fitLog = T.^fitCoefflog(1).*exp(fitCoefflog(2));
fitCoeffexp = polyfit(T(interval2), log(funExp),1);
fitExp = exp(T*fitCoeffexp(1))*exp(fitCoeffexp(2));
disp(['Power law coefficient : ', num2str(fitCoefflog(1))]);
disp(['Exponential timescale : ', num2str(-1/fitCoeffexp(1))]);
loglog(T, fitLog, 'k--', 'linewidth', 4)
loglog(T, fitExp, 'k--', 'linewidth', 4)

legend('Experiment', 'Total tension', 'Total tension  - Total tension at 60s', ...
    'Cytoplasmic tension', 'Viscoelastic tension', 'Active tension', ...
    'No cytoplasmic tension', 'No viscoelastic tension', ...
    'Power law fit', 'Exponential fit')

xlabel('Time [s]'); ylabel('Tension [mN/m]');


figure(3)
loglog(experimentalTime, experimentalData,'d','markerfacecolor','b', 'Markers', 20)
hold on
loglog(T, (Tension)*DimensionalTension, 'linewidth', 4)
xlabel('Time [s]'); ylabel('Tension [mN/m]');
xlim([0.1 120])

figure(4)
loglog(T, (Tension - Tension(idxEnd))*DimensionalTension, 'linewidth', 4)
hold on
loglog(T, fitLog, 'k--', 'linewidth', 4)
loglog(T, fitExp, 'k--', 'linewidth', 4)
xlabel('Time [s]'); ylabel('Tension [mN/m]');
xlim([0.1 55])

figure(5)
semilogy(T, (Tension - Tension(idxEnd))*DimensionalTension, 'linewidth', 4)
hold on
semilogy(T, (NoCytoplasmicTension - NoCytoplasmicTension(idxEnd))*DimensionalTension, 'linewidth', 4)
xlabel('Time [s]'); ylabel('Tension [mN/m]');


toc
