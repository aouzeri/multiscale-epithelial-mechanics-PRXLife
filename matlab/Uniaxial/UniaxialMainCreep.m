%% Creep under uniaxial tension
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

% creep
params.deltat = 0.5e-1/AssumedTaut;
cyclePeriod = 20/AssumedTaut;
tfinal = 10000/AssumedTaut;

%% Fractional elements
params.alphaPower = 0.4;
params.materialConstant = 0.1;

%% Initial conditions
tau0 = InitialTension;
kappa0 = 1;
G0 = 1;
rho0 = kp/kd;                       % reference thickness
y0 = [G0; G0; G0; G0; G0; G0; kappa0; tau0; rho0; rho0; rho0];
yp0 = zeros(size(y0));

%%
tol = 1e-10;             % tolerance value for numerical solver
options = odeset('AbsTol', tol, 'RelTol', tol);

% creating parameter structure for easier implementation
params.xiap = xiap; params.Aap0 = Aap0; params.kp = kp; params.kd = kd; params.rho0 = rho0;
params.eta = eta; params.A_lat_p0 = A_lat_p0; params.A_lat_t0 = A_lat_t0;
params.mu = mu; params.lambda = lambda; params.rho0 = rho0; params.s0 = s0;
params.xilat = xilat; params.C0 = C0; params.V0 = V0;

%% Simulation

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
        
    % relaxing the tissue to the desired starting tension
    if tn < 500/AssumedTaut
        params.appliedTension =  0.7*InitialTension;
        params.materialConstant = 0;
    else
        params.appliedTension =  2.3*0.03;
%         params.appliedTension =  1.3*InitialTension;
        params.materialConstant = 0.1;
    end
    
    % fractional element approximation
    params.passedF = (tnp1 - time').^(-params.alphaPower) * (kappadot./kappa); 
    params.tnp1 = tnp1;

    % solving
    sol = ode15i(@(t,g,gp) ODEuniaxialCreep(t,g,gp,params), [tn tnp1],y0,yp0);%,options);
    T = sol.x;    
    
    if length(T) > 2
                
        [y0,yp0] = deval(sol,sol.x(end)); % updating y0 and yp0 with last

        MetricTensorSharp = [MetricTensorSharp;y0(1:6)'];
        kappa = [kappa; y0(7)];
%         Tau = [Tau;y0(8)];
            Tau = [Tau;params.appliedTension];
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

T = time*AssumedTaut;

%% Plotting creep
h = h0./kappa;

NormalisingTau = rho0*xiap; % apical tension
idx = find(T > 500, 1, 'first');
% idx = 1;
t = T(idx:end) - T(idx);
tau = Tau(idx:end);
stretch = kappa(idx:end) - kappa(idx) + 1;
h = h(idx:end);

figure(2)
hold on
plot(t, tau/NormalisingTau, 'linewidth',2)
xlabel('Time [s]'); ylabel('Tension');

figure(1)
hold on
plot(t, stretch - 1, 'linewidth',2)
plot(t,h, 'linewidth',2)
xlabel('Time [s]'); ylabel('Strain');

figure(3)
loglog(t, stretch - 1, 'linewidth',2)
hold on
xlabel('Time [s]'); ylabel('Strain');


toc