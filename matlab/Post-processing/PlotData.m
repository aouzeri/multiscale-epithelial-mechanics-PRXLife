%% Copyright (c) Adam Ouzeri and Sohan Kale 2025

%% Plotting stress strain data
clc;
clear all;
% close all;

foldername = '[pathtosimulationfolder]/'; tequi = 0;

%%

k_p = 1;
k_d = 1;
c0  = 0.03;
rho_0 = k_p * c0/ k_d;

filename = sprintf('%s/AvgStressStrainData.txt',foldername);
data = dlmread(filename);
nSteps = size(data,1);

Favg      = data(:,1:4);
Pavg      = data(:,5:8)/rho_0;
Savg      = data(:,9:12)/rho_0;
h         = data(:,13);
Fzz       = data(:,14);

% load timesteps data
timefilename = sprintf('%s/TimeData.txt',foldername);
data2    = dlmread(timefilename);
timedata = data2(:,1);

timedata = timedata - tequi;

% Get areal strain
Javg = zeros(nSteps,1);
for i = 1:nSteps
	Fmat = [Favg(i,1), Favg(i,2); Favg(i,3), Favg(i,4)];
	Javg(i) = det(Fmat);
end
eps_A = Javg - 1;


%% Plototing
% figure(1);hold on; plot(eps_A,Savg(:,1),'c-','linewidth',2);
% figure(1);hold on; plot(eps_A,Savg(:,4),'c-','linewidth',2);
figure(1);hold on; plot(eps_A,(Savg(:,1) + Savg(:,4))/2,'linewidth',2);

% figure(2);hold on; plot(timedata(:,1),Savg(:,1),'c--','linewidth',2); xlabel('Time [-]'); ylabel('Tension [-]');
% figure(2);hold on; plot(timedata(:,1),Savg(:,4),'c-','linewidth',2); xlabel('Time [-]'); ylabel('Tension [-]');
figure(2);hold on; plot(timedata(:,1),(Savg(:,1) + Savg(:,4))/2,'linewidth',2); xlabel('Time [-]'); ylabel('Tension [-]');

