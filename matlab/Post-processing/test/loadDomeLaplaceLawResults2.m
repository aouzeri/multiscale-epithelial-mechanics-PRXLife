%% Reading and plotting results saved from DomeLaplaceLaw
% Copyright (c) 2025 Adam Ouzeri

clear all;
close all;

%% calculating quantities
relevantidx = 1:nSteps;

HEIGHT = height(relevantidx);
HEIGHT = HEIGHT - HEIGHT(1);
LUMPRESSURE = lumpressure(relevantidx);
LUMPRESSURE = lumpressure(end);
domeRadius = (HEIGHT.^2 + basalradius.^2)./(2*HEIGHT);
strain = (HEIGHT./basalradius).^2;
DomeTension = domeRadius.*LUMPRESSURE./2/0.03;
STRETCHEDTISSUETENSION = stretchedTissueTension(relevantidx);
TIMEDATA = timedata(relevantidx) - tequi;

% Estimating tissue surface tension from idealised tissue equibiaxial stretching
s0 = 1;
h0 = 2;
Aap0 = 3*sqrt(3)/2*s0^2;            % reference cell area, [micro.m^2]
ah = sqrt(72/3/sqrt(3));
K = ah*h0/sqrt(Aap0)/2;
xi_a = 1.0;
xi_b = 1.0;
xi_l = 0.7;

stretchedTissueTension =  xi_a + xi_b - K*xi_l./(strain + 1).^(3/2);

figure(1)
hold on
plot(strain, DomeTension, 'LineWidth', 2)
title('Tension vs strain')

figure(2)
hold on
plot(strain, domeRadius, 'LineWidth', 2)
title('Radius vs strain')

figure(3)
hold on
plot(TIMEDATA, domeRadius, 'LineWidth', 2)
title('Radius vs time')

figure(9)
hold on
plot(TIMEDATA, strain, 'LineWidth', 2)
title('Strain vs time')

figure(5)
hold on
plot(TIMEDATA, DomeTension, 'LineWidth', 2)
title('Tension vs time')

figure(6)
hold on
plot(TIMEDATA, HEIGHT, 'LineWidth', 2)
title('Height vs time')

figure(7)
hold on
plot(LUMPRESSURE, strain, 'LineWidth', 2)
title('Strain vs pressure')

figure(8)
hold on
plot(TIMEDATA, LUMPRESSURE, 'LineWidth', 2)
title('Pressure vs time')
