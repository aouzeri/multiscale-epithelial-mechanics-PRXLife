%% Plotting results from AnalyzeVTK.m for (cyclic) domes 
% Copyright (c) Adam Ouzeri 2025

clear all;
close all;

%% Run AnalyzeVTK.m first to obtain the values
AnalyzeVTK.m


%% Windows & Ubuntu
timedatafoldername = '[pathtosimulationfolder]'; tequi = 10;

tau_turnover = 38; % fast, middle, slow

foldername = [timedatafoldername, '/AnalyzeVTKpacket/'];


load([foldername,'facetype.mat']); % 'facetype'
load([foldername,'FaceArea.mat']); % 'faceareadata'
load([foldername, 'Gflat_per_face_and_step.mat' ]); % 'faceGfdata'
load([foldername, 'Rho_per_face_and_step.mat']); % 'facerhodata'

timefilename = sprintf('%sTimeData.txt',timedatafoldername);
timedata     = load(timefilename);
% TIME = [0; timedata(1:size(faceareadata,2)-1,1)];
% time = TIME;

%% Shifting time such that t0 = tequi
Nsteps = size(timedata,1);
relevantidx = (find(timedata >= tequi,1) - 1) : Nsteps;
time = timedata(relevantidx - 1);

faceareadata = faceareadata(:,relevantidx);
facerhodata = facerhodata(:,relevantidx);
faceGfdata = faceGfdata(:,:,relevantidx);
% time = TIME(relevantidx); time = time - time(1);


basalFaceID = 985 + 1; % for domes cyclic pressure
apicalFaceID = 986 + 1; % for domes cyclic pressure


%% Calculating trace and determinant of G
traceG = zeros(size(faceGfdata,1),size(faceGfdata,3));
detG = zeros(size(faceareadata,1),length(time));
for i = 1:length(time)
    traceG(:,i) = sum(faceGfdata(:,[1,4,6],i),2);
    for face = 1:size(faceareadata,1)
        G = [faceGfdata(face,1,i), faceGfdata(face,2,i), faceGfdata(face,3,i); ...
            faceGfdata(face,2,i), faceGfdata(face,4,i), faceGfdata(face,5,i); ...
            faceGfdata(face,3,i), faceGfdata(face,5,i), faceGfdata(face,6,i)];
        detG(face,i) = det(G(1:2,1:2));
    end
    
end

%% Plotting
time = (time - tequi) * tau_turnover;

figure(1)
yyaxis left
plot(time,faceareadata(basalFaceID,:)./faceareadata(basalFaceID,1) - 1, 'linewidth', 2, 'DisplayName', '"Top" cell basal area change')
xlabel('Time [s]'); ylabel('Area change')
hold on
plot(time,faceareadata(apicalFaceID,:)./faceareadata(apicalFaceID,1) - 1, 'linewidth', 2, 'DisplayName', '"Top" cell apical area change')

figure(2)
plot(time,traceG(basalFaceID,:), 'linewidth', 2, 'DisplayName', 'Tr(G)')
xlabel('Time [s]'); ylabel('traceG(G)')
hold on

figure(3)
plot(time,facerhodata(basalFaceID,:)/facerhodata(basalFaceID,1), 'linewidth', 2, 'DisplayName', 'Density change basal' )
xlabel('Time [s]'); ylabel('\rho/\rho_0')
hold on
plot(time,facerhodata(apicalFaceID,:)/facerhodata(apicalFaceID,1), 'linewidth', 2, 'DisplayName', 'Density change apical' )

 
figure(5)
yyaxis left
% plot(time,faceareadata(basalFaceID,:)./faceareadata(basalFaceID,1), 'linewidth', 2, 'DisplayName', 'Apical area change basal')
xlabel('Time [s]'); ylabel('Area change, sqrt(det(G))')
hold on
plot(time,faceareadata(basalFaceID,:)./faceareadata(basalFaceID,1), 'k-', 'linewidth', 2, 'DisplayName', '"Top" cell basal area change')
plot(time,faceareadata(apicalFaceID,:)./faceareadata(apicalFaceID,1), 'k-', 'linewidth', 2, 'DisplayName', '"Top" cell apical area change')
plot(time,sqrt(detG(basalFaceID,:))./sqrt(detG(basalFaceID,1)), 'b--', 'linewidth', 2, 'DisplayName', 'det(G) basal')
plot(time,sqrt(detG(apicalFaceID,:))./sqrt(detG(apicalFaceID,1)), 'b--', 'linewidth', 2, 'DisplayName', 'det(G) apical')

% plot(time,traceG(basalFaceID,:)/2, 'linewidth', 2, 'DisplayName', 'Tr(G) basal')
% plot(time,traceG(apicalFaceID,:)/2, 'linewidth', 2, 'DisplayName', 'Tr(G) apical')
% plot(time,facerhodata(basalFaceID,:)/facerhodata(basalFaceID,1), 'linewidth', 2, 'DisplayName', 'Density change basal' )
% plot(time,facerhodata(apicalFaceID,:)/facerhodata(apicalFaceID,1), 'linewidth', 2, 'DisplayName', 'Density change apical' )

%% Plotting Laplace law results for entire dome
load([timedatafoldername, 'workspace.mat'], 'height', 'lumpressure', 'basalradius')
relevantidx2 = (find(lumpressure > 0,1) - 1) : Nsteps;

HEIGHT = height(relevantidx);
HEIGHT = HEIGHT - HEIGHT(1);
LUMPRESSURE = lumpressure(relevantidx)*8000;
domeRadius = (HEIGHT.^2 + basalradius.^2)./(2*HEIGHT);
strain = (HEIGHT./basalradius).^2;
DomeTension = domeRadius.*LUMPRESSURE./2/8000/0.03;


figure(1)
yyaxis left
plot(time, strain, 'k', 'LineWidth', 2,  'DisplayName', 'Dome strain')

figure(4)
hold on
plot(strain, DomeTension, 'LineWidth', 2, 'DisplayName', 'Laplace law tension')
% plot(strain, stretchedTissueTension, 'LineWidth', 2)
title('Tension vs strain')

figure(5)
hold on
yyaxis right
plot(time, LUMPRESSURE, 'LineWidth', 1, 'DisplayName', 'Pressure')

figure(1)
yyaxis right
plot(time, LUMPRESSURE, 'LineWidth', 1, 'DisplayName', 'Pressure')

figure(6)
hold on
yyaxis left
plot(time, DomeTension, 'LineWidth', 2, 'DisplayName', 'Tension')
yyaxis right
plot(time, LUMPRESSURE, 'LineWidth', 2, 'DisplayName', 'Pressure')
title('Tension vs time')

figure(7)
hold on
yyaxis left
plot(relevantidx, DomeTension, 'LineWidth', 2,  'DisplayName', 'Tension')
yyaxis right
plot(relevantidx, LUMPRESSURE, 'LineWidth', 2,  'DisplayName', 'Pressure')

figure(8)
hold on
yyaxis left
plot(relevantidx, domeRadius, 'LineWidth', 2,  'DisplayName', 'Radius')
yyaxis right
plot(relevantidx, LUMPRESSURE, 'LineWidth', 2,  'DisplayName', 'Pressure')
