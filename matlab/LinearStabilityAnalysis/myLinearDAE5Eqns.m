%% Copyright (c) Adam Ouzeri 2025

%% Observing the dyanmics of the linear system
% Cs are all flat
% Gs are all sharp
clear
close all
% clc;
tic
%% Creating the symbolic variables
syms Gxy Gxylat Gz Nap NCorlat NAdhlat epsilon Cxy Cz
syms Iap Jap Ilat Jlat
syms kp kd sigma eta Aap0 lambda mu xiCorlat xiAdhlat h0 ah kon koff
% @(Gap, Gap, Gxy, Gz, Nap, Nlat, epsilon )

Cxy = epsilon + 1;
Cz = 1/(epsilon  + 1)^2;
Aap = Aap0*(epsilon + 1);
Alat = h0*Aap0*ah/sqrt(Aap);

%% Calculating the matrix
GapMat = [Gxy, 0; 0, Gxy];
CapMat = [Cxy, 0 ; 0 , Cxy];

GlatMat = [Gxylat, 0; 0, Gz];
ClatMat = [Cxy, 0 ; 0 , Cz];

Iap = sum(sum(CapMat.*GapMat));
Jap = sqrt(det(CapMat)*det(GapMat));

Ilat = sum(sum(ClatMat.*GlatMat));
Jlat = sqrt(det(ClatMat)*det(GlatMat));

dPsidIap = mu/Jap;
dPsidJap = 2*lambda*(Jap - 1) - mu*Iap/Jap^2;
dPsidIlat = mu/Jlat;
dPsidJlat = 2*lambda*(Jlat - 1) - mu*Ilat/Jlat^2;
alpha = h0*ah/sqrt(Aap0);

f1 = -Jap*Gxy/eta*(dPsidIap + dPsidJap/2);
f2 = -Gxylat/eta*((epsilon + 1) * Gxylat*dPsidIlat + Jlat/2*dPsidJlat);
f3 = -Gz*Gz/eta*(mu/Jlat*Cz + (-mu*Ilat/Jlat^2 + 2*lambda*(Jlat-1))*Jlat/Gz/2);
f4 = kp*Aap - kd*Nap;
f5 = kp*Alat - kd*NCorlat;
f6 = kon*Alat - koff*NAdhlat;

g = (2*Nap/Aap - alpha*(NCorlat*xiCorlat + NAdhlat*xiAdhlat)/Alat/(2*(epsilon + 1)^(3/2))) + 2*Nap/Aap*Jap*(2*dPsidIap + dPsidJap)...
    + alpha*NCorlat/Alat*((sqrt(Gxylat/Gz) - 2/(epsilon + 1)^3*sqrt(Gz/Gxylat))*dPsidIlat - dPsidJlat/2/(epsilon + 1)^(3/2)) - sigma;

f = [f1,f2,f3,f4,f5,f6];

% General form of DAE at fixed point
F = [f,g];
DFx01 = diff(F,Gxy);
DFx02 = diff(F,Gxylat);
DFx03 = diff(F,Gz);
DFx04 = diff(F,Nap);
DFx05 = diff(F,NCorlat);
DFx06 = diff(F,NAdhlat);
DFx07 = diff(F,epsilon);
Dx0F = [DFx01', DFx02', DFx03', DFx04', DFx05', DFx06', DFx07'];


% Building the matrix dfx - dfy*(dgy)^-1*dgx
dfx1 = diff(f,Gxy);
dfx2 = diff(f,Gxylat);
dfx3 = diff(f,Gz);
dfx4 = diff(f,Nap);
dfx5 = diff(f,NCorlat);
dfx6 = diff(f,NAdhlat);


dxf = [dfx1', dfx2', dfx3', dfx4', dfx5', dfx6'];

dyf = diff(f,epsilon);
dyf = dyf';
dyg =  diff(g,epsilon);

dgx1 = diff(g,Gxy);
dgx2 = diff(g,Gxylat);
dgx3 = diff(g,Gz);
dgx4 = diff(g,Nap);
dgx5 = diff(g,NCorlat);
dgx6 = diff(g,NAdhlat);
dxg = [dgx1',dgx2',dgx3',dgx4',dgx5', dgx6'] ;

% symMat = dxf - dyf/dyg*dxg;
% symEig =  eig(symMat);

dxf = matlabFunction(dxf);
dyf = matlabFunction(dyf);
dxg = matlabFunction(dxg);
dyg = matlabFunction(dyg);

f = matlabFunction(f);
g = matlabFunction(g);
F = matlabFunction(F);
Dx0F = matlabFunction(Dx0F);



%% Parameters
% normalising parameters
kd = 1;                 % depolymerization rate [1/s]
s0 = 1;                 % cell dimension [micro.m]
xiap = 1;               % apical stress in reference configuraion, [N/micro.m2]


%% Turnover parameters
kp = 0.03;           % polymerisation rate; kp/kd = rho0 = 0.03; [micro.m]; need to keep it two orders of magnitude smaller than cell dimension
koff = 0.1;          % dictates ahdesion turnover timescale : koff / kd = 0.1 (an order of magnitude smaller)
kon = koff;          % arbitrary. Setting the adhesion density to 1 and then changing the adhesion tension in consequence;
rho0 = kp/kd;
adhesionDens0 = kon/koff; % equilibrium adhesion concentration

%% Viscoelastic parameters
% lambdaV = linspace(0.01,0.05,500);  etaV = 0.2*ones(size(lambdaV));           % bulk modulus; mu/lambda = 1; [N/micro.m2]
etaV = linspace(0.001,1,1000);  lambdaV = linspace(0.001,1.0,1); % viscous remodeling;  [N.s/m2]; need to keed remodelling timescale one order of magnitude smaller to kd
% save('eta500.mat', 'etaV')
muV = lambdaV;               % shear modulus [N/micro.m2]

%% Geometric parameters
h0 = 2;
Aap0 = 3*sqrt(3)/2*s0^2;            % reference cell area, [micro.m^2]
Peri0 = 6*s0;
Apar0 = s0*h0;
Atra0 = s0*h0;
Alat0 = 2*Apar0 + 4*Atra0;
Acell_0 = 2*Aap0 + Alat0;			% reference total cell area, [micro.m^2]
ah = sqrt(72/3/sqrt(3));

%% Tensional parameters
TotalLateralTension = 0.5;
xiCorlat = TotalLateralTension;
xiAdhlat = (TotalLateralTension - xiCorlat)*rho0/adhesionDens0;

%% calculating the corresponding steady-state epsilon
sigmaV = linspace(0.01*rho0,1.99*rho0,1000); % sigmaV = 0.033;
% save('sigma500.mat','sigmaV')
epsilonV = zeros(size(sigmaV));
options = optimoptions('fsolve','Display','iter', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 10000, ...
    'OptimalityTolerance', 1e-18, 'Tolfun', 1e-18);

for n = 1:length(epsilonV)
    fun = @(x) solveSSforStrain(x,rho0,adhesionDens0,h0,ah,Aap0,xiCorlat, xiAdhlat, sigmaV(n));
    x0 = 0;
    epsilonV(n) = fsolve(fun,x0,options);
end
%     save('epsilon100.mat','epsilonV')
%     load('epsilon100.mat')

%% calculating the corresponding steady-state tension

%     epsilonV = 0;
%     sigmaV = zeros(size(epsilonV));
%     options = optimoptions('fsolve','Display','iter', 'MaxFunctionEvaluations', 50000, 'MaxIterations', 10000, ...
%         'OptimalityTolerance', 1e-18, 'Tolfun', 1e-18);
%
%     for n = 1:length(epsilonV)
%         fun = @(x) solveSSforTension(x,epsilonV,rho0,h0,ah,Aap0,xilatV);
%         x0 = 0;
%         sigmaV(n) = fsolve(fun,x0,options);
%         toc
%     end


%% performing the analysis for each parameter value
iter = 0;
bifpointsReduced = zeros(length(sigmaV),length(etaV),length(lambdaV));
Eig = zeros(6,1);
for l = 1:length(lambdaV)
    for i = 1:length(etaV)
        display(['Step (l,i) : ', num2str(l), ', ', num2str(i)])

        for j = 1:length(sigmaV)
            iter = iter + 1;
            
            %% Calculating the fixed point
            %             lambda = lambdaV(j); mu = muV(j);
            lambda = lambdaV(l); mu = muV(l);
            eta = etaV(i);
            sigma = sigmaV(j); epsilon = epsilonV(j);
            %             sigma = sigmaV(1); epsilon = epsilonV(1);
            
            Cxy = epsilon + 1;
            Cz = 1/(epsilon  + 1)^2;
            
            Gxy = 1/Cxy;
            Gxylat = 1/Cxy;
            Gz = 1/Cz;
            
            % equilibrium amounts for a given strain
            Nap = double(rho0*(Aap0*(epsilon + 1)));
            NCorlat = double(rho0*(h0*Aap0*ah/sqrt(Aap0*(epsilon + 1))));
            NAdhlat = double(adhesionDens0*(h0*Aap0*ah/sqrt(Aap0*(epsilon + 1))));
            
            Dxf = double(dxf(Gxy,Gxylat,Gz,epsilon,eta,kd,koff,lambda,mu));
            Dyf = double(dyf(Aap0,Gxy,Gxylat,Gz,ah,epsilon,eta,h0,kon,kp,lambda,mu));
            Dyg = double(dyg(Aap0,Gxy,Gxylat,Gz,NAdhlat,NCorlat,Nap,epsilon,lambda,mu,xiAdhlat,xiCorlat));
            Dxg = double(dxg(Aap0,Gxy,Gxylat,Gz,NCorlat,Nap,epsilon,lambda,mu,xiAdhlat,xiCorlat));
            Mat = Dxf - Dyf/Dyg*Dxg;
            [V,D] = eig(Mat);
            eigen = diag(D);
            Eig(:,iter) = sort(eigen);
            
            %                 symMat = sym(Mat);
            %                 [symVJ, symJordan] = jordan(symMat);
            %                 [VJ, Jordan] = jordan(Mat);
            
            % for bifurcation diagram
            spiraling = 0;
            stable = 0;
            
            % doing checks on these complex eigenvalues
            imagIdxes = find(abs(imag(eigen)) > 1e-10);
            if ~isempty(imagIdxes)
                if length(imagIdxes) == 2
                    complexConjugates = abs(imag(eigen(imagIdxes(1))) + imag(eigen(imagIdxes(2)))) < 1e-10;
                    if ~complexConjugates
                        error('The two eigen values are NOT complex conjugates')
                    end
                elseif rem(length(imagIdxes), 2) == 1
                    error('You have ODD number of complex eigenvalues')
                end
                spiraling = 1;
            end
            
            if ~ismember(1,sign(real(eigen)))
                stable = 1;
            end
            
            if abs(real(eigen)) < 1e-10
                error('You have a zero real part eigen value.')
            end
            
            %% Creating phase diagram
            
            if spiraling
                if stable
                    bifpointsReduced(j,i,l) = 1;
                else
                    bifpointsReduced(j,i,l) = 3;
                end
            else
                if stable
                    bifpointsReduced(j,i,l) = 2;
                else
                    bifpointsReduced(j,i,l) = 4;
                end
            end
            
            %% checking that we have the fixed point and that the system can be reduced
            
            fMat = double(f(Aap0,Gxy,Gxylat,Gz,NAdhlat,NCorlat,Nap,ah,epsilon,eta,h0,kd,koff,kon,kp,lambda,mu));
            gMat = double(g(Aap0,Gxy,Gxylat,Gz,NAdhlat,NCorlat,Nap,epsilon,lambda,mu,sigma,xiAdhlat,xiCorlat));
            dgyMat = double(dyg(Aap0,Gxy,Gxylat,Gz,NAdhlat,NCorlat,Nap,epsilon,lambda,mu,xiAdhlat,xiCorlat));
            
            if find(abs(fMat) > 1e-9 | abs(gMat) > 1e-9)
                error('The given steady state point not a fixed point of the system')
            end
            if abs(dgyMat) < 1e-7
                warning('DgDy is close to 0')
            end
            
            %% checking that we have regular DAE
            DFx0Mat = double(Dx0F(Aap0,Gxy,Gxylat,Gz,NAdhlat,NCorlat,Nap,ah,epsilon,eta,h0,kd,koff,kon,kp,lambda,mu,xiAdhlat,xiCorlat));
            if rank(DFx0Mat) ~= length(DFx0Mat)
                error('The given DAE is not regular')
            end
        end
    end
end

%% Creating the isosurface
% changing the matrix to 0 and 1
if length(lambdaV) > 1
    %     save('bifmatrix3D.mat','bifpointsReduced')
    binaryBifPoints = bifpointsReduced;
    binaryBifPoints(binaryBifPoints < 3) = 0;
    binaryBifPoints(binaryBifPoints >= 3) = 1;
    modifiedBifPoints = binaryBifPoints;
    
%     for y = 1:length(etaV)
%         for z = 1:length(lambdaV)
%             modifiedBifPoints(:,y,z) = [abs(diff(binaryBifPoints(:,y,z)));0] ;
%         end
%     end
    
    
    [X,Y,Z] = meshgrid(etaV,sigmaV/rho0,lambdaV);
    figure(3)
    %     surf(X,Y,Z,modifiedBifPoints)
    pa = patch(isosurface(X,Y,Z,modifiedBifPoints,0));
    isonormals(X,Y,Z,modifiedBifPoints,pa)
    pa.FaceColor = 'red';
    pa.EdgeColor = 'none';
    daspect([1 1 1])
    view(3);
    axis tight
    camlight
    lighting gouraud
    
%     load('bifmatrix_froml269_l500_i500_j500.mat')
%     binaryBifPoints = bifpointsReduced;
%     binaryBifPoints(binaryBifPoints < 3) = 0;
%     binaryBifPoints(binaryBifPoints >= 3) = 1;
%     modifiedBifPoints = binaryBifPoints;
%     
%     for y = 1:length(etaV)
%         for z = 1:length(lambdaV)
%             modifiedBifPoints(:,y,z) = [abs(diff(binaryBifPoints(:,y,z)));0] ;
%         end
%     end
%     
%     
%     [X,Y,Z] = meshgrid(etaV,sigmaV/rho0,lambdaV);
%     figure(3)
%     %     surf(X,Y,Z,modifiedBifPoints)
%     pa = patch(isosurface(X,Y,Z,modifiedBifPoints,0));
%     isonormals(X,Y,Z,modifiedBifPoints,pa)
%     pa.FaceColor = 'red';
%     pa.EdgeColor = 'none';
%     daspect([1 1 1])
%     view(3);
%     axis tight
%     camlight
%     lighting gouraud
    
    %     % 2D verification
    %     Dpoint = modifiedBifPoints(:,:,230);
    %     [x2,y2] = meshgrid(etaV,sigmaV/rho0);
    %     figure(10)
    %     surf(x2,y2,Dpoint)
end
toc

%% creating plot

% load('sigma500.mat')
% load('epsilon500.mat')
% load('bifmatrixReduced5004regions.mat')

%         close all
%     figure(100)
%     hold on
%     % ytickformat(gca, 'percentage')
%     contourf(etaV,lambdaV,bifpointsReduced')
%     % imagesc('XData', etaV, 'YData',sigmaV/0.03, 'CData', bifpointsReduced')
%     colormap([1 0 0; 0.5 0 0; 0 0 1; 0 0 0.5])
%     % colormap([0 0 1; 0 0 0.5])

figure(1)
hold on
%     subplot(3,1, k)
% ytickformat(gca, 'percentage')
contourf(etaV,sigmaV/rho0,bifpointsReduced(:,:,1))
%     contourf(etaV,lambdaV,bifpointsReduced')
%     % imagesc('XData', etaV, 'YData',sigmaV/0.03, 'CData', bifpointsReduced')
colormap([255 150 0; 255 0 0; 166 217 106; 0 130 0]/255)
%     colormap([0.5 0 0; 0 0 0.5])
% title(['Tension : \xi_l = ', num2str(xilatV)])
%     title(['h = ', num2str(h0)])
%
% % close all;
% % save('bifmatrixReduced5004regions.mat','bifpointsReduced')
% subplot(3,3, 3*(k-1) + 2)
% % ytickformat(gca, 'percentage')
% contourf(etaV,epsilonV,bifpointsReduced')
% % imagesc('XData', etaV, 'YData',sigmaV/0.03, 'CData', bifpointsReduced')
% colormap([1 0 0; 0.5 0 0; 0 0 1; 0 0 0.5])
% % title(['Strain : \xi_l = ', num2str(xilatV)])
%
% subplot(3,3, 3*k)
% % ytickformat(gca, 'percentage')
% contourf(etaV,TangModV,bifpointsReduced')
% % imagesc('XData', etaV, 'YData',sigmaV/0.03, 'CData', bifpointsReduced')
% colormap([1 0 0; 0.5 0 0; 0 0 1; 0 0 0.5])
% % title(['Stiffness : \xi_l = ', num2str(xilatV)])

%% Getting trajectories
% k = V\ones(4,1)*1e-4;
% time = linspace(0,100,10000);
% xSol = zeros(4,length(time));
% phit = xSol;
% for i = 1:length(time)
%     t = time(i);
%     xSol(:,i) = k(1)*exp(l(1)*t)*V(:,1) + k(2)*exp(l(2)*t)*V(:,2)+ k(3)*exp(l(3)*t)*V(:,3)+ ...
%         k(4)*exp(l(4)*t)*V(:,4);
% end
% plot(time,xSol,'linewidth', 3);

% syms x1(t) x2(t) x3(t) x4(t) x5(t)
% X = [x1;x2;x3;x4;x5];
% initcond = X(0) == xinit;
% odes = diff(X) == Mat*X + FP;
% [x1Sol(t), x2Sol(t),x3Sol(t),x4Sol(t),x5Sol(t)] = dsolve(odes, initcond);

%% Quasistatic analysis
figure(2)
strainrange = -0.5:0.01:5;
tauquasistatic = (2 - xiCorlat*Aap0*h0*ah/2./(Aap0*(strainrange + 1)).^1.5);
%     subplot(3,3,8)
plot(strainrange, tauquasistatic, 'linewidth', 3)

% Plotting the fixed points stability
eta = 0.3;
[~,closestIdx] = min(abs(etaV - eta));
e1 = epsilonV(find(bifpointsReduced(:,closestIdx) == 2, 1, 'last'));
e2 = epsilonV(find(bifpointsReduced(:,closestIdx) == 1, 1, 'last'));
e3 = epsilonV(find(bifpointsReduced(:,closestIdx) == 3, 1, 'last'));
strainStable = linspace(-0.5,e1,100);
strainStableOsc = linspace(e1,e2,100);
strainUnstableOsc = linspace(e2,e3,100);
strainUnstable = linspace(e3,10,100);
V0 = Aap0*h0;

tauquasistaticStable = (2 - xiCorlat*V0*ah/2./(Aap0*(strainStable + 1)).^1.5);
tauquasistaticStableOsc = (2 - xiCorlat*V0*ah/2./(Aap0*(strainStableOsc + 1)).^1.5);
tauquasistaticUnstableOsc = (2 - xiCorlat*V0*ah/2./(Aap0*(strainUnstableOsc + 1)).^1.5);
tauquasistaticUnstable = (2 - xiCorlat*V0*ah/2./(Aap0*(strainUnstable + 1)).^1.5);

figure(10)
hold on
xtickformat(gca, 'percentage');
plot(100*strainStable, tauquasistaticStable, 'linewidth', 3)
plot(100*strainStableOsc, tauquasistaticStableOsc, 'linewidth', 3)
plot(100*strainUnstableOsc, tauquasistaticUnstableOsc, 'linewidth', 3)
plot(100*strainUnstable, tauquasistaticUnstable, 'linewidth', 3)
%     ylim([0 2])
%     xlim([0 max(strainUnstable*100)])

% figure
% plot(etaV, real(Eig(1,:)), 'linewidth', 3)
% 
% figure
% plot(etaV, imag(Eig(1,:)), 'linewidth', 3)


%% Getting steady-state tension or strain

function F = solveSSforStrain(x,rho0,c0,h0,ah,Aap0,xiCorlat, xiAdhlat, sigma)
F = sigma/rho0 - 2 + (xiCorlat + xiAdhlat*c0/rho0) *h0*Aap0*ah/(2*(Aap0*(x + 1))^(3/2));
end

function F = solveSSforTension(x,strain,rho0,h0,ah,Aap0,xilat)

F = x/rho0 - 2 + xilat*h0*Aap0*ah/(2*(Aap0*(strain + 1))^(3/2));
end
