%% Computing stress-strain curve for domes using Laplace's law
% We first need to extract the pressure and the height of the dome from the
% vtk files. To obtain the height estimate, we take the heighest
% point when the dome is formed (last timestep) and extract the height at
% each timestep.

% Adam Ouzeri (c) April 2025

% clc;
clear;
% close all; 

%Input folder where the simulations were done
foldername = ['[pathtosimulationfolder]']; 

domeCenter = [-0.33717, -0.19041];

% load timesteps data
timefilename = sprintf('%s/TimeData.txt',foldername);
data    = dlmread(timefilename);
timedata = data(:,1);
nSteps = size(data,1);
tequi = 10;
relevantidx = (find(timedata >= tequi,1)) : nSteps; relevantidx = 1:nSteps;

%% Obtaining number of point (all our VTK have the same format: POINTS nPoins float)
filename = sprintf('%ssolution%s.vtk',foldername, num2str(nSteps - 1));
fileID = fopen(filename, 'r');
fulltext = textscan(fileID, '%s', 'delimiter', '\n');
fclose(fileID);
fulltextStr = string(fulltext{1});
pointLineString = fulltextStr{6};
index_of_f = strfind(pointLineString, 'f');
% ending two characters before the letter "f"
nPoints = str2num(pointLineString(8:index_of_f-2));

%% Obtaining lineID for point to calculate height and an estimate of the dome basal radius
heightLineStarts = strfind(fulltextStr, 'DOF_2');
xLineStarts = strfind(fulltextStr, 'DOF_0');
yLineStarts = strfind(fulltextStr, 'DOF_1');
pressureLineStarts = strfind(fulltextStr, 'AUX_7');
densityLineStarts = strfind(fulltextStr, 'AUX_0');


for i = 1:length(heightLineStarts)
    if heightLineStarts{i} == 1
        linenumHeight = i;
        break;
    end
end

for i = 1:length(xLineStarts)
    if xLineStarts{i} == 1
        linenumX = i;
        break;
    end
end

for i = 1:length(yLineStarts)
    if yLineStarts{i} == 1
        linenumY = i;
        break;
    end
end

for i = 1:length(pressureLineStarts)
    if pressureLineStarts{i} == 1
        linenumPressure = i;
        break;
    end
end

for i = 1:length(densityLineStarts)
    if densityLineStarts{i} == 1
        linenumDensity = i;
        break;
    end
end

% storing heigh of all points at the last time step
pointsheight = zeros(nPoints,1);
for i=1:nPoints
    pointsheight(i) = str2double(fulltextStr{linenumHeight + i});
   % pointsheight(i)
end

% getting maximum height 
[~, maxIdx] = max(pointsheight);
pointlineID = maxIdx;

% getting corresponding point on basal face
idxofspaces = find(isspace(fulltextStr{6 + pointlineID}));
xyzCoordsofPointtoFollow = strfind(fulltextStr,fulltextStr{6 + pointlineID}(1:idxofspaces(2)-1));
for i = 1:length(xyzCoordsofPointtoFollow)
    if xyzCoordsofPointtoFollow{i} == 1
        basalpointCounterPart_lineID = i;
        if(~str2num(fulltextStr{basalpointCounterPart_lineID}(idxofspaces(2):end)))
            break;
        end
    end
end
basalpointCounterPart_lineID = basalpointCounterPart_lineID - 6;

% get all basal nodes coordinates
basalCoords = zeros(nPoints,2);
for i = 1:nPoints
    basalCoords(i,:) = [str2double(fulltextStr{linenumX + i}), str2double(fulltextStr{linenumY + i})];
end
tissueCenterCoords = mean(basalCoords);

% getting the fixed basal nodes
idx_0 = find(pointsheight < 0.0001);
fixedBasalNodescoords = zeros(length(idx_0),2);
(i)
for i=1:length(idx_0)
    fixedBasalNodescoords(i,:) = [str2double(fulltextStr{linenumX + idx_0(i)}), str2double(fulltextStr{linenumY + idx_0(i)})];
end

% taking the minimum distance center as basal radius;
distancesFromDomeCenter = sort(sqrt((fixedBasalNodescoords(:,1) - domeCenter(1)).^2 + (fixedBasalNodescoords(:,2) - domeCenter(2)).^2));
xDistancesFromTissueCenter = sort(sqrt((fixedBasalNodescoords(:,1) - domeCenter(1)).^2));
yDistancesFromTissueCenter = sort(sqrt((fixedBasalNodescoords(:,2) - domeCenter(2)).^2));
disp('Basal radius from current measurement')
basalradius = min(distancesFromDomeCenter) 
disp('Basal radius from mesh generation')
basalradius = 7.2676 % for 15x15 at 0.6*l/2

%% Obtaining height and luminal pressure for all timesteps
apicalheight = zeros(nSteps,1);
basalheight = zeros(nSteps,1);
delta_p = zeros(nSteps,1);
density = zeros(nSteps,1);
disp('Extracting results at each time step...')

% Storing relevant lines
basHeightLine = linenumHeight + basalpointCounterPart_lineID - 1
apHeightLine = linenumHeight + pointlineID - 1
pressLine = linenumPressure + pointlineID - 1;
densLine = linenumDensity + pointlineID - 1;

apHeightLine_basHeightLine = apHeightLine - basHeightLine - 1
pressLine_apHeightLine = pressLine - apHeightLine - 1;
densLine_pressLine = densLine - pressLine - 1;

for step = 1:nSteps
    fileStep = num2str(step - 1);
    filename2 = sprintf('%s/solution%s.vtk',foldername, fileStep);  
    
    % collecting basal height, then apical height, then pressure, then
    % density
    fileID = fopen(filename2, 'r');
    filetext = textscan(fileID, '%s',1,'delimiter','\n', 'headerlines', basHeightLine);
    basalheight(step) = str2double(filetext{1});
    
    filetext = textscan(fileID, '%s',1,'delimiter','\n', 'headerlines', apHeightLine_basHeightLine);
    apicalheight(step) = str2double(filetext{1});
    
    filetext = textscan(fileID, '%s',1,'delimiter','\n', 'headerlines', pressLine_apHeightLine);
    delta_p(step) = str2double(filetext{1});
    
%     filetext = textscan(fileID, '%s',1,'delimiter','\n', 'headerlines', densLine_pressLine);
%     density(step) = str2double(filetext{1});  
    
    fclose(fileID) ;
   
    disp(['Step ', num2str(step)])
    
end

% Adjusting values to tequi = 10
zPoint = (apicalheight + basalheight)/2;
% zPoint = basalheight;
% zPoint = apicalheight;

% density = density(relevantidx);


%% Post-processing
height = zPoint(relevantidx); % height = height - height(1);
lumpressure = -delta_p(relevantidx);

% Computing radius
radius = (height.^2 + basalradius.^2)./(2*height);

% Computing dome strain
domestrain = (height./basalradius).^2;

% Computing tension
tension = radius.*lumpressure/2;

% Estimating tissue surface tension from idealised tissue equibiaxial stretching
s0 = 1;
h0 = 2;
Aap0 = 3*sqrt(3)/2*s0^2;            % reference cell area, [micro.m^2]
ah = sqrt(72/3/sqrt(3));
K = ah*h0/sqrt(Aap0)/2;

% if(strcmpi('NoDiffContract', specificFolder))
%     xi_a = 1; 
%     xi_b = 1;
%     xi_l = 0.2;
% elseif(strcmpi('DiffContract', specificFolder))
    xi_a = 1.0; 
    xi_b = 1.0;
    xi_l = 0.5;
% else
%     xi_a = 0.2; 
%     xi_b = 1.8;
%     xi_l = 0.2;
% end

stretchedTissueTension =  xi_a + xi_b - K*xi_l./(domestrain + 1).^(3/2);

% reverse engineering the thing
% assuming P is correct
inferredRadius = 2*stretchedTissueTension./lumpressure;

% assuming that R is correct
inferredPressure = 2*stretchedTissueTension./radius;
correctedTension = tension*inferredPressure(end)/lumpressure(end);

disp("Saving workspace")
save('[pathtosimulationfolder]/workspace')
