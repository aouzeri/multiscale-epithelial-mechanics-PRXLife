% Copyright (c) Adam Ouzeri 2025

clc;
clear;

foldername = '[pathtosimulationfolder]/'; 
savingfoldername = foldername;

timefilename = sprintf('%sTimeData.txt',foldername);
timedata     = load(timefilename);
Nsteps       = length(timedata) + 1; % +1 for initial mesh

%% Read the mesh files
% Identify which elements belong to which face to evaluate data on specific
% cells and to extract further information from the similation data
nFaces = load(sprintf('%svertexmesh_info.txt',foldername));

faceType       = zeros(nFaces,1);
cellfaces      = cell(nFaces,1);
faceNvertices  = zeros(nFaces,1);
faceNelems     = zeros(nFaces,1);

for faceID = 1:nFaces
    fprintf(1, 'Reading facemesh %d\n', faceID);
    
    meshfilename = sprintf('%svertexmesh_%d.txt',foldername,faceID-1);
    fid = fopen(meshfilename,'r');
    if( fid==-1 )
        error('Can''t open the meshfile: %s',meshfilename);
    end
    
    str = fgets(fid);
    faceType(faceID) = sscanf(str,'%d', 1);
    
    str = fgets(fid);
    cellfaces{faceID} = sscanf(str,'%d')';
    
    str = fgets(fid);
    faceNvertices(faceID) = sscanf(str,'%d')';
    
    str = fgets(fid);
    faceNelems(faceID) = sscanf(str,'%d')';
    
    fclose(fid);
end
save([savingfoldername, 'facetype.mat'], 'faceType')

Ncells = nnz(faceType == 0);

% Identify which elements belong to which faces
elemsinface = cell(nFaces,1);
nodesinface = cell(nFaces,1);

elemcount = 0;
for faceID = 1:nFaces
    elemsinface{faceID} = elemcount+1:1:elemcount+faceNelems(faceID);
    elemcount = elemcount + faceNelems(faceID);
end
nElems = elemcount;

elemrhodata    = zeros(nElems, Nsteps);
elemtype       = zeros(nElems, 1);

faceareadata   = zeros(nFaces, Nsteps);
projectedfaceareadata = zeros(nFaces, Nsteps);
facerhodata    = zeros(nFaces, Nsteps);
faceLxdata     = zeros(nFaces, Nsteps);
faceLydata     = zeros(nFaces, Nsteps);
faceGdata      = zeros(nFaces, 6, Nsteps);
faceGfdata     = zeros(nFaces, 6, Nsteps);
facebarydata   = zeros(nFaces, 3); % Stores face barycenters

%% Read the vtk files and extract face data
havenodesinface = 0;
for step = 1:Nsteps
    fprintf(1,'Reading data for Step %d\n',step);
    vtkfilename = sprintf('%ssolution%d.vtk',foldername,step-2);
    [vertexR,connec,fieldnamesmat,fielddatamat] = read_vtk(vtkfilename);
    fielddatamat = cell2mat(fielddatamat);
    Nvertices    = length(vertexR);
    Nelems       = length(connec);
    
    nodeDOFs = fielddatamat(:,2:4);
    AUXDOFs  = fielddatamat(:,5:end);
    vertex_current = reshape(nodeDOFs,Nvertices,3); 
    
    % Identify which element belongs to which face
    elemtoface = zeros(Nelems, 1);
    for faceID = 1:nFaces
        elemIDSinface =  elemsinface{faceID};
        elemtoface(elemIDSinface) = faceID;
    end
    
    % Get elemental areas
    [elemArea, elemBary] = getTriangulationArea(vertex_current,connec);
    projectedvertex_current = vertex_current;
    projectedvertex_current(:,3) = 0;
    [projectedelemArea, ~] = getTriangulationArea(projectedvertex_current,connec);
    
    % Get nodes that belong to a given face (need to be done only once)
    if havenodesinface == 0
        for faceID = 1:nFaces
            elemIDSinface =  elemsinface{faceID};
            nodesinface{faceID} = unique(connec(elemIDSinface(:),1:3));
            elemtype(elemIDSinface) = faceType(faceID);
        end
        apicalelemIDs = elemtype == 0;
        havenodesinface = 1;
    end
    
    % Get face area
    for faceID = 1:nFaces
        elemIDSinface =  elemsinface{faceID};
        for ii = 1:length(elemIDSinface)
            elemID = elemIDSinface(ii);
            faceareadata(faceID,step) = faceareadata(faceID,step) ...
                + elemArea(elemID);
            projectedfaceareadata(faceID,step) = projectedfaceareadata(faceID,step) ...
                + projectedelemArea(elemID);
        end
    end
    
    % Assign data to the face from one of the elements
    % assuming that the constraints imply all nodes on a given face have the
    % same density values, so just take value from the first node. For G
    % and barycenter, the mean of the G value and elemental barycenters
	% is used
    for faceID = 1:nFaces
        elemIDSinface  = elemsinface{faceID};
        nodeIDsinface  = nodesinface{faceID};

		elemrhodata(elemIDSinface,step) = AUXDOFs(nodeIDsinface(1),1);
        facerhodata(faceID,step) = AUXDOFs(nodeIDsinface(1),1);
	% faceLxdata(faceID,step) = AUXDOFs(nodeIDsinface(1),13);
	% faceLydata(faceID,step) = AUXDOFs(nodeIDsinface(1),14); 
        faceGdata(faceID,:,step) = mean(AUXDOFs(nodeIDsinface,2:7));
	        
		% Assign barycenters
		facebarydata(faceID, :)  = mean(elemBary(elemIDSinface,:));
		if faceType(faceID)==0 % Offset apical points to make them visible
			facebarydata(faceID, 3) = facebarydata(faceID, 3) + 0.4;
		end
		
		if faceType(faceID)==1 % Offset basal points to make them visible
			facebarydata(faceID, 3) = facebarydata(faceID, 3) - 0.4;
		end
		
		if faceType(faceID)>1 % move lateral barycenter in the direction of the outward normal
			% get normal direction using any element on the lateral face
			elemID = elemIDSinface(1);
			verticesmat = connec(elemID,:);
			P1 = vertex_current(verticesmat(1),:);
			P2 = vertex_current(verticesmat(2),:);
			P3 = vertex_current(verticesmat(3),:);
			normadir = cross3d(P2 - P1,P3 - P1);
			normadir  = normadir/norm(normadir);
			% Trasnlate the lateral barycenter along this normal
			facebarydata(faceID, :) = facebarydata(faceID, :) + ...
				0.2 * normadir;
		end
    end
    
    % Save density data on apical elements for correlation analysis
%    elemBaryapical = elemBary(apicalelemIDs,:);
%    elemRhoapical  = elemrhodata(apicalelemIDs,step);
%    saveRhoapfilename = sprintf('%sapicalrho%d.txt',savingfoldername,step-1);
%    fp = fopen(saveRhoapfilename,'w+');
%    fprintf(fp,'%.4e %.4e %.4e %.4e\n', [elemBaryapical,elemRhoapical]');
%    fclose(fp);
%    fprintf(1,'====> Saving apical rho data in %s\n',saveRhoapfilename);
    
%    % Save the Gs in a seperate vtk using the barycentric vertices as nodes
%    saveGfilename = sprintf('%ssolutionGbary%d.vtk',savingfoldername,step-1);
%    fp = fopen(saveGfilename,'w+');
    
%    fprintf(fp,'# vtk DataFile Version 2.0\n');
%    fprintf(fp,'vtk output\n');
%    fprintf(fp,'ASCII\n');
%    fprintf(fp,'DATASET POLYDATA\n');
    
    % Write barycentric coordinates
%    fprintf(fp,'POINTS %d float\n',nFaces);
%    fprintf(fp,'%.4e %.4e %.4e\n', [facebarydata(:, 1), facebarydata(:, 2), facebarydata(:, 3)]');
    
    % Write Gsharp components as point data
%    fprintf(fp,'POINT_DATA %d\n',nFaces);
%    writeScalarPointData(fp,faceGdata(:,1,step),nFaces,'GsXX');
%    writeScalarPointData(fp,faceGdata(:,2,step),nFaces,'GsXY');
%    writeScalarPointData(fp,faceGdata(:,3,step),nFaces,'GsXZ');
%    writeScalarPointData(fp,faceGdata(:,4,step),nFaces,'GsYY');
%    writeScalarPointData(fp,faceGdata(:,5,step),nFaces,'GsYZ');
%    writeScalarPointData(fp,faceGdata(:,6,step),nFaces,'GsZZ');
	
	% Evaluate and write Gflat components as point data
	for faceID = 1:nFaces
		GsXX = faceGdata(faceID,1,step);
		GsXY = faceGdata(faceID,2,step);
		GsXZ = faceGdata(faceID,3,step);
		GsYY = faceGdata(faceID,4,step);
		GsYZ = faceGdata(faceID,5,step);
		GsZZ = faceGdata(faceID,6,step);
		Gsharp = [GsXX,GsXY,GsXZ; GsXY, GsYY, GsYZ; GsXZ, GsYZ, GsZZ];
		
		%Gflat  = inv(Gsharp);
		% Since Gsharp is rank-2, the inverse has to be taken by inverting
		% the relevant, non-zero eigenvalues and then transforming with the
		% eigenbasis again
		[V,D] = eig(Gsharp);
		Dinv = zeros(3,3);
		for i = 1:3
			if abs(D(i,i)) > 1e-10
				Dinv(i,i) = 1/D(i,i);
			end
		end
		Gflat = V * Dinv * V';
%	Gflat = Gsharp ; % I think at one point I  was calculating Gflat directly from hiperlife	
		faceGfdata(faceID,1,step) = Gflat(1,1);
		faceGfdata(faceID,2,step) = Gflat(1,2);
		faceGfdata(faceID,3,step) = Gflat(1,3);
		faceGfdata(faceID,4,step) = Gflat(2,2);
		faceGfdata(faceID,5,step) = Gflat(2,3);
		faceGfdata(faceID,6,step) = Gflat(3,3);
	end
	
%	writeScalarPointData(fp,faceGfdata(:,1,step),nFaces,'GfXX');
%     writeScalarPointData(fp,faceGfdata(:,2,step),nFaces,'GfXY');
%    writeScalarPointData(fp,faceGfdata(:,3,step),nFaces,'GfXZ');
%    writeScalarPointData(fp,faceGfdata(:,4,step),nFaces,'GfYY');
%    writeScalarPointData(fp,faceGfdata(:,5,step),nFaces,'GfYZ');
%    writeScalarPointData(fp,faceGfdata(:,6,step),nFaces,'GfZZ');
	
%    fclose(fp);
%    fprintf(1,'====> Saving Gdata data in %s\n',saveGfilename);
end

save([savingfoldername,'FaceArea.mat'], 'faceareadata')
save([savingfoldername,'ProjectedFaceArea.mat'], 'projectedfaceareadata')
save([savingfoldername, 'Gflat_per_face_and_step.mat' ], 'faceGfdata')
save([savingfoldername, 'Rho_per_face_and_step.mat'], 'facerhodata')
save([savingfoldername, 'FaceLx.mat'], 'faceLxdata' )
save([savingfoldername, 'FaceLy.mat'], 'faceLydata')
% Plotting mesh (color by numbers)
% patch('Faces',connec,'Vertices',vertex_current,'CData',elemtoface/Nelems,...
%     'FaceColor','flat');
% colormap(lines);

% Plotting mesh (color by cortical density)
% patch('Faces',connec,'Vertices',vertex_current,'CData',elemrhodata(:,step),...
%     'FaceColor','flat');
% colormap(winter);

% Plot barycenters
% hold on;
% plot3(facebarydata(:, 1), facebarydata(:, 2), facebarydata(:, 3),'ko');

