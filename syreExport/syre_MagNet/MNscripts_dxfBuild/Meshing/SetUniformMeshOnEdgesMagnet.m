function SetUniformMeshOnEdgesMagnet(magData,edgesList,numOfSeg)

% SetUniformMeshOnEdgesMagnet.m [v1.00.00 (30-11-2012)]
% Sets a uniform mesh distribution on a list of edges
% =========================================================================
% Syntax: SetUniformMeshOnEdgesMagnet(magData,edgesList,numOfSeg)
% Input:
%         - magData:    Magnet's data structure
%         - edgeList:   list of edges {face name,edge number}
%         - numOfSeg:   number of segments
% =========================================================================

mh = magData.magnetHandler;

% Creates a visual basic array of elements like:
% ArrayOfValues(0)= "BACK IRON,Face#2,Edge#3"
% ArrayOfValues(1)= "AIRGAP ROTOR,Face#2,Edge#1"
% Those elements are needed to identify edges

nEdges = length(edgesList);

% Selection of involved edges
invoke(mh, 'processCommand',['ReDim ArrayOfValues(',num2str(nEdges-1),')']);
for k = 1 : nEdges
    objName    = edgesList{k}{1};
    faceNumber = edgesList{k}{2};
    edgeNumber = edgesList{k}{3};
    
    faceString = ['Face#',num2str(faceNumber)];
    edgeString = ['Edge#',num2str(edgeNumber)];
    
    finalString = ['"',objName,',',faceString,',',edgeString,'"'];
    invoke(mh, 'processCommand',['ArrayOfValues(',num2str(k-1),')=',finalString]);
end

% Set mesh properties
invoke(mh,'processCommand',['Call getDocument().assignMeshEdgeSubdivisions(ArrayOfValues, "Type=Uniform;Subdivisions=',num2str(numOfSeg),';DensityRatio=0.5")']);
