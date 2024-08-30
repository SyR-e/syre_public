function primal = createPrimal2d(P,T,M,varargin)

% CREATEPRIMAL2D creates the struct of primal complex of cells starting from
% the information provided by the 3D mesh generator
%
% USE:
% primal = createprimal2D(P,T,M)
%
% INPUT:
% 'P': node coordinates
% 'T': tetrahedra/hexahedra connectivity matrix
% 'M': vector with material codes
% 'varargin': variable inputs
%   * 'scale': scaling factor (default: 1)
%
% OUTPUT:
% 'primal': struct of primal mesh
%
% VERSION:
% Date: 19.12.2009
% Copyright(C) 2009-2018: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 14.01.2010: fixed an incongruence on primal.face_edge by flipping 2 and 3 columns
% 11.01.2010: fixed bug on face area orientation
% 13.01.2010: fixed bug on volume calculation for hexahedra
% 20.01.2011: added SCALE as optional input
% 26.05.2011: removed unused nodes
% 17.08.2011: removed calculation of edge lengths, face areas and cell volumes
% 20.09.2011: new refurbished routines CREATEPRIMAL3D_XXX
% 24.05.2012: calculation of incidence matrices (backward compatible)
% 03.01.2012: added input check
% 21.12.2015: added SYM field
% 15.01.2016: new data structure
% 31.10.2018: added mesh info

% 14.12.2015: added re-orientation of PRIMAL.FxE according to PRIMAL.SORT

% set defaults
scale = 1;

% check varargin
for i = 1:2:length(varargin)
    vin = varargin{i};
    if strcmpi(vin,'scale')
        scale = varargin{i+1};
    end
end

% input consistency check
if size(P,2) ~= 3
    P = [P zeros(size(P,1),1)];
end
if length(M) ~= size(T,1)
    error('DualLab:createprimal2d','Size mismatch between T and M');
end

% define structure array
primal = struct(...
    'Node',[],...        % data from mesh generator
    'Fac2Nod',[],...     % connectivity
    'Object',[],...         % object codes
    'Symmetry',[0 0 0],...    % symmetry info
    'Scale',scale,...    % scaling factor
    'ElementType',[],... % element type
    'Sort',[],...        % sorting index
    'Edg2Nod',[],...     % derived structures
    'Fac2Edg',[]);

% fix connectivity with 0 node
if min(T(:)) == 0
    T = T+1;
end

% removed unused nodes
[pix,~,jx] = unique(T);
T = reshape(jx,size(T));
P = P(pix,:);

% nodes
primal.Node = P;
% material code
primal.Object = int32(M);

if size(T,2) == 3
    % soted connectivity 
    [primal.Fac2Nod,IX] = sort(int32(T),2);

    % find sorted elements
    primal.Sort = -ones(size(T,1),1,'int32');
    idx = (IX(:,1) == 1 & IX(:,2) == 2 & IX(:,3) == 3) | ...
        (IX(:,1) == 2 & IX(:,2) == 3 & IX(:,3) == 1) | ...
        (IX(:,1) == 3 & IX(:,2) == 1 & IX(:,3) == 2);
    primal.Sort(idx) = 1;

    primal = createPrimal2dTri(primal);
else
    primal.Fac2Nod = int32(T);
    primal.Sort = ones(size(T,1),1,'int32');

    primal = createPrimal2dQuad(primal);
end

% evaluate mesh quality
q = meshQuality2d(primal);

% print mesh stats
% fprintf('---------------------------\n');
% fprintf('Mesh info\n');
% fprintf('---------------------------\n');
% fprintf('Faces: %d\n',size(primal.Fac2Nod,1));
% fprintf('Edges: %d\n',size(primal.Edg2Nod,1));
% fprintf('Nodes: %d\n',size(primal.Node,1));
% fprintf('Objects: %d\n',length(unique(primal.Object)));
% fprintf('min q2: %e\n',min(q));
% fprintf('max q2: %e\n',max(q));
% fprintf('ave q2: %e\n',mean(q));
% fprintf('---------------------------\n');


end

function q = meshQuality2d(primal)

if strcmpi(primal.ElementType,'Tri')
    q = triMeshQuality(primal.Node,primal.Fac2Nod);
else
    q = 0;
end


end