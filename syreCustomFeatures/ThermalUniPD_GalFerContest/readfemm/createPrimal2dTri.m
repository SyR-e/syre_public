function primal = createPrimal2dTri(primal)

% CREATEPRIMAL2DTRI creates the struct of primal ctriangular mesh starting
% from the information provided by the mesh generator
%
% USE:
% primal = createPrimal2dTri(P,T,M)
%
% INPUT:
% 'primal': partially filled PRIMAL structure
%
% OUTPUT:
% 'primal': struct of primal mesh
%
% NOTE:
% Subfunction specialized for triangles called by CREATEPRIMAL2D
%
% VERSION:
% Date: 21.01.2011
% Copyright(C) 2011-2014: Fabio Freschi (fabio.freschi@polito.it)
%
% HISTORY:
% 26.07.2012: refurbished routine
% 22.08.2012: added sign to face_edge field
% 08.01.2013: node coordinates can have 3 dimensions
% 12.04.2013: changed name from CREATEPRIMAL2D
% 12.04.2013: nodes are not sorted
% 14.08.2013: fixed face_edge sign
% 14.08.2013: added VARARGIN
% 20.08.2013: added field SORT
% 22.12.2014: added EXPR as optional input
% 19.01.2015: fixed orientation bug for quadrangular patches
% 08.02.2016: new data structure

% fprintf('creating primal data structure... ');
% tic;

% element type
primal.ElementType = 'Tri';

% edge-to-node with repetition
e2n = reshape(primal.Fac2Nod(:,[1 2 1 3 2 3])',2,3*size(primal.Fac2Nod,1))';

% remove repetitions and get edges
[primal.Edg2Nod,~,f2e] = unique(e2n,'rows');
primal.Fac2Edg = reshape(int32(f2e),3,size(primal.Fac2Nod,1))';

% add sign to edges
primal.Fac2Edg(:,2) = -primal.Fac2Edg(:,2);

% edge orientation: from node with lower to node with higher index
primal.Edg2Nod(:,1) = -primal.Edg2Nod(:,1);

% fprintf('done %6.2f sec\n',toc);