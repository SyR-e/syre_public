function primal = createPrimal2dQuad(primal)

% CREATEPRIMAL2DQUAD creates the struct of primal quadrilateral mesh
% starting from the information provided by the mesh generator
%
% USE:
% primal = createPrimal2dQuad(P,T,M)
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
% Date: 11.07.2016
%
% HISTORY:

fprintf('creating primal data structure... ');
tic;

% element type
primal.ElementType = 'Quad';

% edge-to-node with repetition
[T,IX] = sort(primal.Fac2Nod,2);
e2n = reshape(T(:,[1 2 1 3 2 4 3 4])',2,4*size(primal.Fac2Nod,1))';

% remove repetitions and get edges
[primal.Edg2Nod,~,f2e] = unique(e2n,'rows');
primal.Fac2Edg = reshape(int32(f2e),4,size(primal.Fac2Nod,1))';

% add sign to edges
iSign = -ones(size(primal.Fac2Nod),'int32');
% 1st edge
idx = (IX(:,1) == 1 & IX(:,2) == 2) | (IX(:,1) == 2 & IX(:,2) == 3) | (IX(:,1) == 3 & IX(:,2) == 4) | (IX(:,1) == 4 & IX(:,2) == 1);
iSign(idx,1) = 1;
% 2nd edge
idx = (IX(:,1) == 1 & IX(:,3) == 2) | (IX(:,1) == 2 & IX(:,3) == 3) | (IX(:,1) == 3 & IX(:,3) == 4) | (IX(:,1) == 4 & IX(:,3) == 1);
iSign(idx,2) = 1;
% 3rd edge
idx = (IX(:,2) == 1 & IX(:,4) == 2) | (IX(:,2) == 2 & IX(:,4) == 3) | (IX(:,2) == 3 & IX(:,4) == 4) | (IX(:,2) == 4 & IX(:,4) == 1);
iSign(idx,3) = 1;
% 3rd edge
idx = (IX(:,3) == 1 & IX(:,4) == 2) | (IX(:,3) == 2 & IX(:,4) == 3) | (IX(:,3) == 3 & IX(:,4) == 4) | (IX(:,3) == 4 & IX(:,4) == 1);
iSign(idx,4) = 1;
primal.Fac2Edg = primal.Fac2Edg.*iSign;

% edge orientation: from node with lower to node with higher index
primal.Edg2Nod(:,1) = -primal.Edg2Nod(:,1);

fprintf('done %6.2f sec\n',toc);