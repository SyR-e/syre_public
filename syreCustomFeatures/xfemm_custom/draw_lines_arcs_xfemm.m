% Copyright 2026
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function FemmProblem = draw_lines_arcs_xfemm(FemmProblem,Mat,group,res)

tol = 1e-6;

debug=0;

[nrig,~] = size(Mat);

nodes = [];
% get all the nodes
for ii=1:nrig
    if Mat(ii,7)==0
        if isempty(nodes)
            nodes = [Mat(ii,1),Mat(ii,2)];
        else
            xy = [Mat(ii,3) Mat(ii,4)];
            dist = (sum(nodes-repmat(xy,[size(nodes,1),1]),2).^2).^0.5;
            if min(dist)>tol
                nodes = [nodes;xy];
            end
        end
        xy = Mat(ii,3:4);
        dist = (sum(nodes-repmat(xy,[size(nodes,1),1]),2).^2).^0.5;
        if min(dist)>tol
            nodes = [nodes;xy];
        end
    else
        if isempty(nodes)
            nodes = [Mat(ii,3), Mat(ii,4)];
        else
            xy = Mat(ii,3:4);
            dist = (sum(nodes-repmat(xy,[size(nodes,1),1]),2).^2).^0.5;
            if min(dist)>tol
                nodes = [nodes;xy];
            end
        end
        xy = Mat(ii,5:6);
        dist = (sum(nodes-repmat(xy,[size(nodes,1),1])).^2).^0.5;
        if min(dist)>tol
            nodes = [nodes;xy];
        end
    end
end

FemmProblem = addnodes_mfemm(FemmProblem,nodes(:,1),nodes(:,2),'InGroup',group);

nodesNum = numel(nodes);

for ii=1:nrig
    if Mat(ii,7)==0 %line
        xy1 = Mat(ii,1:2);
        xy2 = Mat(ii,3:4);
        n0 = findnode_mfemm(FemmProblem,xy1);
        n1 = findnode_mfemm(FemmProblem,xy2);
        if debug
            [a,b,c] = retta_per_2pti(xy1(1),xy1(2),xy2(1),xy2(2));
            d=nan(1,nodesNum);
            for jj=1:nodesNum
                d(jj) = calc_distanza_retta_punto(a,b,c,nodes(jj,1),nodes(jj,2));
            end
            flagInt = d<=tol;
            nodesInt = nodes(flagInt==1,:);
            d1 = sum(nodesInt-repmat(xy1,[size(nodesInt,1),1]),2).^0.5;
            [~,index] = sort(d1);
            nodesInt = nodesInt(index,:);
            d2 = (sum(xy1-xy2).^2)^0.5;
            
        end

        FemmProblem = addsegments_mfemm(FemmProblem,n0,n1,'InGroup',group,'MaxSideLength',res);
    else % arc
        [maxsegdeg,raggio,ang1,ang,xy]=Disegna_Arco(Mat(ii,:),res);
        xy1 = xy(1,:);
        xy2 = xy(2,:);
        n0 = findnode_mfemm(FemmProblem,xy1);
        n1 = findnode_mfemm(FemmProblem,xy2);
        FemmProblem = addarcsegments_mfemm(FemmProblem,n0,n1,ang,'InGroup',group,'MaxSegDegrees',maxsegdeg);
    end
end



