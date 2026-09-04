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

function [optRes] = ORE_loadData(pathname,filename)

if nargin()==0
    [filename,pathname] = uigetfile('*.mat','Load optimization results');
end

load(checkPathSyntax([pathname filename]))

[dataSet,~,~,~,~] = back_compatibility(dataSet,geo0,per,0);
[bounds,objs,~,~,~] = data0(dataSet);

optRes.varNames       = OUT.Param.geo0.RQnames;
optRes.varBounds      = bounds;
optRes.varParetoFront = front(:,1:length(optRes.varNames));

optRes.objNames  = OUT.Param.geo0.OBJnames(objs(:,2)==1);
optRes.objLimits = objs(objs(:,2)==1,1);
optRes.objParetoFront = front(:,length(optRes.varNames)+1:end);

index = 1:1:numel(idx);
optRes.motNum = index(idx==1);

optRes.filt = ones(size(optRes.motNum));

optRes.dataSet = dataSet;
optRes.geo     = geo0;
optRes.per     = per;
optRes.mat     = mat;


% unpenalize Pareto Front

for ii=1:length(optRes.objLimits)
    limit = optRes.objLimits(ii);
    if limit<0
        optRes.objParetoFront(optRes.objParetoFront(:,ii)>limit/10,ii) = optRes.objParetoFront(optRes.objParetoFront(:,ii)>limit,ii)*10;
    else
        optRes.objParetoFront(optRes.objParetoFront(:,ii)>limit*10,ii) = optRes.objParetoFront(optRes.objParetoFront(:,ii)>limit,ii)/10;
    end
end


