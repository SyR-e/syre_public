% Copyright 2014
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

function BarCenter = defineBlockCenters(temp,fem,geo)

% coordinates of center points of all FEMM blocks,
% where block labels will be eventually placed

res=fem.res;
p = geo.p;

% Material codes: refer to BLKLABELS.materials positions
materialCodes;  % load material codes

xc   = temp.xc;
yc   = temp.yc(~isnan(xc));
xmag = temp.xmag(~isnan(xc(1,:)));
ymag = temp.ymag(~isnan(xc(1,:)));
zmag = temp.zmag(~isnan(xc(1,:)));
xc   = xc(~isnan(xc));

if isfield(temp,'xair')
    xair = temp.xair;
    yair = temp.yair(~isnan(xair));
    xair = xair(~isnan(xair));
else
    xair = [];
    yair = [];
end

if ~isempty(xair)
    xc   = [xc';xc'];
    yc   = [yc';-yc'];
    xair = [xair';xair'];
    yair = [yair';-yair'];
    xmag = [xmag';xmag'];
    ymag = [ymag';-ymag'];
    zmag = [zmag';zmag'];
    BarCenter = [
        xc,yc,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag;
        xair,yair,codMatAirRot*ones(length(xair),1),res*ones(length(xair),1),1*ones(length(xair),1),zeros(length(xair),1),zeros(length(xair),1),zeros(length(xair),1)
        ];
    
else
    xc   = [xc';xc'];
    yc   = [yc';-yc'];
    xmag = [xmag';xmag'];
    ymag = [ymag';-ymag'];
    zmag = [zmag';zmag'];
    BarCenter = [xc,yc,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag];
end
