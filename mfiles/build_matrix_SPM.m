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

function rotore = build_matrix_SPM(temp,geo)

xPMco = temp.xPMco;
yPMco = temp.yPMco;
xPMo  = temp.xPMo;
yPMo  = temp.yPMo;
xPMci = temp.xPMci;
yPMci = temp.yPMci;
xPMi  = temp.xPMi;
yPMi  = temp.yPMi;

xArccenter = temp.xArccenter;
yArccenter = temp.yArccenter;

xAiri = temp.xAiri;
yAiri = temp.yAiri;
xAiro = temp.xAiro;
yAiro = temp.yAiro;

xPMso = temp.xPMso;
yPMso = temp.yPMso;
xPMsi = temp.xPMsi;
yPMsi = temp.yPMsi;


rotore = [];
Mag = [];

materialCodes;
% codMatAirRot = 1;
% codMatBar  = 6;

indexEle = 1;

rotore = [rotore
    xPMci yPMci xPMco yPMco NaN   NaN    0 codMatAirRot indexEle
    0     0     xPMco yPMco xAiro yAiro +1 codMatAirRot indexEle
    xAiro yAiro xAiri yAiri NaN   NaN    0 codMatAirRot indexEle
    0     0     xAiri yAiri xPMci yPMci -1 codMatAirRot indexEle
    ];
indexEle = indexEle+1;

Mag = [Mag
    xPMci      yPMci      xPMco yPMco NaN   NaN    0 codMatBar indexEle
    xArccenter yArccenter xPMco yPMco xPMo  yPMo  +1 codMatBar indexEle
    xPMo       yPMo       xPMi  yPMi  NaN   NaN    0 codMatBar indexEle
    0          0          xPMi  yPMi  xPMci yPMci -1 codMatBar indexEle
    ];
indexEle = indexEle+1;

if ~isempty(xPMsi)
    for ii=1:length(xPMsi)
        Mag = [Mag
            xPMsi(ii) yPMsi(ii) xPMso(ii) yPMso(ii) NaN NaN 0 codMatBar -1
            ];
    end
end

rotore = [rotore;Mag];

