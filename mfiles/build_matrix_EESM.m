% Copyright 2024
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

function rotore = build_matrix_EESM(temp,geo)
%%
for ii = 1:1:7
    x_var = sprintf('xp%d',ii);
    y_var = sprintf('yp%d',ii);
    eval(sprintf('%s = temp.(x_var);', x_var));
    eval(sprintf('%s = temp.(y_var);', y_var));
end
for ii = 1:1:4
    x_var = sprintf('xc%d',ii);
    y_var = sprintf('yc%d',ii);
    eval(sprintf('%s = temp.(x_var);', x_var));
    eval(sprintf('%s = temp.(y_var);', y_var));
end
xM = temp.xM;
yM = temp.yM;
narcs = temp.narcs;
theta = temp.theta;
if geo.r_fillet ~= 0
    xc_fillet = temp.xc_fillet;
    yc_fillet = temp.yc_fillet;
end
if geo.r_bfillet ~= 0
    xc_bfillet = temp.xc_bfillet;
    yc_bfillet = temp.yc_bfillet;
end
rotore = [];
Mag = [];
% Cu <- Mag || Mag <- Cu

materialCodes;
indexEle = 1;

% if geo.ps == 1
%     rotore = [rotore
%         %0    0    xp1     yp1     xp2  yp2   1 codMatFeRot indexEle
%         xp2  yp2  xp3     yp3     NaN  NaN   0 codMatFeRot indexEle
%         ];
% end
rotore = [rotore
    0    0    xp1     yp1    xp2  yp2   1 codMatFeRot indexEle
    %xp2  yp2  xp3     yp3    NaN  NaN   0 codMatFeRot indexEle
    %0    0    xp3     yp3     xp4  yp4  -1 codMatFeRot indexEle
    xp3  yp3  xp4(1)  yp4(1)  NaN  NaN   0 codMatFeRot indexEle
    ];
if geo.r_bfillet == 0
    rotore = [rotore
    xp4  yp4  xp5(1)  yp5(1)  NaN  NaN   0 codMatFeRot indexEle
    ];
else
    rotore = [rotore
    xc_bfillet yc_bfillet xp4(1) yp4(1) xp4(2) yp4(2) 1 codMatFeRot indexEle
    xp4(2) yp4(2) xp5(1)    yp5(1)    NaN    NaN    0 codMatFeRot indexEle
    ];
end

if geo.r_fillet == 0
    rotore = [rotore
    xp5  yp5  xp6     yp6     NaN  NaN   0 codMatFeRot indexEle
    ];
else
    Mag = [Mag
    xc_fillet yc_fillet xp5(1) yp5(1) xp5(2) yp5(2) 1 codMatFeRot indexEle
    xp5(2) yp5(2) xp6    yp6    NaN    NaN    0 codMatFeRot indexEle
    ];
end
rotore = [rotore
    
    xp6  yp6  xp7(1)  yp7(1)  NaN  NaN   0 codMatFeRot indexEle
    ];

for ii = 1:1:narcs-1
    if ii == narcs-1
        [xo(ii), yo(ii)] = calc_center_given_3pts(xp7(ii), yp7(ii), xp7(ii+1), yp7(ii+1), xM, yM);
    else
        [xo(ii), yo(ii)] = calc_center_given_3pts(xp7(ii), yp7(ii), xp7(ii+1), yp7(ii+1), xp7(ii+2), yp7(ii+2));
    end
    rotore = [rotore
        xo(ii) yo(ii) xp7(ii) yp7(ii) xp7(ii+1) yp7(ii+1) -1 codMatFeRot indexEle
        ];
end

[xoM, yoM] = calc_center_given_3pts(xp7(narcs), yp7(narcs), xM, yM, xp7(narcs), -yp7(narcs));
rotore = [rotore
    xoM yoM xp7(narcs) yp7(narcs) xM yM -1 codMatFeRot indexEle
    ];

indexEle = indexEle+1;

if geo.r_bfillet == 0 || length(xc1) == 1
    Mag = [Mag
        xc1 yc1 xc2(1) yc2(1) NaN NaN 0 codMatCuRot indexEle
        ];
else
    Mag = [Mag
        xc_bfillet, yc_bfillet xc1(1) yc1(1) xc1(2) yc1(2) 1 codMatCuRot indexEle
        xc1(2) yc1(2) xc2(1) yc2(1) NaN NaN 0 codMatCuRot indexEle
        ];
end
if geo.r_fillet == 0
    Mag = [Mag
    xc2 yc2 xc3 yc3 NaN NaN 0 codMatCuRot indexEle
    %xc2 yc2 xc3 yc3 NaN NaN 0 codMatFeRot indexEle
    ];
else
    Mag = [Mag
    xc_fillet yc_fillet xc2(1) yc2(1) xc2(2) yc2(2) 1 codMatCuRot indexEle
    xc2(2) yc2(2) xc3    yc3    NaN    NaN    0 codMatCuRot indexEle
    % xc_fillet yc_fillet xc2(1) yc2(1) xc2(2) yc2(2) 1 codMatFeRot indexEle
    % xc2(2) yc2(2) xc3    yc3    NaN    NaN    0 codMatFeRot indexEle
    ];
end
Mag = [Mag
    xc3 yc3 xc4 yc4 NaN NaN 0 codMatCuRot indexEle
    xc4 yc4 xc1(1) yc1(1) NaN NaN 0 codMatCuRot indexEle
    ];

rotore = [rotore;Mag];

% figure
% figSetting
% axis equal
% plot(xp7,yp7,'or')
% hold on
% plot(xM,yM,'or')

% test testa polo
% r = temp.rcirc;
% x_center = 0;
% y_center = 0;
% theta = linspace(0, theta, 100);
% x = r * cos(theta) + x_center;
% y = r * sin(theta) + y_center;
% hold on
% plot([x,0],[y,0], 'b-');
