% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [A,G] = calcAreaShape(Mat)

nEle = size(Mat,1);

nPointArc = 201;

X = [];
Y = [];

for ii=1:nEle
    if Mat(ii,7)==0 % line
        X = [X Mat(ii,1) Mat(ii,3)];
        Y = [Y Mat(ii,2) Mat(ii,4)];
    else % arc
        x0 = Mat(ii,1);
        y0 = Mat(ii,2);
        if Mat(ii,7)>0
            x1 = Mat(ii,3);
            y1 = Mat(ii,4);
            x2 = Mat(ii,5);
            y2 = Mat(ii,6);
        else
            x1 = Mat(ii,5);
            y1 = Mat(ii,6);
            x2 = Mat(ii,3);
            y2 = Mat(ii,4);
        end
        r = sqrt((x0 - x1)^2 + (y0 - y1)^2);
        ang1 = atan2((y1-y0),(x1-x0));
        ang2 = atan2((y2-y0),(x2-x0));
        if ang1 < 0
            ang1 = ang1 + 2*pi;
        end
        if ang2 < 0
            ang2 = ang2 + 2*pi;
        end
        theta = linspace(ang1,ang2,nPointArc);
        if (ang2 - ang1) < 0
            ang1 = ang1 - 2*pi;
            theta = linspace(ang1,ang2,201);
        end
        x_n = r*cos(theta(2:end-1)) + x0;
        y_n = r*sin(theta(2:end-1)) + y0;
        if ((X(end)==x2)&&(Y(end)==y2))
            X = [X x2 fliplr(x_n) x1];
            Y = [Y y2 fliplr(y_n) y1];
        else
            X = [X x1 x_n x2];
            Y = [Y y1 y_n y2];
        end
    end
end
warning('off')
shape = polyshape(X,Y);
warning('on')
% shape = convhull(shape);
A     = area(shape,1);
[x,y] = centroid(shape,1);
G = (x^2+y^2)^0.5;

