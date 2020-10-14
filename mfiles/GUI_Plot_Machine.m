% Copyright 2019
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

function GUI_Plot_Machine(h,Mat,debug)

if nargin==2
    debug=0;
end

[nrig,ncol] = size(Mat);
hold(h,'on')
grid(h,'on')

colorCode{1} = [0 0 0]; % Air - rotor
colorCode{2} = [0 0 0]; % Air - stator
colorCode{3} = [0 0 0]; % Slot conductor
colorCode{4} = [0 0 0]; % Fe stator
colorCode{5} = [0 0 0]; % Fe rotor
colorCode{6} = [1 0 0]; % PM
colorCode{7} = [0 0 0]; % Shaft

if ncol==7
    Mat=[Mat ones(nrig,2)];
end


for ii = 1 : nrig
    if debug
        keyboard
    end
    if Mat(ii,7) == 0
        
        % draw lines
        x1 = Mat(ii,3); y1 = Mat(ii,4);
        x2 = Mat(ii,1); y2 = Mat(ii,2);
%         grid on
        plot(h,[x1,x2],[y1,y2],'Linewidth',2,'Color',colorCode{Mat(ii,8)});
%         grid minor, axis equal
        
    elseif abs(Mat(ii,7)) == 1 % || Mat(ii,ncol) == -1
        % draw arcs
        dati = Mat(ii,:);
        % centro
        x0 = dati(1); y0 = dati(2);
        if dati(7) > 0
            % punto 1
            x1 = dati(3); y1 = dati(4);
            %  punto 2
            x2 = dati(5); y2 = dati(6);
        else
            % punto 1
            x2 = dati(3); y2 = dati(4);
            %  punto 2
            x1 = dati(5); y1 = dati(6);
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
        theta = linspace(ang1,ang2,201);
        if (ang2 - ang1) < 0
            ang1 = ang1 - 2*pi;
            theta = linspace(ang1,ang2,201);
        end
        x_n = r*cos(theta(2:end-1)) + x0;
        y_n = r*sin(theta(2:end-1)) + y0;
        x = [x1 x_n x2];
        y = [y1 y_n y2];
        %grid on
        plot(h,x,y,'Linewidth',2,'Color',colorCode{Mat(ii,8)});
%         grid minor, axis equal
        
    elseif abs(Mat(ii,7))==eps % || Mat(i,ncol)== -eps   % tracciamento linee per area inserimento magneti nella geometria Seg/ISeg
        x1=Mat(ii,1);
        x2=Mat(ii,3);
        y1=Mat(ii,2);
        y2=Mat(ii,4);
        plot(h,[x1,x2],[y1,y2],'Linewidth',2,'Color',colorCode{Mat(ii,8)});
    elseif abs(Mat(ii,7))==(1+eps)
        % draw arcs
        dati = Mat(ii,:);
        % centro
        x0 = dati(1); y0 = dati(2);
        if dati(7) > 0
            % punto 1
            x1 = dati(3); y1 = dati(4);
            %  punto 2
            x2 = dati(5); y2 = dati(6);
        else
            % punto 1
            x2 = dati(3); y2 = dati(4);
            %  punto 2
            x1 = dati(5); y1 = dati(6);
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
        theta = linspace(ang1,ang2,201);
        if (ang2 - ang1) < 0
            ang1 = ang1 - 2*pi;
            theta = linspace(ang1,ang2,201);
        end
        x_n = r*cos(theta(2:end-1)) + x0;
        y_n = r*sin(theta(2:end-1)) + y0;
        x = [x1 x_n x2];
        y = [y1 y_n y2];
        %grid on
        plot(h,x,y,'Linewidth',2,'Color',colorCode{Mat(ii,8)});
    end
end
%hold off % Perchè???



