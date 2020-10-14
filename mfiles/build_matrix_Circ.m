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

function rotore = build_matrix_Circ(temp,geo)

x0=geo.x0;

XBan1dx=temp.XBanqdx;
XBan1sx=temp.XBanqsx;
xc=temp.xc;
yc=temp.yc;
XBan3dx=temp.X3;
YBan3dx=temp.Y3;
XBan3sx=temp.X4;
YBan3sx=temp.Y4;

XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;

YBan1dx=temp.YpontRadBarDx;
YBan1sx=temp.YpontRadBarSx;
error_mex=temp.error_mex;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (geo.BarFillFac==0)
    rotore=[];
    
    for ii=1:geo.nlay
        
        if (YpontRadSx(ii)~=0)
            rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
                XpontRadSx(ii),YpontRadSx(ii),XBan1sx(ii),YBan1sx(ii),NaN,NaN,0;...
                XpontRadDx(ii),YpontRadDx(ii),XBan1dx(ii),YBan1dx(ii),NaN,NaN,0];
        end
        if(error_mex(ii)==0)
            rotore=[ rotore;
                xc(ii) yc(ii) XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1;
                x0 0 XBan3sx(ii) YBan3sx(ii) XBan1sx(ii) YBan1sx(ii) 1];
            rotore=[rotore;
                x0 0 XBan3dx(ii) YBan3dx(ii) XBan1dx(ii) YBan1dx(ii) 1];
        else
            rotore=[rotore;
                xc(ii) 0 XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1];
            
        end
        
        if (YpontRadSx(ii)==0)
            rotore=[rotore;
                XBan1sx(ii) YBan1sx(ii) XBan1dx(ii) YBan1dx(ii) NaN NaN 0];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Xmag5dx=temp.X5;
    Ymag5dx=temp.Y5;
    Xmag6sx=temp.X6;
    Ymag6sx=temp.Y6;
    rotore=[];
    
    for ii=1:geo.nlay
        
        if (YpontRadSx(ii)~=0)
            rotore=[rotore;XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;...
                Xmag6sx(ii) Ymag6sx(ii) Xmag5dx(ii) Ymag5dx(ii) NaN NaN 0;...
                XpontRadSx(ii),YpontRadSx(ii),XBan1sx(ii),YBan1sx(ii),NaN,NaN,0;...
                XpontRadDx(ii),YpontRadDx(ii),XBan1dx(ii),YBan1dx(ii),NaN,NaN,0];
        end
        if(error_mex(ii)==0)
            rotore =[ rotore;
                xc(ii) yc(ii) XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1;
                x0 0 XBan3sx(ii) YBan3sx(ii) Xmag6sx(ii) Ymag6sx(ii) 1;
                x0 0 Xmag6sx(ii) Ymag6sx(ii)  XBan1sx(ii) YBan1sx(ii) 1];
            
            rotore=[rotore;
                x0 0 XBan3dx(ii) YBan3dx(ii) Xmag5dx(ii) Ymag5dx(ii) 1;
                x0 0 Xmag5dx(ii) Ymag5dx(ii) XBan1dx(ii) YBan1dx(ii) 1];
            
        else
            rotore=[rotore;
                xc(ii) 0 XBan3dx(ii) YBan3dx(ii) XBan3sx(ii) YBan3sx(ii) 1];
            
        end
        
        if (YpontRadSx(ii)==0)
            rotore=[rotore;
                XBan1sx(ii) YBan1sx(ii) XBan1dx(ii) YBan1dx(ii) NaN NaN 0;
                Xmag6sx(ii) Ymag6sx(ii) Xmag5dx(ii) Ymag5dx(ii) NaN NaN 0];
        end
    end
end

