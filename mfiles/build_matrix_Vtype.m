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

%% Definizione matrice "rotore" per la costruzione della geometria di rotore della macchina
%% rev.Gallo 20/04/2018

function rotore = build_matrix_Vtype(temp,geo)

x0=geo.x0;

B1k=temp.B1k;
B2k=temp.B2k;
xpont=temp.xpont;
ypont=temp.ypont;
xxD1k=temp.xxD1k;
yyD1k=temp.yyD1k;
xxD2k=temp.xxD2k;
yyD2k=temp.yyD2k;
XcRibTraf1=temp.XcRibTraf1;
YcRibTraf1=temp.YcRibTraf1;
XcRibTraf2=temp.XcRibTraf2;
YcRibTraf2=temp.YcRibTraf2;

XpontRadSx=temp.XpontRadSx;
YpontRadSx=temp.YpontRadSx;
XpontRadDx=temp.XpontRadDx;
YpontRadDx=temp.YpontRadDx;
XpontRadBarDx=temp.XpontRadBarDx;
YpontRadBarDx=temp.YpontRadBarDx;
XpontRadBarSx=temp.XpontRadBarSx;
YpontRadBarSx=temp.YpontRadBarSx;

%Punti caratteristici aggiuntivi per geometria Vtype - rev.Gallo
XcRaccpontRadSx=temp.XcRaccpontRadSx; %raggio raccordi su ponticello radiale
YcRaccpontRadSx=temp.YcRaccpontRadSx;
XcRaccpontRadDx=temp.XcRaccpontRadDx;
YcRaccpontRadDx=temp.YcRaccpontRadDx; 

Xmag5dx=temp.XMag5;  % lato superiore area rettangolare magnete
Ymag5dx=temp.YMag5;
Xmag6sx=temp.XMag6;
Ymag6sx=temp.YMag6;

XMagpontRadSx=temp.XMagpontRadSx; %lato inferiore area rettangolare magnete
YMagpontRadSx=temp.YMagpontRadSx;
XMagpontRadDx=temp.XMagpontRadDx;
YMagpontRadDx=temp.YMagpontRadDx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (geo.BarFillFac==0)
    rotore=[];
    
    for ii=1:geo.nlay
        
        if (YpontRadSx(ii)~=0) % ponticello radiale
            rotore=[rotore;
                XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;
                XcRaccpontRadSx(ii),YcRaccpontRadSx(ii),XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii),YpontRadBarSx(ii),-1;
                XcRaccpontRadDx(ii),YcRaccpontRadDx(ii),XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii),YpontRadBarDx(ii),1];
        else
            rotore=[rotore;
                XpontRadBarSx(ii), YpontRadBarSx(ii), XpontRadBarDx(ii), YpontRadBarDx(ii), NaN, NaN, 0];
        end
        % archi al traferro
        rotore=[ rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
        % rette laterali barriere
        rotore=[rotore;
            xxD1k(ii), yyD1k(ii), XpontRadBarSx(ii), YpontRadBarSx(ii),NaN,NaN,0;
            xxD2k(ii), yyD2k(ii), XpontRadBarDx(ii), YpontRadBarDx(ii),NaN,NaN,0];
        
        if (YpontRadSx(ii)==0) % bisettrice polo in barriera
            rotore=[rotore;
                xxD1k(ii) yyD1k(ii) xxD1k(ii) yyD1k(ii) NaN NaN 0];
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    
    rotore=[];
    
    for ii=1:geo.nlay
        
        if (YpontRadSx(ii)~=0) % ponticello radiale
            rotore=[rotore;
                XpontRadSx(ii),YpontRadSx(ii),XpontRadDx(ii),YpontRadDx(ii),NaN,NaN,0;
                XcRaccpontRadSx(ii),YcRaccpontRadSx(ii),XpontRadSx(ii),YpontRadSx(ii),XpontRadBarSx(ii),YpontRadBarSx(ii),-1;
                XcRaccpontRadDx(ii),YcRaccpontRadDx(ii),XpontRadDx(ii),YpontRadDx(ii),XpontRadBarDx(ii),YpontRadBarDx(ii),1];
        else
            rotore=[rotore;
                XpontRadBarSx(ii), YpontRadBarSx(ii), XpontRadBarDx(ii), YpontRadBarDx(ii), NaN, NaN, 0];
        end
        % archi al traferro
        rotore=[rotore;
            XcRibTraf1(ii) YcRibTraf1(ii) xpont(ii) ypont(ii) xxD1k(ii) yyD1k(ii) 1;
            XcRibTraf2(ii) YcRibTraf2(ii) xxD2k(ii) yyD2k(ii) xpont(ii) ypont(ii) 1];
       
        %rette  laterali barriere
        rotore=[rotore;
            xxD1k(ii), yyD1k(ii), XpontRadBarSx(ii), YpontRadBarSx(ii),NaN,NaN,0;
            xxD2k(ii), yyD2k(ii), XpontRadBarDx(ii), YpontRadBarDx(ii),NaN,NaN,0]; 
        %rette area rettangolare magnete
        rotore=[rotore;
            Xmag6sx(ii), Ymag6sx(ii), Xmag5dx(ii), Ymag5dx(ii), NaN, NaN, 0;
            XMagpontRadSx(ii), YMagpontRadSx(ii), XMagpontRadDx(ii), YMagpontRadDx(ii), NaN, NaN, 0];

    end
end

