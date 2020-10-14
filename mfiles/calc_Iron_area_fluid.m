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

%% Calcolo delle aree di ferro di barriera

A=zeros(1,nlay);
xBarrier1=[];
yBarrier1=[];
xBarrier2=[];
yBarrier2=[];


for i=1:nlay
       Dangle=0.01;

    if (i==1)
        [xc2,yc2,r2,angleA2,angleB2]=circonferenza_per_3_pti(B2k(i),0,xTraf02(i),yTraf02(i),xxB2k_mean(i),yyB2k_mean(i));
        if (angleA2<angleB2)
            angleB2=-2*pi+angleB2;
        end
        
        xBarrier2{i}=xc2+r2*cos([angleB2:Dangle:angleA2]);
        yBarrier2 {i}=yc2+r2*sin([angleB2:Dangle:angleA2]);
        
        [xc1,yc1,r1,angleA1,angleB1]=circonferenza_per_3_pti(B1k(i),0,xTraf01(i),yTraf01(i),xxB1k_mean(i),yyB1k_mean(i));
        
        if (angleA1<angleB1)
            angleB1=-2*pi+angleB1;
        end
        xBarrier1{i}=xc1+r1*cos([angleB1:Dangle:angleA1]);
        yBarrier1{i}=yc1+r1*sin([angleB1:Dangle:angleA1]);

Dx=(r-xpont(1))/5;  %per avere almeno 5 divisioni;
xcir_plot=[r:-Dx:r*cos(pi/2/p)];
ycir_plot=sqrt(r^2-xcir_plot.^2);
VectCir=find(xcir_plot>=xpont(1));
x_ext_rot=xcir_plot(VectCir);
y_ext_rot=ycir_plot(VectCir);

        X=[xBarrier2{i},x_ext_rot];
        Y=[yBarrier2{i},y_ext_rot];

%         figure(100);hold on;fill(X,Y,'r');hold off;
        A(i)=polyarea(X,Y);
clear X Y;
        %%
    else
        [xc1,yc1,r1,angleA1,angleB1]=circonferenza_per_3_pti(B1k(i-1),0,xTraf01(i-1),yTraf01(i-1),xxB1k_mean(i-1),yyB1k_mean(i-1));
        [xc2,yc2,r2,angleA2,angleB2]=circonferenza_per_3_pti(B2k(i),0,xTraf02(i),yTraf02(i),xxB2k_mean(i),yyB2k_mean(i));
        
        if (angleA1<angleB1)
            angleB1=-2*pi+angleB1;
        end
        xBarrier1{i}=xc1+r1*cos([angleB1:Dangle:angleA1]);
        yBarrier1{i}=yc1+r1*sin([angleB1:Dangle:angleA1]);
                
        if (angleA2<angleB2)
            angleB2=-2*pi+angleB2;
        end
        xBarrier2{i}=xc2+r2*cos([angleB2:Dangle:angleA2]);
        yBarrier2{i}=yc2+r2*sin([angleB2:Dangle:angleA2]);

        X=[fliplr(xBarrier1{i}),xBarrier2{i}];
        Y=[fliplr(yBarrier1{i}),yBarrier2{i}];

%         figure(100);hold on;fill(X,Y,'r');hold off;
        A(i)=polyarea(X,Y);
        clear X Y;
    end
    
end
%% Last flux barrier it'snt previosly calculated
[xc1,yc1,r1,angleA1,angleB1]=circonferenza_per_3_pti(B1k(nlay),0,xTraf01(nlay),yTraf01(nlay),xxB1k_mean(nlay),yyB1k_mean(nlay));

if (angleA1<angleB1)
    angleB1=-2*pi+angleB1;
end
xBarrier1{nlay+1}=xc1+r1*cos([angleB1:Dangle:angleA1]);
yBarrier1{nlay+1}=yc1+r1*sin([angleB1:Dangle:angleA1]);

%% Total area calculation: 
Afe=cumsum(A);  % area of semipole

rG=([r B1k(1:nlay-1)]+B2k)/2;
M_Fe = 2*Afe*l * 1e-9 * rhoFE ;   % massa ferro appeso ai ponticelli

F_centrifuga = M_Fe .* rG/1000 *  (nmax * pi/30)^2;
%sigma_max = sigma_max;    % N/mm2 - snervamento lamierino

if geo.radial_ribs_eval == 0
    pont = F_centrifuga/(sigma_max * l);    % mm
else
    pont = geo.pontR;
end

pont(pont < pont0) = 0; % NOTA BENE: Elimino i ponticelli troppo sottili
hpont=pont;
racc_pont=abs(B1k-B2k)/4;
racc_pont=pont0*ones(size(pont));

for ii=1:nlay
    
    if hpont(ii)>0
        
        YpontRadBarSx(ii)=hpont(ii)+racc_pont(ii);
        XpontRadBarSx(ii)=interp1(yBarrier1{ii+1},xBarrier1{ii+1},YpontRadBarSx(ii));
        YpontRadBarDx(ii)=hpont(ii)+racc_pont(ii);
        XpontRadBarDx(ii)=interp1(yBarrier2{ii},xBarrier2{ii},YpontRadBarDx(ii));
        XpontRadDx(ii)=B2k(ii)-racc_pont(ii);
        YpontRadDx(ii)=hpont(ii);
        XpontRadSx(ii)=B1k(ii)+racc_pont(ii);
        YpontRadSx(ii)=hpont(ii);

    else
        
        
        XpontRadBarSx(ii)=B1k(ii);
        YpontRadBarSx(ii)=0;
        XpontRadBarDx(ii)=B2k(ii);
        YpontRadBarDx(ii)=0;
        XpontRadDx(ii)=NaN;
        YpontRadDx(ii)=0;
        XpontRadSx(ii)=NaN;
        YpontRadSx(ii)=0;
    end
    
end

