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

function [dqtMap,fdfq] = skew_dqt(hEdit,motorModel)

debug=0;

% Load data
axisType  = motorModel.data.axisType;
motorType = motorModel.data.motorType;
p         = motorModel.data.p;
Lld       = motorModel.data.Lld;
Llq       = motorModel.data.Llq;

dqtMap    = motorModel.dqtMap;

ang_sk_m = motorModel.skew.thSkw;
nSlice   = motorModel.skew.nSlice;
nPoints  = motorModel.skew.nPoints;

ang_sk = ang_sk_m*p*pi/180; % elt rad
k = 1:1:nSlice;
k = k-mean(k);
alfa_k = k*ang_sk/(nSlice);

% k = 0:1:nSlice-1;
% alfa_k = (1/2+ k ) * ang_sk / nSlice - ang_sk / 2;

rot_alfa_k = exp(-1i*alfa_k);

IdMax = max(dqtMap.data.Id,[],'all');
IdMin = min(dqtMap.data.Id,[],'all');
IqMax = max(dqtMap.data.Iq,[],'all');
IqMin = min(dqtMap.data.Iq,[],'all');

% max current amplitude (diagonal of the id iq rectangle)
% diagonal_raw = abs(IdMax + 1j*IqMax);
diagonal_angle = angle(IdMax + 1j*IqMax);
% reduced diagonal (tentative)
diagonal_temp = IdMax/(cos(diagonal_angle-ang_sk/2));
id_span_new = diagonal_temp * cos(diagonal_angle);
iq_span_new = diagonal_temp * sin(diagonal_angle);
% if no good, try the other side
if (iq_span_new > IqMax)
    diagonal_temp = IqMax/(cos(diagonal_angle-ang_sk/2));
    id_span_new = diagonal_temp * cos(diagonal_angle);
    iq_span_new = diagonal_temp * sin(diagonal_angle);
end

% 3D matrixes of skewed machine
if IdMin==-IdMax
    id = linspace(-id_span_new,id_span_new,nPoints);
else
    id = linspace(IdMin,id_span_new,nPoints);
end
if IqMin==-IqMax
    iq = linspace(-iq_span_new,iq_span_new,nPoints);
else
    iq = linspace(IqMin,iq_span_new,nPoints);
end
th = dqtMap.th;
[data.Id,data.Iq,data.th]=ndgrid(id,iq,th);
[dS,qS,tS]=size(data.Id);
data.Fd = zeros(dS,qS,tS);
data.Fq = zeros(dS,qS,tS);
data.T  = zeros(dS,qS,tS);

fInt = dqtMap.fInt;
% The 3-d inductances must be removed from dqtMap.data and then fInt must
% be recomputed

data0.th = dqtMap.data.th;
data0.Id = dqtMap.data.Id;
data0.Iq = dqtMap.data.Iq;
data0.Fd = dqtMap.data.Fd-Lld*data0.Id;
data0.Fq = dqtMap.data.Fq-Llq*data0.Iq;
data0.T  = dqtMap.data.T;

fInt.Fd = griddedInterpolant(data0.Id,data0.Iq,data0.th,data0.Fd,'spline');
fInt.Fq = griddedInterpolant(data0.Id,data0.Iq,data0.th,data0.Fq,'spline');
fInt.T  = griddedInterpolant(data0.Id,data0.Iq,data0.th,data0.T,'spline');


for dd=1:dS
    for qq=1:qS
        FdTmp = zeros(1,1,tS);
        FqTmp = zeros(1,1,tS);
        TTmp  = zeros(1,1,tS);
        thTmp = data.th(dd,qq,:);
        if debug
            figure()
            figSetting
            title(['dd=' int2str(dd) ' / qq=' int2str(qq)])
            xlabel('$\theta_e$ [$^\circ$]')
            ylabel('$T$ [Nm]')
        end
        for ss=1:length(alfa_k)
            % position
            thSlice=thTmp+alfa_k(ss)*180/pi;
            thSlice(thSlice<0)=thSlice(thSlice<0)+360;
            thSlice(thSlice>360)=thSlice(thSlice>360)-360;
            % current
            IdSlice = data.Id(dd,qq,:);
            IqSlice = data.Iq(dd,qq,:);
            Idq     = (IdSlice+j*IqSlice)*rot_alfa_k(ss);
            IdSlice = real(Idq);
            IqSlice = imag(Idq);
            % flux and torque interpolation
            FdSlice = fInt.Fd(IdSlice,IqSlice,thSlice);
            FqSlice = fInt.Fq(IdSlice,IqSlice,thSlice);
            TSlice  = fInt.T(IdSlice,IqSlice,thSlice);
            
            % controllo se sono fuori dalle mappe originali (solo se
            % mancano i quadranti)
            
            thSliceInv=-thSlice;
            %thSliceInv=-thTmp+alfa_k(ss)*180/pi;
            thSliceInv(thSliceInv<0)=thSliceInv(thSliceInv<0)+360;
            thSliceInv(thSliceInv>360)=thSliceInv(thSliceInv>360)-360;
            thSliceInv(thSliceInv<0)=thSliceInv(thSliceInv<0)+360;
            thSliceInv(thSliceInv>360)=thSliceInv(thSliceInv>360)-360;
            
            if strcmp(motorType,'SR')
                % symmetry on both axis
                if (abs(IdMin)<IdMax)
                    FdSlice(real(Idq)<0)=-fInt.Fd(-IdSlice(real(Idq)<0),IqSlice(real(Idq)<0),thSliceInv(real(Idq)<0));
                    FqSlice(real(Idq)<0)=+fInt.Fq(-IdSlice(real(Idq)<0),IqSlice(real(Idq)<0),thSliceInv(real(Idq)<0));
                    TSlice(real(Idq)<0)=-fInt.T(-IdSlice(real(Idq)<0),IqSlice(real(Idq)<0),thSliceInv(real(Idq)<0));
                end
                if (abs(IqMin)<IqMax)
                    FdSlice(imag(Idq)<0)=+fInt.Fd(IdSlice(imag(Idq)<0),-IqSlice(imag(Idq)<0),thSliceInv(imag(Idq)<0));
                    FqSlice(imag(Idq)<0)=-fInt.Fq(IdSlice(imag(Idq)<0),-IqSlice(imag(Idq)<0),thSliceInv(imag(Idq)<0));
                    TSlice(imag(Idq)<0)=-fInt.T(IdSlice(imag(Idq)<0),-IqSlice(imag(Idq)<0),thSliceInv(imag(Idq)<0));
                end
            elseif strcmp(motorType,'PM')
                if strcmp(axisType,'SR')
                    %symmetry just on d axis
                    if (abs(IdMin)<IdMax)
                        FdSlice(real(Idq)<0)=-fInt.Fd(-IdSlice(real(Idq)<0),IqSlice(real(Idq)<0),thSliceInv(real(Idq)<0));
                        FqSlice(real(Idq)<0)=+fInt.Fq(-IdSlice(real(Idq)<0),IqSlice(real(Idq)<0),thSliceInv(real(Idq)<0));
                        TSlice(real(Idq)<0)=-fInt.T(-IdSlice(real(Idq)<0),IqSlice(real(Idq)<0),thSliceInv(real(Idq)<0));
                    end
                    if (abs(IqMin)<IqMax)
                        FdSlice(imag(Idq)<0)=NaN;
                        FqSlice(imag(Idq)<0)=NaN;
                        TSlice(imag(Idq)<0)=NaN;
                    end
                else
                    % symmetry just on q axis
                    if (IdMax<abs(IdMin))
                        FdSlice(real(Idq)>0)=NaN;
                        FqSlice(real(Idq)>0)=NaN;
                        TSlice(real(Idq)>0)=NaN;
                    end
                    if (abs(IqMin)<IqMax)
                        FdSlice(imag(Idq)<0)=+fInt.Fd(IdSlice(imag(Idq)<0),-IqSlice(imag(Idq)<0),thSliceInv(imag(Idq)<0));
                        FqSlice(imag(Idq)<0)=-fInt.Fq(IdSlice(imag(Idq)<0),-IqSlice(imag(Idq)<0),thSliceInv(imag(Idq)<0));
                        TSlice(imag(Idq)<0)=-fInt.T(IdSlice(imag(Idq)<0),-IqSlice(imag(Idq)<0),thSliceInv(imag(Idq)<0));
                    end 
                end
            end
            
            
            FdTmp=FdTmp+FdSlice;
            FqTmp=FqTmp+FqSlice;
            TTmp=TTmp+TSlice;
            if debug
                thPlot=reshape(thTmp,[tS,1,1]);
                TPlot=reshape(TSlice,[tS,1,1]);
                plot(thPlot,TPlot,'-o')
            end
        end
        
        FdTmp=FdTmp/length(alfa_k);
        FqTmp=FqTmp/length(alfa_k);
        TTmp=TTmp/length(alfa_k);
        if debug
            thPlot=reshape(thTmp,[tS,1,1]);
            TPlot=reshape(TTmp,[tS,1,1]);
            plot(thPlot,TPlot,'-ro','LineWidth',2);
            keyboard
        end
        data.Fd(dd,qq,:)=FdTmp;
        data.Fq(dd,qq,:)=FqTmp;
        data.T(dd,qq,:)=TTmp;
        
%         % figura
%         xData=get(hplot,'XData');
%         yData=get(hplot,'yData');
%         xData=[xData data.Id(dd,qq,1)];
%         yData=[yData data.Iq(dd,qq,1)];
%         set(hplot,'XData',xData,'yData',yData);
%         pause(0.01)

        indexTime = (qq+(dd-1)*qS)/(dS*qS);
%         disp(['dqtMap skewing - ' int2str(qq/qS*100) '% of ' int2str(dd/dS*100) '%'])
%         disp(['dqtMap skewing - ' int2str(indexTime*100) '%'])
        if ~isempty(hEdit)
            set(hEdit,'Value',[int2str(indexTime*100) '%']);
            drawnow();
        else
            disp(['  dqtMap skewing  - ' int2str(indexTime*100) '%'])
        end
    end
end

data.Fd = data.Fd+Lld*data.Id;
data.Fq = data.Fq+Llq*data.Iq;

data0 = dqtMap.data;
dqtMap.data = data;

fInt.Id = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Id,'linear','none');
fInt.Iq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Iq,'linear','none');
fInt.th = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.th,'linear','none');
fInt.Fd = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fd,'linear','none');
fInt.Fq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fq,'linear','none');
fInt.T  = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.T,'linear','none');

dqtMap.fInt=fInt;

% Save output data
Id   = mean(dqtMap.data.Id,3);
Iq   = mean(dqtMap.data.Iq,3);
Fd   = griddedInterpolant(Id,Iq,mean(dqtMap.data.Fd,3),'linear','none');
Fq   = griddedInterpolant(Id,Iq,mean(dqtMap.data.Fq,3),'linear','none');
T    = griddedInterpolant(Id,Iq,mean(dqtMap.data.T,3),'linear','none');
dT   = griddedInterpolant(Id,Iq,std(dqtMap.data.T,0,3,'omitnan'),'linear','none');
dTpp = griddedInterpolant(Id,Iq,max(dqtMap.data.T,[],3)-min(dqtMap.data.T,[],3),'linear','none');

[Id,Iq]=ndgrid(linspace(min(id),max(id),256),linspace(min(iq),max(iq),256));

Fd   = Fd(Id,Iq);
Fq   = Fq(Id,Iq);
T    = T(Id,Iq);
dT   = dT(Id,Iq);
dTpp = dTpp(Id,Iq);

fdfq.Id   = Id';
fdfq.Iq   = Iq';
fdfq.Fd   = Fd';
fdfq.Fq   = Fq';
fdfq.T    = T';
fdfq.dT   = dT';
fdfq.dTpp = dTpp';





