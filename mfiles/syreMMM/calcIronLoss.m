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


function [Pfe,Pfes_h,Pfes_c,Pfer_h,Pfer_c,Ppm]=calcIronLoss(IronLossModel,fdfq,FreqElet)

Id=fdfq.Id;
Iq=fdfq.Iq;
Fd=fdfq.Fd;
Fq=fdfq.Fq;

if strcmp(IronLossModel.type,'0')
    Pfe=zeros(size(Id));
elseif strcmp(IronLossModel.type,'map')
    IdMap=IronLossModel.Id;
    IqMap=IronLossModel.Iq;
    
    if isfield(IronLossModel,'Pfeh')
        Pfeh = interp2(IdMap,IqMap,IronLossModel.Pfe_h,Id,Iq,'cubic',1e50);
        Pfec = interp2(IdMap,IqMap,IronLossModel.Pfe_c,Id,Iq,'cubic',1e50);
        Ppm  = interp2(IdMap,IqMap,IronLossModel.Ppm,Id,Iq,'cubic',1e50);
        
        Pfeh = Pfeh*(FreqElet./IronLossModel.f0).^IronLossModel.expH;
        Pfec = Pfec*(FreqElet./IronLossModel.f0).^IronLossModel.expC;
        Ppm  = IronLossModel.segPM*Ppm*(FreqElet./IronLossModel.f0).^IronLossModel.expPM;
    else
        Pfes_h = interp2(IdMap,IqMap,IronLossModel.Pfes_h,Id,Iq,'cubic',NaN);
        Pfes_c = interp2(IdMap,IqMap,IronLossModel.Pfes_c,Id,Iq,'cubic',NaN);
        Pfer_h = interp2(IdMap,IqMap,IronLossModel.Pfer_h,Id,Iq,'cubic',NaN);
        Pfer_c = interp2(IdMap,IqMap,IronLossModel.Pfer_c,Id,Iq,'cubic',NaN);
        Ppm    = interp2(IdMap,IqMap,IronLossModel.Ppm,Id,Iq,'cubic',NaN);
        
        Pfes_h = Pfes_h*(FreqElet./IronLossModel.f0).^IronLossModel.expH;
        Pfes_c = Pfes_c*(FreqElet./IronLossModel.f0).^IronLossModel.expC;
        Pfer_h = Pfer_h*(FreqElet./IronLossModel.f0).^IronLossModel.expH;
        Pfer_c = Pfer_c*(FreqElet./IronLossModel.f0).^IronLossModel.expC;
        Ppm  = IronLossModel.segPM*Ppm*(FreqElet./IronLossModel.f0).^IronLossModel.expPM;
        Pfeh = Pfes_h+Pfer_h;
        Pfec = Pfes_c+Pfer_c;
    end
    
    
    Pfe=Pfeh+Pfec+Ppm;
    
elseif strcmp(IronLossModel.type,'fitIM1')
    a  = IronLossModel.param.a;
    b  = IronLossModel.param.b;
    k1 = IronLossModel.param.k1;
    k2 = IronLossModel.param.k2;
    k3 = IronLossModel.param.k3;
    f0 = IronLossModel.param.f0;
    c1 = IronLossModel.param.c1;
    c2 = IronLossModel.param.c2;
    c3 = IronLossModel.param.c3;

    Pfe=k1*(fs/f0).^c1.*Fs.^a+k2*(fs/f0).^c2.*Iq.^b-k3*(fs/f0).^c3.*Fs.*Iq;
    
    Pfes_h = NaN;
    Pfes_c = NaN;
    Pfer_h = NaN;
    Pfer_c = NaN;
    Ppm    = NaN;
    
elseif strcmp(IronLossModel.type,'point')
    Pfe0=IronLossModel.Pfe0;
    F0=IronLossModel.F0;
    f0=IronLossModel.f0;
    expFlux=IronLossModel.expFlux;
    expFreq=IronLossModel.expFreq;
    
    Fabs=abs(Fd+j*Fq);
    
    Pfe=Pfe0*(Fabs./F0).^expFlux*(FreqElet./f0).^expFreq;
    Pfes_h = NaN;
    Pfes_c = NaN;
    Pfer_h = NaN;
    Pfer_c = NaN;
    Ppm    = NaN;
    
elseif strcmp(IronLossModel.type,'fitExpSyRMTPA')
    k1=IronLossModel.param.k1;
    k2=IronLossModel.param.k2;
    k3=IronLossModel.param.k3;
    e1=IronLossModel.param.e1;
    e2=IronLossModel.param.e2;
    
    Iabs=abs(Id+j*Iq);
    Pfe=k1*FreqElet.^e1+k2*Iabs.^e2+k3*FreqElet.*Iabs;
    Pfe(Pfe<0)=0;
    
    Pfes_h = NaN;
    Pfes_c = NaN;
    Pfer_h = NaN;
    Pfer_c = NaN;
    Ppm    = NaN;
end




