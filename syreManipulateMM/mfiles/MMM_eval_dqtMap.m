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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dqtMap] = MMM_eval_dqtMap(pathname,filename)

% if nargin()==0
%     pathname = [cd '\'];
% else
%     pathname = motorModel.data.pathname;
% end

if nargin()==0
    pathname = [cd '\'];
    [filename,pathname] = uigetfile([pathname '*.mat'],'Load F_map file for dqtMap evaluation');
elseif nargin()==1
    [filename,pathname] = uigetfile([pathname '*.mat'],'Load F_map file for dqtMap evaluation');
    
end

if filename
    load([pathname filename]);
    if ~(exist('OUT','var')||exist('SOL','var'))
        error('SOL or OUT structure not present. Impossible to build the 3D map')
    end
    
    %% dqtMap
    % notes on back compatibility.
    % There are three cases possible for dqtMap working.
    % 1) last version is based on rev417, with the correct name of the structure (SOL). It is a cell matrix and each cell is the out.SOL strucrure from FEMMfitness
    % 2) past version: the name of the variable was wrong (OUT instead of SOL), but the content was exactly the same
    % 3) first version: the cell matrix OUT was based on the out.SOL matrix (and not the structure), so each columns of out.SOL had a meaning. To identify this case, the flag_struct variable is used.
    
    
    if exist('OUT','var')
        if isstruct(OUT{1,1})
            flag_struct=1;
        else
            flag_struct=0;
        end
        SOL = OUT;
    else
        flag_struct = 1;
    end
    
    if flag_struct
        thTmp=SOL{1,1}.th;
    else
        thTmp=SOL{1,1}(:,1)';
    end
    nsim=length(thTmp);
    xdeg=floor(abs(thTmp(end)-thTmp(1))/(nsim-1)*nsim);
    
    nRep=(360/xdeg); % number of repetition
    
    thVect=[];
    
    for ii=1:nRep
        thPar=thTmp+xdeg*(ii-1);
        thVect=[thVect thPar];
    end
    
    thVect(thVect<0)=360+thVect(thVect<0);
    [th,index]=sort(thVect);
    
    dqtMap.th=th;
    
    for indQ=1:size(SOL,1)
        for indD=1:size(SOL,2)
            if flag_struct
                IdTmp = SOL{indQ,indD}.id;
                IqTmp = SOL{indQ,indD}.iq;
                FdTmp = SOL{indQ,indD}.fd;
                FqTmp = SOL{indQ,indD}.fq;
                TTmp  = (SOL{indQ,indD}.T);
                if isfield(SOL{indQ,indD},'fa')
                    FaTmp = SOL{indQ,indD}.fa;
                    FbTmp = SOL{indQ,indD}.fb;
                    FcTmp = SOL{indQ,indD}.fc;
                else
                    FaTmp = zeros(size(th));
                    FbTmp = zeros(size(th));
                    FcTmp = zeros(size(th));
                end
                if isfield(SOL{indQ,indD},'ia')
                    IaTmp = SOL{indQ,indD}.ia;
                    IbTmp = SOL{indQ,indD}.ib;
                    IcTmp = SOL{indQ,indD}.ic;
                else
                    IaTmp = zeros(size(th));
                    IbTmp = zeros(size(th));
                    IcTmp = zeros(size(th));
                end
            else
                IdTmp = SOL{indQ,indD}(:,2)';
                IqTmp = SOL{indQ,indD}(:,3)';
                FdTmp = SOL{indQ,indD}(:,4)';
                FqTmp = SOL{indQ,indD}(:,5)';
                TTmp  = (SOL{indQ,indD}(:,6))';
                FaTmp = zeros(size(th));
                FbTmp = zeros(size(th));
                FcTmp = zeros(size(th));
                IaTmp = zeros(size(th));
                IbTmp = zeros(size(th));
                IcTmp = zeros(size(th));
            end
            
            Id = repmat(IdTmp,1,nRep);
            Iq = repmat(IqTmp,1,nRep);
            Fd = repmat(FdTmp,1,nRep);
            Fq = repmat(FqTmp,1,nRep);
            T  = repmat(TTmp,1,nRep);
            
            switch xdeg
                case 360
                    Fa = FaTmp;
                    Fb = FbTmp;
                    Fc = FcTmp;
                    
                    Ia = IaTmp;
                    Ib = IbTmp;
                    Ic = IcTmp;
                case 60
                    Fa = [FaTmp -FbTmp FcTmp -FaTmp FbTmp -FcTmp];
                    Fb = [FbTmp -FcTmp FaTmp -FbTmp FcTmp -FaTmp];
                    Fc = [FcTmp -FaTmp FbTmp -FcTmp FaTmp -FbTmp];
                    
                    Ia = [IaTmp -IbTmp IcTmp -IaTmp IbTmp -IcTmp];
                    Ib = [IbTmp -IcTmp IaTmp -IbTmp IcTmp -IaTmp];
                    Ic = [IcTmp -IaTmp IbTmp -IcTmp IaTmp -IbTmp];
                case 120
                    Fa = [FaTmp FcTmp FbTmp];
                    Fb = [FbTmp FaTmp FcTmp];
                    Fc = [FcTmp FbTmp FaTmp];
                    
                    Ia = [IaTmp IcTmp IbTmp];
                    Ib = [IbTmp IaTmp IcTmp];
                    Ic = [IcTmp IbTmp IaTmp];
                case 180
                    Fa = [FaTmp -FaTmp];
                    Fb = [FbTmp -FbTmp];
                    Fc = [FcTmp -FcTmp];
                    
                    Ia = [IaTmp -IaTmp];
                    Ib = [IbTmp -IbTmp];
                    Ic = [IcTmp -IcTmp];
                    
            end
            
            Id = Id(index);
            Iq = Iq(index);
            Fd = Fd(index);
            Fq = Fq(index);
            T  = T(index);
            Fa = Fa(index);
            Fb = Fb(index);
            Fc = Fc(index);
            Ia = Ia(index);
            Ib = Ib(index);
            Ic = Ic(index);
            
            for indT=1:length(th)
                dqtMap.data.Fd(indD,indQ,indT) = Fd(indT);
                dqtMap.data.Fq(indD,indQ,indT) = Fq(indT);
                dqtMap.data.T(indD,indQ,indT)  = T(indT);
                dqtMap.data.Fa(indD,indQ,indT) = Fa(indT);
                dqtMap.data.Fb(indD,indQ,indT) = Fb(indT);
                dqtMap.data.Fc(indD,indQ,indT) = Fc(indT);
                dqtMap.data.Ia(indD,indQ,indT) = Ia(indT);
                dqtMap.data.Ib(indD,indQ,indT) = Ib(indT);
                dqtMap.data.Ic(indD,indQ,indT) = Ic(indT);
                
                dqtMap.Id(indD) = Id(indT);
                dqtMap.Iq(indQ) = Iq(indT);
            end
        end
    end
    
    % ndmesh because matlab is stupid...
    [dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th]=ndgrid(dqtMap.Id,dqtMap.Iq,dqtMap.th);
    
    % F_map
    Fmap.Id   = mean(dqtMap.data.Id,3);
    Fmap.Iq   = mean(dqtMap.data.Iq,3);
    Fmap.Fd   = mean(dqtMap.data.Fd,3);
    Fmap.Fq   = mean(dqtMap.data.Fq,3);
    Fmap.T    = mean(dqtMap.data.T,3);
    Fmap.dT   = std(dqtMap.data.T,0,3);
    Fmap.dTpp = max(dqtMap.data.T,[],3)-min(dqtMap.data.T,[],3);
    
    % fdfq_idiq_n256
    fdfq.Id=linspace(min(min(Fmap.Id)),max(max(Fmap.Id)),256);
    fdfq.Iq=linspace(min(min(Fmap.Iq)),max(max(Fmap.Iq)),256);
    [fdfq.Id,fdfq.Iq]=meshgrid(fdfq.Id,fdfq.Iq);
    fdfq.Fd   = interp2(Fmap.Id(:,1),Fmap.Iq(1,:),Fmap.Fd',fdfq.Id,fdfq.Iq,'spline');
    fdfq.Fq   = interp2(Fmap.Id(:,1),Fmap.Iq(1,:),Fmap.Fq',fdfq.Id,fdfq.Iq,'spline');
    fdfq.T    = interp2(Fmap.Id(:,1),Fmap.Iq(1,:),Fmap.T',fdfq.Id,fdfq.Iq,'spline');
    fdfq.dT   = interp2(Fmap.Id(:,1),Fmap.Iq(1,:),Fmap.dT',fdfq.Id,fdfq.Iq,'spline');
    fdfq.dTpp = interp2(Fmap.Id(:,1),Fmap.Iq(1,:),Fmap.dTpp',fdfq.Id,fdfq.Iq,'spline');
    
    % Interpolant
    fInt.Id = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Id,'spline');
    fInt.Iq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Iq,'spline');
    fInt.th = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.th,'spline');
    fInt.Fd = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fd,'spline');
    fInt.Fq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fq,'spline');
    fInt.T  = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.T,'spline');
    fInt.Fa = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fa,'spline');
    fInt.Fb = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fb,'spline');
    fInt.Fc = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fc,'spline');
    
    % dqtMap.Fmap=Fmap;
    dqtMap.fInt=fInt;
    % dqtMap.fdfq=fdfq;
    
    % Gli assi delle correnti sono scambiati tra 2D e 3D:
    %  2D --> (Iq,Id)
    %  3D --> (Id,Iq,th)
    
else
    dqtMap = [];
end