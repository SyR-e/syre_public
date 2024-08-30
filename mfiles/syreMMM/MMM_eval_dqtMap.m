% Copyright 2023
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

    if exist('dataSet','var')
        n3phase = dataSet.Num3PhaseCircuit;
        th0 = geo.th0;
    else
        n3phase = 1;
    end

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

    if rem(nRep,1)==0

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
                    if isfield(SOL{indQ,indD},'we')
                        WeTmp = SOL{indQ,indD}.we;
                        WcTmp = SOL{indQ,indD}.wc;
                    else
                        WeTmp = zeros(size(th));
                        WcTmp = zeros(size(th));
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
                    WeTmp = zeros(size(th));
                    WcTmp = zeros(size(th));
                end

                Id = repmat(IdTmp,1,nRep);
                Iq = repmat(IqTmp,1,nRep);
                Fd = repmat(FdTmp,1,nRep);
                Fq = repmat(FqTmp,1,nRep);
                T  = repmat(TTmp,1,nRep);
                We = repmat(WeTmp,1,nRep);
                Wc = repmat(WcTmp,1,nRep);

                % single set elaboration
                if n3phase>1
                    for ii=1:n3phase
                        Iph = phaseQuantityDecoding(IaTmp,IbTmp,IcTmp,xdeg);
                        Ia = Iph.a(ii,:);
                        Ib = Iph.b(ii,:);
                        Ic = Iph.c(ii,:);

                        Fph = phaseQuantityDecoding(FaTmp,FbTmp,FcTmp,xdeg);
                        Fa = Fph.a(ii,:);
                        Fb = Fph.b(ii,:);
                        Fc = Fph.c(ii,:);

                        setsTmp(ii).ia = Ia;
                        setsTmp(ii).ib = Ib;
                        setsTmp(ii).ic = Ic;
                        setsTmp(ii).fa = Fa;
                        setsTmp(ii).fb = Fb;
                        setsTmp(ii).fc = Fc;

                        % dq transformation (single set)
                        idq = abc2dq(Ia,Ib,Ic,(thVect+th0(ii)-th0(1))*pi/180);
                        setsTmp(ii).id = idq(1,:);
                        setsTmp(ii).iq = idq(2,:);
                        fdq = abc2dq(Fa,Fb,Fc,(thVect+th0(ii)-th0(1))*pi/180);
                        setsTmp(ii).fd = fdq(1,:);
                        setsTmp(ii).fq = fdq(2,:);
                    end
                end

                Iph = phaseQuantityDecoding(IaTmp,IbTmp,IcTmp,xdeg);
                Ia = Iph.a(1,:);
                Ib = Iph.b(1,:);
                Ic = Iph.c(1,:);

                Fph = phaseQuantityDecoding(FaTmp,FbTmp,FcTmp,xdeg);
                Fa = Fph.a(1,:);
                Fb = Fph.b(1,:);
                Fc = Fph.c(1,:);

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
                We = We(index);
                Wc = Wc(index);

                if n3phase>1
                    for ii=1:n3phase
                        setsTmp(ii).id = setsTmp(ii).id(index);
                        setsTmp(ii).iq = setsTmp(ii).iq(index);
                        setsTmp(ii).fd = setsTmp(ii).fd(index);
                        setsTmp(ii).fq = setsTmp(ii).fq(index);
                        setsTmp(ii).ia = setsTmp(ii).ia(index);
                        setsTmp(ii).ib = setsTmp(ii).ib(index);
                        setsTmp(ii).ic = setsTmp(ii).ic(index);
                        setsTmp(ii).fa = setsTmp(ii).fa(index);
                        setsTmp(ii).fb = setsTmp(ii).fb(index);
                        setsTmp(ii).fc = setsTmp(ii).fc(index);
                    end
                end

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
                    dqtMap.data.We(indD,indQ,indT) = We(indT);
                    dqtMap.data.Wc(indD,indQ,indT) = Wc(indT);

                    dqtMap.Id(indD) = Id(indT);
                    dqtMap.Iq(indQ) = Iq(indT);

                    if n3phase>1
                        for ii=1:n3phase
                            dqtMap.sets(ii).Id(indD,indQ,indT) = setsTmp(ii).id(indT);
                            dqtMap.sets(ii).Iq(indD,indQ,indT) = setsTmp(ii).iq(indT);
                            dqtMap.sets(ii).Fd(indD,indQ,indT) = setsTmp(ii).fd(indT);
                            dqtMap.sets(ii).Fq(indD,indQ,indT) = setsTmp(ii).fq(indT);
                            dqtMap.sets(ii).Ia(indD,indQ,indT) = setsTmp(ii).ia(indT);
                            dqtMap.sets(ii).Ib(indD,indQ,indT) = setsTmp(ii).ib(indT);
                            dqtMap.sets(ii).Ic(indD,indQ,indT) = setsTmp(ii).ic(indT);
                            dqtMap.sets(ii).Fa(indD,indQ,indT) = setsTmp(ii).fa(indT);
                            dqtMap.sets(ii).Fb(indD,indQ,indT) = setsTmp(ii).fb(indT);
                            dqtMap.sets(ii).Fc(indD,indQ,indT) = setsTmp(ii).fc(indT);
                        end
                    end
                end
            end
        end

        % ndmesh for compatibility
        [dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th]=ndgrid(dqtMap.Id,dqtMap.Iq,dqtMap.th);
    else
        dqtMap = [];
    end


    % Interpolant
    %     fInt.Id = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Id,'spline');
    %     fInt.Iq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Iq,'spline');
    %     fInt.th = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.th,'spline');
    %     fInt.Fd = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fd,'spline');
    %     fInt.Fq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fq,'spline');
    %     fInt.T  = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.T,'spline');
    %     fInt.Fa = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fa,'spline');
    %     fInt.Fb = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fb,'spline');
    %     fInt.Fc = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fc,'spline');
    %     fInt.Ia = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Ia,'spline');
    %     fInt.Ib = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Ib,'spline');
    %     fInt.Ic = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Ic,'spline');
    %     fInt.We = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.We,'spline');
    %     fInt.Wc = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Wc,'spline');

    % dqtMap.Fmap=Fmap;
    %dqtMap.fInt = fInt;
    % dqtMap.fdfq=fdfq;

    % Gli assi delle correnti sono scambiati tra 2D e 3D:
    %  2D --> (Iq,Id)
    %  3D --> (Id,Iq,th)

else
    dqtMap = [];
end