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
% For the VFM dTpp is set to NaN (change it if necessary) 

function [AOA] = MMM_eval_AOA(motorModel,method,factors)


if nargin()<2
    method = 'LUT';
end

if nargin()<3
    factors.kRefine   = 5;
    factors.kRefine3D = 3;
    factors.levels    = 201;
end

kRefine   = factors.kRefine;
kRefine3D = factors.kRefine3D;
numLevels = factors.levels;


%% Load data
fdfq = motorModel.FluxMap_dq;
axisType = motorModel.data.axisType;



if strcmp(motorModel.data.motorType,'EE')
    ThirdDim = 'Ir';
elseif strcmp(motorModel.data.motorType,'VFM')
    ThirdDim = 'MS';
end

if ~(strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType,'VFM'))
        fdfq = mapsReInterpolation(fdfq,'Id','Iq',kRefine*size(fdfq.Id,1),'linear');
else
    nmax = 255;
    disp('Control Trajectories Evaluation...')
    if kRefine3D*size(fdfq.Id,1) <= nmax && kRefine3D*size(fdfq.Id,3)<= nmax
        fdfq = mapsReInterpolation3D(fdfq,'Id','Iq', ThirdDim,kRefine3D*size(fdfq.Id,1),'linear');
    elseif size(fdfq.Id,3)~=size(fdfq.Id,1)
        kRefine3D = nmax/size(fdfq.Id,1);                    % set the maximum number of iterpolation points to nmax to avoid 'out of memory' error
        fdfq = mapsReInterpolation3D(fdfq,'Id','Iq', ThirdDim,kRefine3D*size(fdfq.Id,1),'linear');
    end
end

Id   = fdfq.Id;
Iq   = fdfq.Iq;
Fd   = fdfq.Fd;
Fq   = fdfq.Fq;

if strcmp(motorModel.data.motorType,'VFM') % Here the MS demag and remag limit are taken into accounts
    if any(fdfq.MS>10,"all")               % If the range of MS is between 0 and 100
        fdfq.MS=fdfq.MS/100;               % Correction to MS between 0 and 1
    end
    Remag = motorModel.VFMdata.data.Remag;
    Demag = motorModel.VFMdata.data.Demag;
    index = 1:numel(Remag.id);
    idRem_vect = reshape(Remag.id(index),[],1);
    iqRem_vect = reshape(Remag.iq(index),[],1);
    MsRem_vect = reshape(Remag.MS(index),[],1);
    index = 1:numel(Demag.id);
    idDem_vect = reshape(Demag.id(index),[],1);
    iqDem_vect = reshape(Demag.iq(index),[],1);
    MsDem_vect = reshape(Demag.MS(index),[],1);

    F_idDemLim = scatteredInterpolant(iqDem_vect, MsDem_vect, idDem_vect, 'linear', 'none');    %Id fucntion limit for demag
    F_idRemLim = scatteredInterpolant(iqRem_vect, MsRem_vect, idRem_vect, 'linear', 'none');    %Id function  Limit for remag
    PF = nan(size(fdfq.Id,2),size(fdfq.Iq,1),size(fdfq.MS,3));
    dTpp = nan(size(fdfq.Id,2),size(fdfq.Iq,1),size(fdfq.MS,3));
    T = nan(size(fdfq.Id,2),size(fdfq.Iq,1),size(fdfq.MS,3));
    Fd_ =  nan(size(fdfq.Id,2),size(fdfq.Iq,1),size(fdfq.MS,3));     % Flux d to be used in MTPF, considering demag remag limit
    Fq_ =  nan(size(fdfq.Id,2),size(fdfq.Iq,1),size(fdfq.MS,3));     % Flux q to be used in MTPF, considering demag remag limit
    for kk = 1:size(fdfq.MS,3)      % Torque computation taking into account the id Lim 
        for jj = 1:size(fdfq.Iq,1)     
            for ii = 1:size(fdfq.Id,2)
                if (fdfq.Id(jj,ii,kk) < 0 && fdfq.Id(jj,ii,kk) < F_idDemLim(fdfq.Iq(jj,ii,kk),fdfq.MS(jj,ii,kk))) || (fdfq.Id(jj,ii,kk) >0 && fdfq.Id(jj,ii,kk) > F_idRemLim(fdfq.Iq(jj,ii,kk),fdfq.MS(jj,ii,kk))) || (fdfq.Id(jj,ii,kk) < 0 && isnan(F_idDemLim(fdfq.Iq(jj,ii,kk),fdfq.MS(jj,ii,kk)))) || (fdfq.Id(jj,ii,kk) >0 && isnan(F_idRemLim(fdfq.Iq(jj,ii,kk),fdfq.MS(jj,ii,kk)))) 
                    %Empty on purpose: the value are alreaddy initialized as NaN 
                else
                    T(jj,ii,kk) = 1.5*motorModel.data.p*(fdfq.Fd(jj,ii,kk)*fdfq.Iq(jj,ii,kk)-fdfq.Fq(jj,ii,kk)*fdfq.Id(jj,ii,kk));
                    Fd_(jj,ii,kk) = fdfq.Fd(jj,ii,kk);
                    Fq_(jj,ii,kk) = fdfq.Fq(jj,ii,kk);
                    PF(jj,ii,kk) = sin(atan2(fdfq.Iq(jj,ii,kk),fdfq.Id(jj,ii,kk))-atan2(fdfq.Fq(jj,ii,kk),fdfq.Fd(jj,ii,kk)));
                end
            end
        end
    end
else
    T    = fdfq.T;
    dTpp = fdfq.dTpp;
    PF   = sin(atan2(Iq,Id)-atan2(Fq,Fd));
end


%% Extract curves

%MTPA
if ~(strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType,'VFM'))
    Ilevels = linspace(0,max(abs(Id(:)+j*Iq(:))),numLevels);
    [id,iq] = calcOptCtrl(Id,Iq,T,abs(Id+j*Iq),axisType,Ilevels);
else
    dim_3D = size(fdfq.(ThirdDim),3);
    max_opt = 501;
    id = NaN(dim_3D,max_opt);
    iq = NaN(dim_3D,max_opt);
    Ivect = linspace(0,max(abs(Id(:)+j*Iq(:))),max_opt);
    for ii = 1:dim_3D
        % if strcmp(motorModel.data.motorType,'VFM')
            [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii),Iq(:,:,ii),T(:,:,ii),abs(Id(:,:,ii)+j*Iq(:,:,ii)),axisType,Ivect,1);
        % else
        %     [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii)',Iq(:,:,ii)',T(:,:,ii)',abs(Id(:,:,ii)'+j*Iq(:,:,ii)'),axisType,Ivect,1);
        % end
    end
end

% Check MTPA
if strcmp(motorModel.data.motorType,'VFM')
    id_max = max(id,[],"all");
    iq_max = max(iq,[],"all");
    for kk = 1:dim_3D
        index = find(hypot(id(kk,:),iq(kk,:)) > min(id_max,iq_max));
        id(kk,index)= NaN;
        iq(kk,index)= NaN;
    end
elseif strcmp(motorModel.data.axisType,'SR')
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    % index = find(id==min(Id,[],'all'));
    % id(index) = NaN;
    % iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    % index = find(iq<0);
    % id(index) = NaN;
    % iq(index) = NaN;
else
    % index = find(id==max(Id,[],'all'));
    % id(index) = NaN;
    % iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    % index = find(id>0);
    % id(index) = NaN;
    % iq(index) = NaN;
end

if ~(strcmp(motorModel.data.motorType,'EE')||strcmp(motorModel.data.motorType,'VFM'))
    index = ~isnan(id);
    id = id(index);
    iq = iq(index);
    
    index = ~isnan(iq);
    id = id(index);
    iq = iq(index);
    
    if ((id(1)~=0)&&(iq(1)~=0))
        id = [0 id];
        iq = [0 iq];
    end
end

MTPA.id   = id;
MTPA.iq   = iq;

if strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType, 'VFM')
    temp_ThirdDim   = fdfq.(ThirdDim)(1,1,:);
    MTPA.([ThirdDim '_layers']) = temp_ThirdDim(:);
    MTPV.([ThirdDim '_layers']) = temp_ThirdDim(:);
end


% MTPA.fd   = interp2(Id,Iq,Fd,id,iq);
% MTPA.fq   = interp2(Id,Iq,Fq,id,iq);
% MTPA.T    = interp2(Id,Iq,T,id,iq);
% MTPA.dTpp = interp2(Id,Iq,dTpp,id,iq);

% MTPV
if ~(strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType,'VFM'))
    Flevels = linspace(0,max(abs(Fd(:)+j*Fq(:))),numLevels);
    [id,iq] = calcOptCtrl(Id,Iq,T,abs(Fd+j*Fq),axisType,Flevels);
else
    dim_3D = size(fdfq.(ThirdDim),3);
    max_opt = 501;
    id = NaN(dim_3D,max_opt);
    iq = NaN(dim_3D,max_opt);
    Fvect = linspace(0,max(abs(Fd(:)+j*Fq(:))),max_opt);
    for ii = 1:dim_3D
        % if strcmp(motorModel.data.motorType,'VFM')
            [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii),Iq(:,:,ii),T(:,:,ii),abs(Fd(:,:,ii)+j*Fq(:,:,ii)),axisType,Fvect,1);
        % else
        %     [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii)',Iq(:,:,ii)',T(:,:,ii)',abs(Fd(:,:,ii)'+j*Fq(:,:,ii)'),axisType,Fvect,1);
        % end
    end
end


% Check MTPV
if strcmp(motorModel.data.motorType, 'VFM')
    for kk = 1:dim_3D
        mask = false(size(id(1,:)));
        index = find(hypot(id(kk,:),iq(kk,:)) > min(id_max,iq_max));
        id(kk,index)= NaN;
        iq(kk,index)= NaN;
        if kk==30
            a=1;
        end
        [~,consecutive_idx] = findNonNaNGroups(id(kk,:));
        mask(consecutive_idx) = true;  
        id(kk,~mask) = NaN;
    end
elseif strcmp(motorModel.data.axisType,'SR')
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq<0);
    id(index) = NaN;
    iq(index) = NaN;
else
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id>0);
    id(index) = NaN;
    iq(index) = NaN;
end

if ~(strcmp(motorModel.data.motorType,'EE')||strcmp(motorModel.data.motorType,'VFM'))
    index = ~isnan(id);
    id = id(index);
    iq = iq(index);
    
    index = ~isnan(iq);
    id = id(index);
    iq = iq(index);
end

MTPV.id = id;
MTPV.iq = iq;

% MPFPA
if ~(strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType,'VFM'))
    [id,iq] = calcOptCtrl(Id,Iq,PF,abs(Id+j*Iq),axisType,Ilevels);
else
    dim_3D = size(fdfq.(ThirdDim),3);
    max_opt = 501;
    id = NaN(dim_3D,max_opt);
    iq = NaN(dim_3D,max_opt);
    Ivect = linspace(0,max(abs(Id(:)+j*Iq(:))),max_opt);
    for ii = 1:dim_3D
        % if strcmp(motorModel.data.motorType,'VFM')
            [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii),Iq(:,:,ii),PF(:,:,ii),abs(Id(:,:,ii)+j*Iq(:,:,ii)),axisType,Ivect,1);
        % else
        %     [id(ii,:),iq(ii,:)] = calcOptCtrl(Id(:,:,ii)',Iq(:,:,ii)',PF(:,:,ii)',abs(Id(:,:,ii)'+j*Iq(:,:,ii)'),axisType,Ivect,1);
        % end
    end
end


if  strcmp(motorModel.data.motorType, 'VFM')
    for kk = 1:dim_3D
        index = find(hypot(id(kk,:),iq(kk,:)) > min(id_max,iq_max));
        id(kk,index)= NaN;
        iq(kk,index)= NaN;
    end
elseif strcmp(motorModel.data.axisType,'SR')
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
else
    index = find(id==max(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(id==min(Id,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==max(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
    index = find(iq==min(Iq,[],'all'));
    id(index) = NaN;
    iq(index) = NaN;
end

MPFPA.id = id;
MPFPA.iq = iq;

%% interpolate (if needed)

nPoints = 101;


if strcmp(method,'Fit')
    if strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType,'VFM')
        warning('Fit option not implemented yet for EESM and VFM, proceeding with LUT instead!!!')
    else
        ft = fittype('poly8');
        if strcmp(axisType,'SR')
            opts = fitoptions(...
                'Method','LinearLeastSquares',...
                'Lower',[zeros(1,8) 0],...
                'Upper',[inf*ones(1,8),0]);
            [xData,yData] = prepareCurveData(MTPA.id,MTPA.iq);
            fitFun = fit(xData,yData,ft,opts);
            MTPA.id = linspace(0,max(MTPA.id),nPoints)';
            MTPA.iq = fitFun(MTPA.id);
            if ~isempty(MTPV.id)
                opts = fitoptions(...
                    'Method','LinearLeastSquares',...
                    'Lower',[zeros(1,8)],...
                    'Upper',[inf*ones(1,8)]);
                [xData,yData] = prepareCurveData(MTPV.id,MTPV.iq);
                fitFun = fit(xData,yData,ft,opts);
                MTPV.id = linspace(0,max(MTPV.id),nPoints)';
                MTPV.iq = fitFun(MTPV.id);
            end
        else
            opts = fitoptions(...
                'Method','LinearLeastSquares',...
                'Upper',[zeros(1,8) 0],...
                'Lower',[-inf*ones(1,8),0]);
            [xData,yData] = prepareCurveData(MTPA.iq,MTPA.id);
            fitFun = fit(xData,yData,ft,opts);
            MTPA.iq = linspace(0,max(MTPA.iq),nPoints)';
            MTPA.id = fitFun(MTPA.iq);
            if ~isempty(MTPV.iq)
                opts = fitoptions(...
                    'Method','LinearLeastSquares',...
                    'Upper',[zeros(1,8)],...
                    'Lower',[-inf*ones(1,8)]);
                [xData,yData] = prepareCurveData(MTPV.iq,MTPV.id);
                fitFun = fit(xData,yData,ft,opts);
                MTPV.iq = linspace(0,max(MTPV.iq),nPoints)';
                MTPV.id = fitFun(MTPV.iq);
            end
        end
        AOA.method = 'Fit';
        
        % filt for trajectories outside the domain
        index = find(MTPA.id>max(Id,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
        index = find(MTPA.id<min(Id,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
        index = find(MTPA.iq>max(Iq,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
        index = find(MTPA.iq<min(Iq,[],'all'));
        MTPA.id(index) = NaN;
        MTPA.iq(index) = NaN;
    
        index = ~isnan(MTPA.id);
        MTPA.id = MTPA.id(index);
        MTPA.iq = MTPA.iq(index);
        
        index = find(MTPV.id>max(Id,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
        index = find(MTPV.id<min(Id,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
        index = find(MTPV.iq>max(Iq,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
        index = find(MTPV.iq<min(Iq,[],'all'));
        MTPV.id(index) = NaN;
        MTPV.iq(index) = NaN;
    
        index = ~isnan(MTPV.id);
        MTPV.id = MTPV.id(index);
        MTPV.iq = MTPV.iq(index);
    end
else
    AOA.method = 'LUT';
end


%% Interp flux linkages and torque


if ~(strcmp(motorModel.data.motorType,'EE') || strcmp(motorModel.data.motorType,'VFM'))
    MTPA.fd    = interp2(Id,Iq,Fd,MTPA.id,MTPA.iq);
    MTPA.fq    = interp2(Id,Iq,Fq,MTPA.id,MTPA.iq);
    MTPA.T     = interp2(Id,Iq,T,MTPA.id,MTPA.iq);
    MTPA.dTpp  = interp2(Id,Iq,dTpp,MTPA.id,MTPA.iq);
    
    MTPV.fd    = interp2(Id,Iq,Fd,MTPV.id,MTPV.iq);
    MTPV.fq    = interp2(Id,Iq,Fq,MTPV.id,MTPV.iq);
    MTPV.T     = interp2(Id,Iq,T,MTPV.id,MTPV.iq);
    MTPV.dTpp  = interp2(Id,Iq,dTpp,MTPV.id,MTPV.iq);
    
    MPFPA.fd   = interp2(Id,Iq,Fd,MPFPA.id,MPFPA.iq);
    MPFPA.fq   = interp2(Id,Iq,Fq,MPFPA.id,MPFPA.iq);
    MPFPA.T    = interp2(Id,Iq,T,MPFPA.id,MPFPA.iq);
    MPFPA.dTpp = interp2(Id,Iq,dTpp,MPFPA.id,MPFPA.iq);
else
    dim_3D = size(fdfq.(ThirdDim),3);
    for ii=1:dim_3D
        % if strcmp(motorModel.data.motorType,'VFM')
            Id_vect = Id(1,:,1);
            Iq_vect = Iq(:,1,1);
        % else    
        %     Id_vect = Id(:,1,1);
        %     Iq_vect = Iq(1,:,1);
        % end
        [Id_mg, Iq_mg] = meshgrid(Id_vect(:), Iq_vect(:));

        MTPA.fd(ii,:)    = interp2(Id_mg,Iq_mg,Fd(:,:,ii),MTPA.id(ii,:),MTPA.iq(ii,:));
        MTPA.fq(ii,:)    = interp2(Id_mg,Iq_mg,Fq(:,:,ii),MTPA.id(ii,:),MTPA.iq(ii,:));
        MTPA.T(ii,:)     = interp2(Id_mg,Iq_mg,T(:,:,ii),MTPA.id(ii,:),MTPA.iq(ii,:));
        MTPA.dTpp(ii,:)  = interp2(Id_mg,Iq_mg,dTpp(:,:,ii),MTPA.id(ii,:),MTPA.iq(ii,:));
        
        MTPV.fd(ii,:)    = interp2(Id_mg,Iq_mg,Fd(:,:,ii),MTPV.id(ii,:),MTPV.iq(ii,:));
        MTPV.fq(ii,:)    = interp2(Id_mg,Iq_mg,Fq(:,:,ii),MTPV.id(ii,:),MTPV.iq(ii,:));
        MTPV.T(ii,:)     = interp2(Id_mg,Iq_mg,T(:,:,ii),MTPV.id(ii,:),MTPV.iq(ii,:));
        MTPV.dTpp(ii,:)  = interp2(Id_mg,Iq_mg,dTpp(:,:,ii),MTPV.id(ii,:),MTPV.iq(ii,:));
    
        MPFPA.fd(ii,:)   = interp2(Id_mg,Iq_mg,Fd(:,:,ii),MPFPA.id(ii,:),MPFPA.iq(ii,:));
        MPFPA.fq(ii,:)   = interp2(Id_mg,Iq_mg,Fq(:,:,ii),MPFPA.id(ii,:),MPFPA.iq(ii,:));
        MPFPA.T(ii,:)    = interp2(Id_mg,Iq_mg,T(:,:,ii),MPFPA.id(ii,:),MPFPA.iq(ii,:));
        MPFPA.dTpp(ii,:) = interp2(Id_mg,Iq_mg,dTpp(:,:,ii),MPFPA.id(ii,:),MPFPA.iq(ii,:));
    end
end

AOA.MTPA  = MTPA;
AOA.MTPV  = MTPV;
AOA.MPFPA = MPFPA;










