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

function [idiq] = MMM_eval_inverse_dq_Simulink(motorModel)

    %% Direct Flux Maps
    Id = motorModel.FluxMap_dq.Id;
    Iq = motorModel.FluxMap_dq.Iq;
    Fd = motorModel.FluxMap_dq.Fd;
    Fq = motorModel.FluxMap_dq.Fq;
    
    %% Motor type
    
    if strcmp(motorModel.data.axisType,'SR') && strcmp(motorModel.data.motorType,'SR')
        Quad_Maps = 0; %SyR Convention - 1st quadrant maps
    elseif strcmp(motorModel.data.axisType,'SR') && strcmp(motorModel.data.motorType,'PM')
        Quad_Maps = 1; %PM-SyR - 1st and 4th quadrant maps
    elseif strcmp(motorModel.data.axisType,'PM') && strcmp(motorModel.data.motorType,'PM')
        Quad_Maps = 2; %IPM - 1st and 2st quadrant maps
    end


    %% Number of points
    % should be an even number 
    Np_ref = 256;       

%% Compute Inverse Flux Maps

    %Compute the grid limits in the flux domain

    FdMin = min(Fd,[],'all');
    FdMax = max(Fd,[],'all');
    FqMin = min(Fq,[],'all');
    FqMax = max(Fq,[],'all');

    % Compute the regular flux domain 
    fD_vect = linspace(FdMin,FdMax,Np_ref+1);
    fQ_vect = linspace(FqMin,FqMax,Np_ref+1);
    [fD,fQ]=meshgrid(fD_vect,fQ_vect);
    
    %Enumeration of the Elements
    index = 1:1:numel(Id);

    % Create Vectors for scattered interpolant
    % (scattered interpolant wants column vectors)
    
    IdVect = Id(index)';
    IqVect = Iq(index)';
    FdVect = Fd(index)';
    FqVect = Fq(index)';
    
    % Filt NaN
    IdVect = IdVect(~isnan(FdVect));
    IqVect = IqVect(~isnan(FdVect));
    FdVect = FdVect(~isnan(FdVect));
    FqVect = FqVect(~isnan(FdVect));

    % Compute Interpolating Function of current function of fluxes
    intD = scatteredInterpolant(FdVect,FqVect,IdVect,'natural','none');
    intQ = scatteredInterpolant(FdVect,FqVect,IqVect,'natural','none');
    
    % Inverse Flux Maps 
    iD = intD(fD,fQ);
    iQ = intQ(fD,fQ);
    
    % Manual Corrections
    switch(Quad_Maps)
        case{0,2}
            iD(1,:) = iD(2,:);
            iQ(1,:) = iQ(2,:);
        case 1
            iD(:,1) = iD(:,2);
            iQ(:,1) = iQ(:,2);
    end

%% Flux to p.u. Error Maps
    
    Nd_ref = Np_ref+1;
    Nq_ref = Np_ref+1;

%% Remove NaN limits of 

    switch(Quad_Maps)
    case {0,2}
        % Remove NaN limits along d-axis
        range_x = ones(1,length(iD(1,:)));
        for i=1:length(iD(1,:))
            tmp1=iD(:,i);
            % find valid elements
            tmp2 = tmp1(~isnan(tmp1));  % index is given by logical array
            % if elements is NaN
            if(isempty(tmp2))   
                range_x(1,i)=0; % along row all elements are NaN
            end
        end
        
        % Remove NaN limits along q-axis
        range_y = ones(length(iD(:,1)),1);
        for i=1:length(iD(:,1))
            tmp1 = iD(i,:);
            tmp2=tmp1(~isnan(tmp1));
            if(isempty(tmp2))
                range_y(i,1)=0; % along column all elements are NaN
            end
        end
    case 1
     
        % Remove NaN limits along d-axis
        range_x = ones(1,length(iQ(1,:)));
        for i=1:length(iQ(1,:))
            tmp1=iQ(:,i);
            % find valid elements
            tmp2 = tmp1(~isnan(tmp1));  % index is given by logical array
            % if elements is NaN
            if(isempty(tmp2))   
                range_x(1,i)=0; % along row all elements are NaN
            end
        end
        
        % Remove NaN limits along q-axis
        range_y = ones(length(iQ(:,1)),1);
        for i=1:length(iQ(:,1))
            tmp1 = iQ(i,:);
            tmp2=tmp1(~isnan(tmp1));
            if(isempty(tmp2))
                range_y(i,1)=0; % along column all elements are NaN
            end
        end
    end
    % Extract Valid grid point
    % Maps are rescaled eliminating all NaN rows/columns
    tmp_x = find(range_x>0.5);
    tmp_y = find(range_y>0.5);
    fD_rescaled  =  fD(tmp_y,tmp_x);
    fQ_rescaled  =  fQ(tmp_y,tmp_x);
    iD_rescaled  =  iD(tmp_y,tmp_x);
    iQ_rescaled  =  iQ(tmp_y,tmp_x);

%% Compute Inversion based on Motor Type 
    switch(Quad_Maps)
        case {0,2}
            %% Define Maximum q-axis Flux function of d-axis Flux
            
            fDVct = fD_rescaled(1,:);
            fQMax = zeros(1,length(iD_rescaled(1,:)));
            ixMax = zeros(1,length(iD_rescaled(1,:)));
            for i=1:length(iD_rescaled(1,:))
                tmp1 = iD_rescaled(:,i);
                % find index of first occourence of NaN, then minus one for the last
                % valid element
                % tmp2 = find(~isnan(tmp1)<0.5,1)-1; 
                tmp2 = find(~isnan(tmp1)==1);
                if(isempty(tmp2)) % if all elements are valid
                    fQMax(1,i) = fQ_rescaled(end,1);
                    ixMax(1,i) = length(iD_rescaled(:,1));
                else
                    fQMax(1,i) = fQ_rescaled(tmp2(end),1);
                    ixMax(1,i) = tmp2(end);
                end   
            end
        
        
            %% p.u interpolation along q-axis
      
            %Define fq vct in p.u.
            fQVct_pu = linspace(0,1,Nq_ref)';
            
            %Define normalized matrix along q-axis
            iD_pu_q = zeros(Nq_ref,length(iD_rescaled(1,:)));
            iQ_pu_q = zeros(Nq_ref,length(iD_rescaled(1,:)));
            
            % p.u. Interpolation along q-axis
            for i=1:length(iD_rescaled(1,:))
                tmp_fq_unorm = fQ_rescaled(1:ixMax(i),i);
                tmp_id_unorm = iD_rescaled(1:ixMax(i),i);
                tmp_iq_unorm = iQ_rescaled(1:ixMax(i),i);
                tmp_fq_norm  = fQVct_pu*fQMax(1,i);
                tmp_id_unorm = interp1(tmp_fq_unorm,tmp_id_unorm,tmp_fq_norm,'spline');
                tmp_iq_unorm = interp1(tmp_fq_unorm,tmp_iq_unorm,tmp_fq_norm,'spline');
                iD_pu_q(:,i) = tmp_id_unorm;
                iQ_pu_q(:,i) = tmp_iq_unorm;
            end
        
            %% p.u. interpolation along d-axis
            
            fDVct_pu = linspace(0,1,Nd_ref);
            tmp_fd_min = min(fD_rescaled(1,:))*ones(size(fD_rescaled(1,:)));
            tmp_fd_max = max(fD_rescaled(1,:))*ones(size(fD_rescaled(1,:)));
            tmp_fd_unorm = (fD_rescaled(1,:)-tmp_fd_min)./(tmp_fd_max-tmp_fd_min);
            
            %Define normlized matrix
            Fd_vct = fDVct_pu*(max(fD_rescaled(1,:))-min(fD_rescaled(1,:)))+min(fD_rescaled(1,:))*ones(size(fDVct_pu));
            Fq_max_vct = interp1(fDVct,fQMax,Fd_vct,'spline');
            iD_pu = interp2(tmp_fd_unorm,fQVct_pu,iD_pu_q,fDVct_pu,fQVct_pu);
            iQ_pu = interp2(tmp_fd_unorm,fQVct_pu,iQ_pu_q,fDVct_pu,fQVct_pu);
            [fD_pu,fQ_pu] = meshgrid(fDVct_pu,fQVct_pu); 
            %% Rename Variables for saving
                    
            idiq.fD_vct_ref = Fd_vct;
            idiq.fQ_vct_max = Fq_max_vct;
            idiq.fD_pu_norm = fD_pu;
            idiq.fQ_pu_norm = fQ_pu;
            idiq.iD_pu_norm = iD_pu;
            idiq.iQ_pu_norm = iQ_pu;
    case 1
            %% Define Maximum d-axis Flux function of q-axis Flux
           
            fQVct = fQ_rescaled(:,1);
            fDMax = zeros(length(iQ_rescaled(:,1)),1);
            fDMin = zeros(length(iQ_rescaled(:,1)),1);
            
            for i=1:length(iQ_rescaled(:,1))
                tmp1 = iQ_rescaled(i,:);
                % find index of first occourence of NaN, then minus one for the last
                % valid element
                   
                tmp2 = find(~isnan(tmp1)==1);
                if(isempty(tmp2)) % if all elements are valid
                    fDMax(i,1) = fD_rescaled(1,end);
                    ixMax(i,1) = length(iQ_rescaled(:,1));               
                else
                    fDMax(i,1) = fD_rescaled(1,tmp2(end));
                    ixMax(i,1) = tmp2(end);
                 
                end   
            end
        
        
            %% p.u interpolation along d-axis
            
            %Define fq vct in p.u.
            fDVct_pu = linspace(0,1,Nd_ref);
            
            %Define normalized matrix along d-axis
            iD_pu_q = zeros(length(iQ_rescaled(:,1)),Nd_ref);
            iQ_pu_q = zeros(length(iQ_rescaled(:,1)),Nd_ref);
            
            % p.u. Interpolation along d-axis
            for i=1:length(iQ_rescaled(:,1))
                tmp_fd_unorm = fD_rescaled(i,1:ixMax(i));
                tmp_id_unorm = iD_rescaled(i,1:ixMax(i));
                tmp_iq_unorm = iQ_rescaled(i,1:ixMax(i));
                tmp_fd_norm  = fDVct_pu*fDMax(i,1);
                tmp_id_unorm = interp1(tmp_fd_unorm,tmp_id_unorm,tmp_fd_norm,'spline');
                tmp_iq_unorm = interp1(tmp_fd_unorm,tmp_iq_unorm,tmp_fd_norm,'spline');
                iD_pu_q(i,:) = tmp_id_unorm;
                iQ_pu_q(i,:) = tmp_iq_unorm;
            end
        
        
            %% p.u. interpolation along q-axis
            
            fQVct_pu = linspace(0,1,Nq_ref)';
            tmp_fq_min = min(fQ_rescaled(:,1))*ones(size(fQ_rescaled(:,1)));
            tmp_fq_max = max(fQ_rescaled(:,1))*ones(size(fQ_rescaled(:,1)));
            tmp_fq_unorm = (fQ_rescaled(:,1)-tmp_fq_min)./(tmp_fq_max-tmp_fq_min);
            
            %Define normlized matrix
            Fq_vct = fQVct_pu*(max(fQ_rescaled(:,1))-min(fQ_rescaled(:,1)))+min(fQ_rescaled(:,1))*ones(size(fQVct_pu));
            Fd_max_vct = interp1(fQVct,fDMax,Fq_vct,'spline');
            iD_pu = interp2(fDVct_pu,tmp_fq_unorm,iD_pu_q,fDVct_pu,fQVct_pu);
            iQ_pu = interp2(fDVct_pu,tmp_fq_unorm,iQ_pu_q,fDVct_pu,fQVct_pu);
            [fD_pu,fQ_pu] = meshgrid(fDVct_pu,fQVct_pu); 
       
            %% Rename Variables for saving
            
            idiq.fQ_vct_ref = Fq_vct;
            idiq.fD_vct_max = Fd_max_vct;
            idiq.fD_pu_norm = fD_pu;
            idiq.fQ_pu_norm = fQ_pu;
            idiq.iD_pu_norm = iD_pu;
            idiq.iQ_pu_norm = iQ_pu;

    end

end