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

function [geo,gamma,mat] = interpretRQ(RQ,geo,mat)

% INTERPRETRQ interprets the string of the MOOA inputs (RQ)
% and returns dalpha, hc_pu, dx, gamma
% dalpha per unit

if not(isempty(geo.RQnames))
    if ~strcmp(mat.LayerMag.MatName,'Air')
        if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Vtype'))
            geo.PMdim = -ones(1,geo.nlay);
        elseif (strcmp(geo.RotType,'Seg')||strcmp(geo.RotType,'ISeg'))
            geo.PMdim = -ones(2,geo.nlay);
        else
            geo.PMdim = zeros(1,geo.nlay);
        end
    else
        if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Vtype'))
            geo.PMdim = zeros(1,geo.nlay);
        elseif (strcmp(geo.RotType,'Seg')||strcmp(geo.RotType,'ISeg'))
            geo.PMdim = zeros(2,geo.nlay);
        else
            geo.PMdim = zeros(1,geo.nlay);
        end
    end
    
    first_index = 1;
    if strcmp(geo.RQnames{first_index},'dalpha')
        last_index = first_index + geo.nlay-1;
        geo.dalpha_pu = RQ(first_index:last_index);
    else
        last_index = 0;
    end
    
    % the sum of pu angles is rescaled to one
    % alpha1 is not rescaled: only the angles from the second barrier onwards
    if sum(geo.dalpha_pu) > (1-0.05) % 0.05 is the angular space guaranteed for the spider
        if geo.dalpha_pu(1)>(1-0.05*(geo.nlay)) % max geo.dalpha_pu(1)=1-0.05-0.05*(nlay-1)
            geo.dalpha_pu(1)=1-0.05*(geo.nlay);
        end
        geo.dalpha_pu(2:end) = geo.dalpha_pu(2:end)/sum(geo.dalpha_pu(2:end))*(1-0.05-geo.dalpha_pu(1));
    end
    
    first_index = last_index + 1;
    if (length(geo.RQnames)>(last_index+1)||last_index==0)
        if strcmp(geo.RQnames{first_index},'hc')
            % hc per unit
            last_index = first_index + geo.nlay - 1;
            geo.hc_pu = RQ(first_index:last_index);
            % debug
            geo.hc_pu = sort(geo.hc_pu);
        end
    end
    
    first_index = last_index + 1;
    if (length(geo.RQnames)>(last_index+1)||last_index==0)
        if strcmp(geo.RQnames{first_index},'dx')
            % dx per unit
            last_index = first_index + geo.nlay - 1;
            geo.dx = RQ(first_index:last_index);
        end
    end
    
    first_index = last_index + 1;
    if length(geo.RQnames)>(last_index+1)
        if strcmp(geo.RQnames{first_index},'Br')
            % Br [T]
            last_index = first_index + geo.nlay - 1;
            mat.LayerMag.Br = RQ(first_index:last_index);
        end
    end
    
    first_index = last_index + 1;
    if length(geo.RQnames)>(last_index+1)
        if strcmp(geo.RQnames{first_index},'betaPMshape')
            % Br [T]
            last_index = first_index + geo.nlay - 1;
            geo.betaPMshape = RQ(first_index:last_index);
        end
    end
    
    flagPM=0;
    if length(geo.RQnames)>(last_index)
        for k = last_index+1:length(geo.RQnames)
            eval(['geo.' geo.RQnames{k} ' = ' num2str(RQ(k)) ';'])
            if strcmp(geo.RQnames{k},'VanglePM')||strcmp(geo.RQnames{k},'th_FBS')
                eval(['geo.' geo.RQnames{k} ' = ' num2str(RQ(k)) ' *pi/180;']);
            end
            if strcmp(geo.RQnames{k},'PMdim(1)')
                flagPM=1;
            end
        end
    end
    
    if flagPM
        geo.PMdim = -geo.PMdim; % PMdim negative to use the per-unit notation (just during optimization)
    end

    
end

geo.dalpha = geo.dalpha_pu*(90/geo.p);  % [mec degrees]
% current phase angle
if (not(isempty(RQ))&&strcmp(geo.RQnames{end},'gamma'))
    gamma = RQ(end);
else
    if strcmp(geo.RotType,'SPM')
        gamma = 90;
    elseif strcmp(geo.RotType,'Vtype')
        gamma = 135;
    else
        gamma = 45;
    end
end

if strcmp(geo.RotType,'SPM')&&ismember('hc',geo.RQnames)
    geo.lm = geo.hc_pu;
end