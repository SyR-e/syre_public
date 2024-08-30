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
            geo.PMdim = -[ones(1,geo.nlay);zeros(1,geo.nlay)];
        elseif (strcmp(geo.RotType,'Seg')||strcmp(geo.RotType,'ISeg'))
            geo.PMdim = -ones(2,geo.nlay);
        else
            geo.PMdim = zeros(1,geo.nlay);
        end
    else
        if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Vtype'))
            geo.PMdim = zeros(2,geo.nlay);
        elseif (strcmp(geo.RotType,'Seg')||strcmp(geo.RotType,'ISeg'))
            geo.PMdim = zeros(2,geo.nlay);
        else
            geo.PMdim = zeros(1,geo.nlay);
        end
    end
    
    
    
    flagPM = 0;
    flagHC = 0;
    flagPT = 0;
    flagPR = 0;
    for ii=1:length(geo.RQnames)
        eval(['geo.' geo.RQnames{ii} ' = ' num2str(RQ(ii)) ';'])
        if contains(geo.RQnames{ii},'PMdim')
            flagPM = 1;
        end
        if contains(geo.RQnames{ii},'FBS')
            geo.th_FBS = geo.th_FBS*pi/180;
        end
        if contains(geo.RQnames{ii},'hc_pu')
            flagHC = 1;
        end
        if contains(geo.RQnames{ii},'pontT')
            flagPT = 1;
        end
        if contains(geo.RQnames{ii},'pontR')
            flagPR = 1;
        end
    end
    
    if flagPM
        geo.PMdim = -geo.PMdim; % negative during optimization
    end
    
    if flagHC
        geo.hc_pu = sort(geo.hc_pu);
    end

    if flagPT
        geo.pontT = sort(geo.pontT);
    end
    
    if flagPR
        geo.pontR = sort(geo.pontR);
    end
    
    if strcmp(geo.RQnames{ii},'gamma')
        geo = rmfield(geo,'gamma');
    end
    
end

if strcmp(geo.RotType,'SPM-Halbach')        
    geo.dalpha_pu = 1;              
else
    if sum(geo.dalpha_pu) > (1-0.05) % 0.05 is the angular space guaranteed for the spider
        if geo.dalpha_pu(1)>(1-0.05*(geo.nlay)) % max geo.dalpha_pu(1)=1-0.05-0.05*(nlay-1)
            geo.dalpha_pu(1)=1-0.05*(geo.nlay);
        end
        geo.dalpha_pu(2:end) = geo.dalpha_pu(2:end)/sum(geo.dalpha_pu(2:end))*(1-0.05-geo.dalpha_pu(1));
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

if strcmp(geo.RotType,'SPM')
    geo.hc = geo.hc_pu*geo.g;
end

% if strcmp(geo.RotType,'SPM')&&ismember('hc(1)',geo.RQnames)
%     geo.lm = geo.hc_pu;
% end