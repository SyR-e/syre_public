% Copyright 2022
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


function [Id,Iq] = calcSkewCurrentDomain(lim,thSKW,nPoints,axisType)


% ang_sk = ang_sk_m*p*pi/180; % elt rad

IdMin = lim.IdMin;
IdMax = lim.IdMax;
IqMin = lim.IqMin;
IqMax = lim.IqMax;

if strcmp(axisType,'SR')
    if IdMin<0
        numQuad = 4;
    else
        if IqMin<0
            numQuad = 2;
        else
            numQuad = 1;
        end
    end
    Iabs = (IdMax^2+IqMax^2)^0.5;
    Iangle = atan2(IqMax,IdMax);
else
    if IqMin<0
        numQuad = 4;
    else
        if IdMax>0
            numQuad = 2;
        else
            numQuad = 1;
        end
    end
    Iabs = (IdMin^2+IqMax^2)^0.5;
    Iangle = atan2(IqMax,IdMin);
end

switch axisType
    case 'SR'
        lim.IdMax = Iabs*cos(Iangle+thSKW);
        lim.IqMax = Iabs*sin(Iangle-thSKW);
        switch numQuad
            case 1
                lim.IdMin = 0;
                lim.IqMin = 0;
                lim.nD = nPoints;
                lim.nQ = nPoints;
            case 2
                lim.IdMin = 0;
                lim.IqMin = -lim.IqMax;
                lim.nD = nPoints;
                lim.nQ = 2*nPoints-1;
            case 3
                lim.IdMin = -lim.IdMax;
                lim.IqMin = -lim.IqMax;
                lim.nD = 2*nPoints-1;
                lim.nQ = 2*nPoints-1;
        end
    case 'PM'
        lim.IdMin = Iabs*cos(Iangle-thSKW);
        lim.IqMax = Iabs*sin(Iangle+thSKW);
        switch numQuad
            case 1
                lim.IdMax = 0;
                lim.IqMin = 0;
                lim.nD = nPoints;
                lim.nQ = nPoints;
            case 2
                lim.IdMax = -lim.IdMin;
                lim.IqMin = 0;
                lim.nD = 2*nPoints-1;
                lim.nQ = nPoints;
            case 4
                lim.IdMax = -lim.IdMin;
                lim.IqMin = -lim.IqMax;
                lim.nD = 2*nPoints-1;
                lim.nq = 2*nPoints-1;
        end
end

Id = linspace(lim.IdMin,lim.IdMax,lim.nD);
Iq = linspace(lim.IqMin,lim.IqMax,lim.nQ);
% [Id,Iq] = meshgrid(Id,Iq);
