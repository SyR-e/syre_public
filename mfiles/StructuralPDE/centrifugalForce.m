% Copyright 2021
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

function [f] = centrifugalForce(region,y,dataForCF)

kgm3_Fe  = dataForCF.kgm3_Fe;
kgm3_PM  = dataForCF.kgm3_PM;
w        = dataForCF.w;
psMagnet = dataForCF.psMagnet;


x0 = 0;
y0 = 0;

xy = (region.x-x0)+j*(region.y-y0);

kgm3 = kgm3_Fe*ones(size(xy));
if ~isempty(psMagnet)
    for ii=1:length(xy)
        [in,out] = intersect(psMagnet,[real(xy(ii))*[1;1],imag(xy(ii))*[1;1]]);
        if isempty(out)
            kgm3(ii)=kgm3_PM;
        end
    end
end


% plot(real(xy),imag(xy),'r.')

r  = abs(xy);
th = angle(xy);
% f = kgm3*V.*w.^2.*[r.*cos(th);r.*sin(th)];
f = kgm3.*w.^2.*[r.*cos(th);r.*sin(th)];
