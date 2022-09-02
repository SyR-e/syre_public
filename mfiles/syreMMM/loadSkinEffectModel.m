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

function SkinEffModel=loadSkinEffectModel(filename)

if strcmp(filename,'0')
    SkinEffModel.type='0';
else
    load(filename)
    if exist('results','var') % from FEA simulations
        if isfield(results,'T')
            SkinEffModel.type='interpFreqTemp';
            f = results.f;
            T = results.T;
            k = results.k;
            
            if f(1,1)~=0
                f = [zeros(size(f,1),1) f];
                k = [ones(size(k,1),1) k];
                T = [T(:,1) T];
            end
            
            SkinEffModel.f = f;
            SkinEffModel.T = T;
            SkinEffModel.k = k;
        else
            f=results.f(1,:);
            k=results.k(1,:);
            if f(1)~=0
                f = [0 f];
                k = [1 k];
            end

            [p,s] = polyfit(f,k,7);

            SkinEffModel.type='interpFreq';
            SkinEffModel.f=f;
            SkinEffModel.k=k;
            SkinEffModel.p=p;
            SkinEffModel.s=s;
            SkinEffModel.n=7;
        end
    elseif exist('skinEffect','var') % from MMM GUI
        SkinEffModel = skinEffect;
    else
        SkinEffModel = [];
    end
end
        