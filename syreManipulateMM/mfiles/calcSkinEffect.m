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

function [kAC] = calcSkinEffect(skinEffectModel,freq,temp,method)

switch skinEffectModel.type
    case '0'
        kAC = 1;
    case 'interpFreq'
        if strcmp(method,'LUT')
            kAC = interp1(skinEffectModel.f,skinEffectModel.k,freq);
        else
            kAC = polyval(skinEffectModel.p,freq);
        end
    case 'interpFreqTemp'
        if strcmp(method,'LUT')
            kAC = interp2(skinEffectModel.f,skinEffectModel.T,skinEffectModel.k,freq,temp);
        else
            f = skinEffectModel.f(skinEffectModel.T==temp);
            k = skinEffectModel.k(skinEffectModel.T==temp);
            [p,~] = polyfit(f,k,7);
            kAC = polyval(p,freq);
        end
end