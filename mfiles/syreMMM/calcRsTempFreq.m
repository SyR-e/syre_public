% Copyright 2022
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

function [Rs,kAC] = calcRsTempFreq(Rs0,temp0,l,lend,acLossFactor,skinEffectMethod,temp,freq)

% calcolo della resistenza a 20°C
R20 = Rs0/(1+0.004*(temp0-20));

if strcmp(skinEffectMethod,'0')
    kAC = 1;
else
    kAC = calcSkinEffect(acLossFactor,abs(freq),abs(temp),skinEffectMethod);
end

Rs  = R20.*(kAC*l/(lend+l)+lend/(lend+l)).*(1+0.004*(temp-20));
