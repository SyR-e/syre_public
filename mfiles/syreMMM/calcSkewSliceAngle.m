% Copyright 2026
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

function [alpha_k,rot_alpha_k] = calcSkewSliceAngle(motorModel)

ang_sk_m = motorModel.tmpSkew.thSkw;
nSlice   = motorModel.tmpSkew.nSlice;
sk_shape = motorModel.tmpSkew.shape;
p        = motorModel.data.p;

% NB: Herringbone skewing have 2x on ang_skew because it is performed on
% half of the slices of the linear skew

switch sk_shape
    case 'Linear skewing'
        ang_sk = ang_sk_m*p*pi/180; % elt rad
        k = 1:1:nSlice;
        k = k-mean(k);
    case 'Herringbone (V-shape) skewing'
        ang_sk = 2*ang_sk_m*p*pi/180; % elt rad
        if rem(nSlice,2)==0
            k = 1:1:nSlice/2;
            k = [k fliplr(k)];
        else
            k = 1:1:(nSlice+1)/2;
            k = [k fliplr(k(1:end-1))];
        end
        k = k-mean(k);
end

alpha_k = k*ang_sk/(nSlice);
rot_alpha_k = exp(-1i*alpha_k);