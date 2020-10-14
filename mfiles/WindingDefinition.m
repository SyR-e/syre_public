% Copyright 2020
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

function [win] = WindingDefinition(dataSet)
% This function define the winding pattern inside the GUI callbacks.
% It is based on Koil, as before, but supports also multi-3phase windings,
% by extending the 3phase circuit.

% Load data
p    = dataSet.NumOfPolePairs;
q    = dataSet.NumOfSlots;
n3ph = dataSet.Num3PhaseCircuit;
yq   = dataSet.PitchShortFac;

if n3ph>1
    yq=1;
end

Q    = 6*n3ph*p*q;
Q1   = 6*p*q;

% Call Koil_syre for the baseline 3phase winding
path = pwd;
cd(fullfile(path,'koil'));
system(['koil_syre.exe',' ',num2str(Q1),' ',num2str(p),' ',num2str(yq*3*q)]);
cd(path);
win1 = MatrixWin();

% align the win matrix with the phase 1 all at the beginning
if q>1
    tmp  = sum(win1(1,1:q)~=ones(1,q));
    win1 = [win1(:,end-tmp+1:end) win1(:,1:end-tmp)];
end


% Define the correct periodicity of the motor and keep just the
% corresponing slots
t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
Q1s = Q1/t2;
% win1 = win1(:,1:floor(Q1s));
% ps = floor(Q1s/(3*q));


[nL,nS] = size(win1);

win = zeros(nL,n3ph*nS);
Qs = Q1s*n3ph;

if q>1
    indexBase = 1:1:q;
    indexSet  = zeros(1,nS);
    for ii=1:nS/q
        ini=1+q*(ii-1);
        fin=ini+q-1;
        indexSet(ini:fin)=indexBase+(ii-1)*n3ph*q;
    end
    for ii=1:n3ph
        win(:,indexSet+q*(ii-1))=(abs(win1)+3*(ii-1)).*sign(win1);
    end
    
else
    for ii=1:n3ph
        win(:,ii:n3ph:end)=(abs(win1)+3*(ii-1)).*sign(win1);
    end
end


