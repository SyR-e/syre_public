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

function [demagCurve] = MMM_load_demag(filename)

if nargin~=1
    [filename,pathname,~] = uigetfile([cd '\*.mat'],'Select Demagnetization Curve');
    filename = [pathname filename];
end

if isempty(filename)
    demagCurve = [];
else
    load(filename,'demagLimit')
    if isfield(demagLimit,'Idemag')
        demagCurve = demagLimit;
    elseif isfield(demagLimit,'current')
        demagCurve.Idemag = demagLimit.current;
        demagCurve.tempPM = demagLimit.temperature;
        demagCurve.dPM    = demagLimit.dPM;
        demagCurve.Bmin   = demagLimit.Bmin;
        
        test.Idemag = [];
        test.tempPM = [];
        test.dPM    = [];
        test.Bmin   = [];
    
        for ii=1:length(demagLimit.temperature)
            fieldName = ['temp_' int2str(ii)];
            iTmp = demagLimit.iterations.(fieldName).Iiter;
            dTmp = demagLimit.iterations.(fieldName).dPMiter;
            BTmp = demagLimit.iterations.(fieldName).Biter;
    
            test.Idemag = [test.Idemag iTmp];
            test.tempPM = [test.tempPM demagLimit.temperature(ii)*ones(size(iTmp))];
            test.dPM    = [test.dPM dTmp];
            test.Bmin   = [test.Bmin BTmp];
        end
        demagCurve.test = test;
    else
        demagCurve = [];
    end
end







