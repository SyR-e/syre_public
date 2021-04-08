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

function [AOA] = MMM_MTPAload(pathname)

[filename,pathname] = uigetfile([pathname 'AOA\' '*.mat'],'Select AOA folder');

if filename
    tmp = load([pathname 'ktMax_idiq.mat']);
    MTPA.id   = tmp.id_KtMax;
    MTPA.iq   = tmp.iq_KtMax;
    MTPA.T    = tmp.T_KtMax;
    MTPA.dTpp = tmp.dTpp_KtMax;
    MTPA.I    = abs(MTPA.id+j*MTPA.iq);

    tmp = load([pathname 'ktMax_FdFq.mat']);
    MTPA.fd   = tmp.fd_KtMax;
    MTPA.fq   = tmp.fq_KtMax;
    
    tmp = load([pathname 'kvMax_idiq.mat']);
    MTPV.id   = tmp.id_KvMax;
    MTPV.iq   = tmp.iq_KvMax;
    MTPV.T    = tmp.T_KvMax;
    
    tmp = load([pathname 'kvMax_FdFq.mat']);
    MTPV.fd   = tmp.fd_KvMax;
    MTPV.fq   = tmp.fq_KvMax;
    
    AOA.MTPA = MTPA;
    AOA.MTPV = MTPV;
    AOA.method = 'LUT';
    
else
    AOA = [];
end