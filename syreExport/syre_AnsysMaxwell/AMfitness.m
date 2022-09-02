% Copyright 2021
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



function [geo,mat,out,filepath] = AMfitness(geo,per,mat,eval_type,filepath,filename,ipypath)

[SOL]=simulate_xdeg_AM(geo,per,mat,eval_type,filepath,filename,ipypath);

th0 = geo.th0;
%from cls files we obtain column vectors
T     = SOL.T; %[t(ms),T(Nm)]
F     = SOL.F; %[t(ms),f1(Wb),f2(Wb),f3(Wb)]
Theta = geo.p*SOL.Theta; %[t(ms),Theta(deg_e)]
I     = SOL.I; %[t(ms),i1(A),i2(A),i3(A)]
%V     = SOL.V; %[t(ms),v1(V),v2(V),v3(V)]

[id,iq] = abc2dq_AM(I(:,1)',I(:,2)',I(:,3)',Theta',th0);
[fd,fq] = abc2dq_AM(F(:,1)',F(:,2)',F(:,3)',Theta',th0);
IPF     = sin(atan(iq./id)-atan(fq./fd));

%output function's value
out.th   = Theta;
out.T    = T;
out.fd   = fd;
out.fq   = fq;
out.id   = id;
out.iq   = iq;
out.IPF  = IPF;
out.xdeg = SOL.xdeg;

if per.corelossflag==1
    %CoreLoss=SOL.CoreLoss{:,2}; %[t(ms),CoreLoss(W)]
    out.CoreLoss   = SOL.CoreLoss;
    out.CoreLoss_s = SOL.CoreLoss_s;
    out.CoreLoss_r = SOL.CoreLoss_r;
    out.PMloss     = SOL.PMloss;
    
    %  out.CoreLossAvg = SOL.CoreLossAvg{1,15};
    %  out.HysteresisLossAvg = SOL.HysteresisLossAvg{1,15};
    %  out.EddyLossAvg = SOL.EddyLossAvg{1,16};
    %  out.ExcessLossAvg = SOL.ExcessLossAvg{1,15};
    
    %   out.CoreLossAvg = SOL.CoreLossAvg{1,6};
    %  out.HysteresisLossAvg = SOL.HysteresisLossAvg{1,6};
    %  out.EddyLossAvg = SOL.EddyLossAvg{1,6};
    %  out.ExcessLossAvg = SOL.ExcessLossAvg{1,6};
end

end