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

function [Inductance] = MMM_eval_inductanceMap(motorModel)

% Computation of the incremental inductances

Id = motorModel.fdfq.Id;
Iq = motorModel.fdfq.Iq;
Fd = motorModel.fdfq.Fd;
Fq = motorModel.fdfq.Fq;

[dFdd,dFdq] = gradient(Fd);
[dFqd,dFqq] = gradient(Fq);
[dIdd,~]    = gradient(Id);
[~,dIqq]    = gradient(Iq);

Ldd = dFdd./dIdd;
Ldq = dFdq./dIqq;
Lqd = dFqd./dIdd;
Lqq = dFqq./dIqq;

% [dFdd,dFdq] = gradient(Fd);
% [dFqd,dFqq] = gradient(Fq);
% [dIdd,~]    = gradient(Id);
% [~,dIqq]    = gradient(Iq);

% Ldd = dFdd./dIdd;
% Ldq = dFdq./dIqq;
% Lqq = dFqq./dIqq;
% Lqd = dFqd./dIdd;

% elaboration from C_ldlq_idiq.m
% ldd = diff(Fd,1,2)./diff(Id,1,2);
% ldq = diff(Fd,1,1)./diff(Iq,1,1);
%
% lqq = diff(Fq,1,1)./diff(Iq,1,1);
% lqd = diff(Fq,1,2)./diff(Id,1,2);


% % debug
% figure()
% figSetting()
% view(3)
% surf(Id(2:end,2:end),Iq(2:end,2:end),ldd(2:end,:),'FaceColor','r','EdgeColor','none');
% surf(Id,Iq,Ldd,'FaceColor','none','EdgeColor','b')
%
% figure()
% figSetting()
% view(3)
% surf(Id(2:end,2:end),Iq(2:end,2:end),ldq(:,2:end),'FaceColor','r','EdgeColor','none');
% surf(Id,Iq,Ldq,'FaceColor','none','EdgeColor','b')
%
% figure()
% figSetting()
% view(3)
% surf(Id(2:end,2:end),Iq(2:end,2:end),lqd(2:end,:),'FaceColor','r','EdgeColor','none');
% surf(Id,Iq,Lqd,'FaceColor','none','EdgeColor','b')
%
% figure()
% figSetting()
% view(3)
% surf(Id(2:end,2:end),Iq(2:end,2:end),lqq(:,2:end),'FaceColor','r','EdgeColor','none');
% surf(Id,Iq,Lqq,'FaceColor','none','EdgeColor','b')
%
% % end debug

% output data
Inductance.Id = Id;
Inductance.Iq = Iq;
Inductance.Ldd = Ldd;
Inductance.Ldq = Ldq;
Inductance.Lqd = Lqd;
Inductance.Lqq = Lqq;
end