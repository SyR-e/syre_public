% Copyright 2014
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

%% Matteo 31/10/2011
%intersection line and circumference:
function [x,y]=intersezione_retta_circonferenza(xc,yc,r,m2,q2)

A=1+m2.^2;
B=2.*m2.*q2-2.*xc-2.*m2.*yc;
C=xc.^2+q2.^2+yc.^2-2.*q2.*yc-r.^2;
% tmp=roots([A B C]);
% NB: if r>0, the output point is the bigger. If r<0, the output point is
% the lower
ksing = r./abs(r);
x = (-B+ksing.*(B.^2-4.*A.*C).^0.5)./(2*A);
% if r>=0
%     x = (-B+(B.^2-4.*A.*C).^0.5)./(2*A);
% else
%     x = (-B-(B.^2-4.*A.*C).^0.5)./(2*A);
% end
y=m2.*x+q2;
