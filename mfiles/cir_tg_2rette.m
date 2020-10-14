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

% Funzione: circonferenza (incongita) di raggio noto, tangente ad 2 rette (note)
% Input: punti delle rette 
% Output: coordinate centro cerchio, pti di tangenza alle rette
function [x,y,xtan1,ytan1,xtan2,ytan2]=cir_tg_2rette(x1,y1,x2,y2,x3,y3,x4,y4,R)
% keyboard
[a1,b1,c1]=retta_per_2pti(x1,y1,x2,y2);
[a2,b2,c2]=retta_per_2pti(x3,y3,x4,y4);

A1=(a1 - (b1*(a2*(a1^2 + b1^2)^(1/2) - a1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) - b1*(a2^2 + b2^2)^(1/2)))^2;
B1=2*(a1 - (b1*(a2*(a1^2 + b1^2)^(1/2) - a1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) - b1*(a2^2 + b2^2)^(1/2)))*(c1 - (b1*(c2*(a1^2 + b1^2)^(1/2) - c1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) - b1*(a2^2 + b2^2)^(1/2)));
C1=(c1 - (b1*(c2*(a1^2 + b1^2)^(1/2) - c1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) - b1*(a2^2 + b2^2)^(1/2)))^2- R^2*(a1^2 + b1^2);
x1=roots([A1 B1 C1]);
y1=-((sqrt(a2^2+b2^2)*a1-sqrt(a1^2+b1^2)*a2).*x1+(sqrt(a2^2+b2^2)*c1-sqrt(a1^2+b1^2)*c2))/(sqrt(a2^2+b2^2)*b1-sqrt(a1^2+b1^2)*b2);

A2=(a1 - (b1*(a2*(a1^2 + b1^2)^(1/2) + a1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2)+ b1*(a2^2 + b2^2)^(1/2)))^2;
B2=2*(a1 - (b1*(a2*(a1^2 + b1^2)^(1/2) + a1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) + b1*(a2^2 + b2^2)^(1/2)))*(c1 - (b1*(c2*(a1^2 + b1^2)^(1/2) + c1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) + b1*(a2^2 + b2^2)^(1/2)));
C2=(c1 - (b1*(c2*(a1^2 + b1^2)^(1/2) + c1*(a2^2 + b2^2)^(1/2)))/(b2*(a1^2 + b1^2)^(1/2) + b1*(a2^2 + b2^2)^(1/2)))^2-R^2*(a1^2 + b1^2);
x2=roots([A2 B2 C2]);
y2=-((sqrt(a2^2+b2^2)*a1+sqrt(a1^2+b1^2)*a2).*x2+(sqrt(a2^2+b2^2)*c1+sqrt(a1^2+b1^2)*c2))/(sqrt(a2^2+b2^2)*b1+sqrt(a1^2+b1^2)*b2);

x=[x1;x2];
y=[y1;y2];

if(b1==0)
    xtan1=-c1/a1*ones(1,length(x));
    ytan1=y;

else
    xtan1=-(a1*c1/b1^2+a1*y/b1-x).*(1/(1+a1^2/b1^2));
    ytan1=-(a1.*xtan1+c1)./b1;
end

if(b2==0)
    xtan1=-c2/a2*ones(1,length(x));
    ytan1=y;
else
    xtan2=-(a2*c2/b2^2+a2*y/b2-x).*(1/(1+a2^2/b2^2));
    ytan2=-(a2.*xtan2+c2)./b2;

end


