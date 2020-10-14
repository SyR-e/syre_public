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

function [x,y]=interp_flux_barrier(xxBk,yyBk,yEnd)

x=zeros(size(yyBk,1),50); y=zeros(size(yyBk,1),50);
xAux=zeros(1,50); yAux=zeros(1,50);

for k=1:size(yyBk,1)
    pos=find(yyBk(k,:)<=yEnd(k)*20);
    yAux=linspace(0,yEnd(k),50);
    xAux=interp1(yyBk(k,pos),xxBk(k,pos),yAux);
    x(k,:)=xAux;
    y(k,:)=yAux;
end
end