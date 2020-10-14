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

% mean barrier value:
function [xBkMean,yBkMean]=valore_medio_di_barriera(xxBk,yyBk,yTraf)

xAux=zeros(1,50); yAux=zeros(1,50);
yBkMean=zeros(1,size(yyBk,1));
xBkMean=zeros(1,size(yyBk,1));
for k=1:size(yyBk,1)
%     keyboard
    pos=find(yyBk(k,:)>0 & yyBk(k,:)<=yTraf(k));
    if (length(pos)==1|| isempty(pos))
        pos=find(yyBk(k,:)>=0); 
        yAux=linspace(yyBk(k,pos(1)),yTraf(k),100);
        xAux=interp1(yyBk(k,pos),xxBk(k,pos),yAux);
        yBkMean(k)=mean(yAux);
        xBkMean(k)=interp1(yAux,xAux,mean(yAux));
    else
        yAux=linspace(yyBk(k,pos(1)),yyBk(k,pos(end)),100);
        xAux=interp1(yyBk(k,pos),xxBk(k,pos),yAux);
        yBkMean(k)=mean(yAux);
        xBkMean(k)=interp1(yAux,xAux,mean(yAux));
    end
    
end

end