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

function [front,idx] = FastParetoEstimation(x,COST) %pop,info)

pop=[x COST];
[xx,yx]=size(x);
[xC,yC]=size(COST);
info.fitSize=yC;
info.fitIdx = [yx+1:yx+info.fitSize];
[N,aux] = size(pop);
clear aux;

M = info.fitSize;
idx = ones(N,1);
fitP = zeros(1,M);
fitQ = zeros(1,M);
domFlag = 0;

for p=1:N
    if idx(p)
        fitP = pop(p,info.fitIdx);
        for q=p+1:N
            fitQ = pop(q,info.fitIdx);
            
            domFlag = Dominance(fitP,fitQ,0);
            switch domFlag
                case  1 
                    idx(q) = 0;
                case -1 
                    idx(p) = 0;
                otherwise
                        
            end
            
            if ~idx(p)
                break;
            end
        end
    end
end
                  
front = [];
for k=1:N
    if idx(k)
        front = [front;pop(k,:)];
    end
end
