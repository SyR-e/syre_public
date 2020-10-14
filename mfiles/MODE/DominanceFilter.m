% Copyright 2014
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

function [Fout, Cout]=DominanceFilter(F,C)

Xpop=size(F,1);
Nobj=size(F,2);
Nvar=size(C,2);
k=0;

for xpop=1:Xpop
    dom=0;
    
    for i=1:Xpop
        if F(xpop,:)==F(i,:)
            if xpop >= i
                dom=1;
                break;
            end
        else
            if F(xpop,:)>=F(i,:)
                dom=1;
                break;
            end
        end
    end
    
    if dom==0
        k=k+1;
        F(k,:)=F(xpop,:);
        C(k,:)=C(xpop,:);
    end
end
Fout=F;
Cout=C;