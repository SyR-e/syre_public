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

%% evaluateF.m
% J  [OUT] : The objective Vector. J is a matrix with as many rows as
%            trial vectors in X and as many columns as objectives.
% X   [IN] : Decision Variable Vector. X is a matrix with as many rows as
%            trial vector and as many columns as decision variables.
% Dat [IN] : Parameters defined in NNCparam.m
%
%% Main call
function [J,G]=evaluateF(X,Dat)

%numberOfIndividuals = Dat.XPOP;
numberOfIndividuals = size(X,1);

J=zeros(numberOfIndividuals,Dat.NOBJ);
f=Dat.CostProblem;

    parfor i=1:numberOfIndividuals
        disp(['Evaluating ' num2str(i) '/' num2str(numberOfIndividuals) ' solution.'])
        [J(i,:),G{i}]=f(X(i,:)');
    end
end

%% DTLZ2 Benchmark function. Defined in:
% K. Deb, L. Tiele, M. Laummans, and E. Zitzler. Scalable test problems
% for evolutionary multi-objective optimization. Institut fur Technische
% Informatik und Kommunikationsnetze, ETH Zurich, Tech. Rep. TIK-Technical
% Report No. 112, Feb. 2001.
function J=DTLZ2(X,Dat)

Xpop=size(X,1);
Nvar=Dat.NVAR;
M=Dat.NOBJ;
K=Nvar+1-M;
J=ones(Xpop,M);

for xpop=1:Xpop
    Gxm=(X(xpop,M:Nvar)-0.5*ones(1,K))*(X(xpop,M:Nvar)-0.5*ones(1,K))';
    Cos=cos(X(xpop,1:M-1)*pi/2);
    
    J(xpop,1)=prod(Cos)*(1+Gxm);
    for nobj=1:M-1
        J(xpop,nobj+1)=(J(xpop,1)/prod(Cos(1,M-nobj:M-1)))...
            *sin(X(xpop,M-nobj)*pi/2);
    end
end
end