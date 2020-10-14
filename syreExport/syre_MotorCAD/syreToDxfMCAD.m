% Copyright 2015
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



function syreToDxfMCAD(pathname,filename)

% syreToDxf.m - exports a fem model created by syre to dxf
% input: motorname.mat (created by syre along with motorname.fem)
% output: motorname.dxf, into the folder motorname

%  if nargin<4
%     close all; clc;
%    load LastPath.mat
    load([pathname filename]);
    
%     geo.statordxf=geo.stator;
%     
%     xi=geo.R; yi=0; angle=pi/(geo.p*geo.Qs);
%     xos = real((xi + 1j*yi).* exp(1j*angle));
%     yos = imag((xi + 1j*yi).* exp(1j*angle));
%     geo.statordxf(2,5)=xos;    geo.statordxf(2,6)=yos;
%     
%     xi=geo.r+geo.g; yi=0; angle=pi/(geo.p*geo.Qs);
%     xor= real((xi + 1j*yi).* exp(1j*angle));
%     yor= imag((xi + 1j*yi).* exp(1j*angle));
%     geo.statordxf(4,5)=xi;    geo.statordxf(4,6)=yi;
%     geo.statordxf(4,3)=xor;    geo.statordxf(4,4)=yor;
%      
%     geo.statordxf(3,1)=xos;    geo.statordxf(3,2)=yos;
%     geo.statordxf(3,3)=xor;    geo.statordxf(3,4)=yor;

 
   stator = geo.stator(1:17,:);
   stator(:,1:6)=0;
   rotor =  geo.rotor;
%   magneti = geo.Mag;
% end

%%%%%%%%
% export to dxf

pathname_DXF=pathname;

if not(isfolder(pathname_DXF))
    mkdir(pathname_DXF);
end
% nomefile=[pathname_DXF,filename(1:end-4),'.dxf'];
raggi=[];
avvolgimento=[];
magneti=[];

DXFconv(raggi,avvolgimento,rotor,stator,magneti,[pathname_DXF filename(1:end-4),'.dxf']);
