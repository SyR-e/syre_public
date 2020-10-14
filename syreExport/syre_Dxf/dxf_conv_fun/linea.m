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

function linea(a,layer,color,tiplin,k,fid)
%% Scrittura dei punti che caratterizzano una linea
 
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'LINE';
fprintf(fid,'%s\n',testo);

%definisco il Layer
testo = '8';
fprintf(fid,'%s\n',testo);
testo = layer;
fprintf(fid,'%s\n',testo);

%definisco il tipo di linea
testo = '6';
fprintf(fid,'%s\n',testo);
testo = tiplin;
fprintf(fid,'%s\n',testo);

%definisco il colore
testo = '62';
fprintf(fid,'%s\n',testo);
testo = color;
fprintf(fid,'%s\n',testo);

%definizione del primo punto
testo = '10';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,1),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '20';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,2),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '30';
fprintf(fid,'%s\n',testo);
testo = '0';
fprintf(fid,'%s\n',testo);

%definizione del secondo punto
testo = '11';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,3),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '21';
fprintf(fid,'%s\n',testo);
testo = num2str(a(k,4),'%15.8f');
fprintf(fid,'%s\n',testo);
testo = '31';
fprintf(fid,'%s\n',testo);
testo = '0';
fprintf(fid,'%s\n',testo);

