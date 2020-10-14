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

function arco(a,layer,color,tiplin,k,fid)
%% Scrittura dei punti che caratterizzano un arco

if a(k,7)==-1
    b=a(k,:);
    a(k,3:4)=b(5:6);
    a(k,5:6)=b(3:4);
end

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ARC';
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

%definizione del centro
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

%definizione del raggio
r = sqrt((a(k,3)-a(k,1))^2 + (a(k,4)-a(k,2))^2);
%
testo = '40';
fprintf(fid,'%s\n',testo);
testo = num2str(r,'%15.8f');
fprintf(fid,'%s\n',testo);

% Def angolo Iniziale
sinteta = (a(k,4)-a(k,2));
costeta = (a(k,3)-a(k,1));
teta = 180/pi*angle(costeta+1i*sinteta);
%scrittura dell'angolo iniziale
testo = '50';
fprintf(fid,'%s\n',testo);
testo = num2str(teta,'%15.8f');
fprintf(fid,'%s\n',testo);

% Def angolo finale
sinteta = (a(k,6)-a(k,2));
costeta = (a(k,5)-a(k,1));
teta = 180/pi*angle(costeta+1i*sinteta);
% scrittura dell'angolo finale
testo = '51';
fprintf(fid,'%s\n',testo);
testo = num2str(teta,'%15.8f');
fprintf(fid,'%s\n',testo);

