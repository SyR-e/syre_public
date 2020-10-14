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

function [maxsegdeg,raggio,ang1,ang]=Disegna_Arco(dati,ris)
% Disegna_Arco(dati,ris)
% input:
% - dati: una riga della matrice statore o rotore (xy punto 1, xy punto 2, xy centro, verso
% rotazione)
% 21-06-2010 BB
% - ris: risoluzione con cui devo disegnare l'arco [mm]

% mi_drawarc(x1,y1,x2,y2,angle,maxseg). maxseg->numero massimo di segmenti
% in cui divido l'arco da disegnare.


% centro
x0 = dati(1); y0 = dati(2);
if dati(7) > 0
    % punto 1
    x1 = dati(3); y1 = dati(4);
    %  punto 2
    x2 = dati(5); y2 = dati(6);
else
    % punto 1
    x2 = dati(3); y2 = dati(4);
    %  punto 2
    x1 = dati(5); y1 = dati(6);
end

% apertura angolare
ang1 = 180/pi * atan2((y1-y0),(x1-x0));
ang2 = 180/pi * atan2((y2-y0),(x2-x0));

ang=ang2-ang1;
if ang<0
    ang=ang+360;
end

% 21-06-2010 BB
% Calcolo il numero di segmenti in cui spezzo il disegno dell'arco in modo
% da rispettare il valore si risoluzione inserito in ingresso
raggio=sqrt((x0-x1)^2+(y0-y1)^2);
maxseg=ang*pi/180*raggio/ris;
maxsegdeg=0.25*(ang/maxseg);
% 2013/07/05 MG modifica al fine di evitare che circonferenze troppo
% piccole abbiano troppo pochi punti, non importa se viene meno la costanza
% della risoluzione
if (maxsegdeg>20)
    maxsegdeg=20;
end
if isnan(maxsegdeg)
    maxsegdeg=20;
end

% keyboard
% disegna arco
% mi_drawarc(x1,y1,x2,y2,ang,ris);

% mi_addnode(x1,y1);mi_addnode(x2,y2)
mi_drawarc(x1,y1,x2,y2,ang,maxsegdeg);




