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

function DXFconv(raggi,avvolgimento,rotore,statore,magneti,nomefile,LAria)

% Scrittura del file dxf per motori
% versione 1.0
% DXFconv(rotore,statore,magneti,nomefile)
%
% RAGGI descrizione dei raggi di limitazione del disegno
% AVVOLGIMENTO descrizione dei punti dell'avvolgimento
% ROTORE    descrizione dei punti di rotore
% STATORE   descrizione dei punti di statore
% MAGNETI   descrizione dei punti dei magneti
% NOMEFILE  nome del file completo di path (system indipendet)

% Apertura file di destinazione
fid = fopen(nomefile,'w');

% Definizione sezione HEADER
header(fid)

% Definizione sezione CLASSES
classes(fid)

% Definizione sezione TABLES
tables(fid)

% Definizione sezione BLOCKS
blocks(fid)

% Definizione sezione ENTITIES
if nargin<7   
    entities(raggi,avvolgimento,rotore,statore,magneti,fid);
else
    entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria);
end

% Definizione sezione OBJECTS
objects(fid)

% Chiusura del file
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'EOF';
fprintf(fid,'%s\n',testo);

% Chiusura file di lavoro
fclose(fid);


