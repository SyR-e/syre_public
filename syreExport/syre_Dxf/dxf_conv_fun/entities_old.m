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

% 28-05-2010 BB
% Modificati gli input per disegnare eventualmente anche l'albero del
% motore -> NO Modifica annullata
function entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria)
% function entities(raggi,avvolgimento,rotore,statore,magneti,fid,LAria,albero)

%% Scrittura dei punti che caratterizzano i diversi layer
%color 0 = nero
%color 1 = rosso
%color 2 = giallo
%color 3 = verde
%color 4 = ciano
%color 5 = blu
%color 6 = magenta

% generazione sezione ENTITIES
% Scrittura intestazione ENTITIES
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'ENTITIES';
fprintf(fid,'%s\n',testo);

% Riconoscimento della matrice (statore)
if ~isempty(statore)
    dimobj = size(statore);
    righe = dimobj(1);
    layer = '0';
    color = '0';
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if statore(k,7) == 0
            linea(statore,layer,color,TipoLinea,k,fid);
        elseif abs(statore(k,7)) == 1
            arco(statore,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

% Riconoscimento della matrice (rotore)
if ~isempty(rotore)
    dimobj = size(rotore);
    righe = dimobj(1);
    layer = '1';
    color = '1';
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if rotore(k,7) == 0
            linea(rotore,layer,color,TipoLinea,k,fid);
        elseif abs(rotore(k,7)) == 1
            arco(rotore,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

% Riconoscimento della matrice (magneti)
if ~isempty(magneti)
    dimobj = size(magneti);
    righe = dimobj(1);
    layer = '2';
    color = '3';
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if magneti(k,7) == 0
            linea(magneti,layer,color,TipoLinea,k,fid);
        elseif abs(magneti(k,7)) == 1
            arco(magneti,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end


%riconoscimento della matrice (Rame)
if ~isempty(avvolgimento)
    dimobj = size(avvolgimento);
    righe = dimobj(1);
    layer = '3';
    color = '4';
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if avvolgimento(k,7) == 0
            linea(avvolgimento,layer,color,TipoLinea,k,fid);
        elseif abs(avvolgimento(k,7)) == 1
            arco(avvolgimento,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

%riconoscimento della matrice (Raggi)
if ~isempty(raggi)
    dimobj = size(raggi);
    righe = dimobj(1);
    layer = '4';
    color = '5';
    TipoLinea = 'Continuous';
    
    for k = 1:righe
        if raggi(k,7) == 0
            linea(raggi,layer,color,TipoLinea,k,fid);
        elseif abs(raggi(k,7)) == 1
            arco(raggi,layer,color,TipoLinea,k,fid);
        else
            break
        end
    end
end

% 28-05-2010 BB
% Gli input da raggi (il primo) a fid (l'ultimo prima degli eventuali
% LAria e albero) sono 6. Quindi devo controllare se le matrici LAria e albero
% sono oppure no vuote se gli input sono più di sei

% NO -> Modifica annullata

if nargin>7
% if nargin>6
    %riconoscimento della matrice (LAria)
    if ~isempty(LAria)
        dimobj = size(LAria);
        righe = dimobj(1);
        layer = '5';
        color = '6';
        TipoLinea = 'Continuous';
        
        for k = 1:righe
            if LAria(k,7) == 0
                linea(LAria,layer,color,TipoLinea,k,fid);
            elseif abs(raggi(k,7)) == 1
                arco(LAria,layer,color,TipoLinea,k,fid);
            else
                break
            end
        end
    end
    % 28-05-2010 BB
    % Modificati gli input per disegnare eventualmente anche l'albero del
    % motore
    %riconoscimento della matrice (albero)
    % if nargin==8 && (~isempty(albero))
    %     dimobj = size(albero);
    %     righe = dimobj(1);
    %     layer = '6';
    %     color = '2';
    %     TipoLinea = 'Continuous';
        
    %     for k = 1:righe
    %         if albero(k,7) == 0
    %             linea(albero,layer,color,TipoLinea,k,fid);
    %         elseif abs(albero(k,7)) == 1
    %             arco(albero,layer,color,TipoLinea,k,fid);
    %         else
    %             break
    %         end
    %     end
    % end
end

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);


