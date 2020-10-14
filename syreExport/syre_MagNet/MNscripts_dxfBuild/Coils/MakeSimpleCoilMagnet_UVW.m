%% 2013/07/26 MG coil construction

clear strato_pos cava_pos strato_neg cava_neg
[strato_pos,cava_pos]=find(avv == coil_number);
[strato_neg,cava_neg]=find(avv == -coil_number);

%% vettore segno dei lati
dh = h.documentHandler;
mh = h.magnetHandler;

u=[];
% fase X, pos
if ~isempty(strato_pos)
    for jj=1:max(size(strato_pos))
if jj==1
invoke(mh,'processCommand',[CGDGV,'.selectObject("slot_',num2str(strato_pos(jj)),'_',num2str(cava_pos(jj)),'", infoSetSelection)']);
else
    invoke(mh,'processCommand',['CALL getDocument().getView().selectObject("slot_',num2str(strato_pos(jj)),'_',num2str(cava_pos(jj)),'",',ITIS,')']);
end
    u=[u;1];
    end
end

% Fase X neg


if ~isempty(strato_neg)
    for jj=1:max(size(strato_neg))
if jj==1
invoke(mh,'processCommand',[CGDGV,'.selectObject("slot_',num2str(strato_neg(jj)),'_',num2str(cava_neg(jj)),'", infoSetSelection)']);
else
    invoke(mh,'processCommand',['CALL getDocument().getView().selectObject("slot_',num2str(strato_neg(jj)),'_',num2str(cava_neg(jj)),'",',ITIS,')']);
end
    u=[u;-1];
    end
end

%% make coil 
%% selezione lati positivi
string=sprintf('REDIM ArrayOfValues(%d)',(max(size(u))-1));
invoke(mh,'processCommand',string);
clear string

if ~isempty(strato_pos)
    for jj=1:max(size(strato_pos))
       string=sprintf('ArrayOfValues(%d)= "slot_%d_%d"',jj-1,strato_pos(jj),cava_pos(jj)); 
    invoke(mh,'processCommand',string);
    
    end
else
    jj=0;
end

%% se ci sono solo lati negativi (es. fase W)
% if isempty(strato_pos)
%     jj = 0;
% end
%% selezione lati negativi

if ~isempty(strato_neg)
    if max(size(strato_neg))==1
        jj = jj +1; %% contatore totale
        string=sprintf('ArrayOfValues(%d)= "slot_%d_%d"',jj-1,strato_neg,cava_neg);
            invoke(mh,'processCommand',string);
    else
        for kk=1:max(size(strato_neg))
            jj = jj +1; %% contatore totale
            string=sprintf('ArrayOfValues(%d)= "slot_%d_%d"',jj-1,strato_neg(kk),cava_neg(kk));
            invoke(mh,'processCommand',string);
        end
    end
end

string=sprintf('%s.makeSimpleCoil(1, ArrayOfValues)',CGD);
invoke(mh,'processCommand',string);
string=sprintf('%s.renameObject("Coil#%g", "%s")',CGD,coil_number,coil_name);
invoke(mh,'processCommand',string);
string=sprintf('%s.setCoilNumberOfTurns("%s", 6)',CGD,coil_name);
invoke(mh,'processCommand',string);

for jj=1:size(u,1)
    if (rem(jj,2) ~= 0) && (sign(u(jj)) == -1)
        string=sprintf('%s.reverseCoilSide("%s", %d)',CGD,coil_name,jj);
        invoke(mh,'processCommand',string);
    end
    if (rem(jj,2) == 0) && (sign(u(jj)) == 1)
        string=sprintf('%s.reverseCoilSide("%s", %d)',CGD,coil_name,jj);
        invoke(mh,'processCommand',string);
    end
end

string=sprintf('%s.selectObject("%s", infoSetSelection)',CGDGV,coil_name);
invoke(mh,'processCommand',string);
string=sprintf('%s.beginUndoGroup("Set %s Properties", true)',CGD,coil_name);
invoke(mh,'processCommand',string);
string=sprintf('REDIM ArrayOfValues(5)');
invoke(mh,'processCommand',string);
string=sprintf('ArrayOfValues(0)= 0');
invoke(mh,'processCommand',string);
string=sprintf('ArrayOfValues(1)= 33.3333');
invoke(mh,'processCommand',string);
string=sprintf('ArrayOfValues(2)= 50');
invoke(mh,'processCommand',string);
string=sprintf('ArrayOfValues(3)= 0');
invoke(mh,'processCommand',string);
string=sprintf('ArrayOfValues(4)= 0');
invoke(mh,'processCommand',string);
string=sprintf('ArrayOfValues(5)= 0');
invoke(mh,'processCommand',string);
string=sprintf('%s.setSourceWaveform("%s","SIN", ArrayOfValues)',CGD,coil_name);
invoke(mh,'processCommand',string);
string=sprintf('%s.endUndoGroup()',CGD);
invoke(mh,'processCommand',string);
