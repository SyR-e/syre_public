% Assegna materiali alle diverse parti della macchina

%%%%%%%%%%%%%
%% Statore
%%%%%%%%%%%%%

[var1,var2]=size(BLKLABELSstat.names.air_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[BLKLABELSstat.xy(kk,1), BLKLABELSstat.xy(kk,2)],cell2mat(BLKLABELSstat.names.air_slot{kk,1}),l,BLKLABELS.materials(2,:),'None',[0, 0, 1],BLKLABELSstat.xy(kk,4));
    indice_air_slot=kk;
end

[var1,var2]=size(BLKLABELSstat.names.Cu_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[BLKLABELSstat.xy(kk+indice_air_slot,1), BLKLABELSstat.xy(kk+indice_air_slot,2)],cell2mat(BLKLABELSstat.names.Cu_slot{kk,1}),l,BLKLABELS.materials(3,:),'None',[0, 0, 1],-1);
    indice_Cu_slot=kk+indice_air_slot;
end
[var1,var2]=size(BLKLABELSstat.names.FeYoke);
for kk=1:var1
    h = MakeComponentMagnet(h,[BLKLABELSstat.xy(kk+indice_Cu_slot,1), BLKLABELSstat.xy(kk+indice_Cu_slot,2)],cell2mat(BLKLABELSstat.names.FeYoke{kk,1}),l,BLKLABELS.materials(4,:),'None',[0, 0, 1],BLKLABELSstat.xy(kk+indice_Cu_slot,4));
end