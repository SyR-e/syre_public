
%%%%%%%%%%%
%% Rotore  // assegnazione dinamica come in FEMM non commentata !!funziona
%%%%%%%%%%%
% Barrier
% Mac.Br=0.4;
% CENTRI.materials(7,:)='NMF3F';

% if Br==0 || strcmp(geo.RotType,'SPM')
%     [BarSize,var2]=size(BLKLABELSrot.BarName);
%     for kk=1:BarSize
%         h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],cell2mat(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(2,:),'None',[0, 0, 1],BLKLABELSrot.xy(kk,4));
%     end
% else
%     %Barrier Filled with magnet
%     [BarSize,var2]=size(BLKLABELSrot.BarName);
%     [LabelSize,~]=size(BLKLABELSrot.xy);
%     BarSize=LabelSize-2;
%     for kk=1:LabelSize-2
%         if BLKLABELSrot.xy(kk,3)==6
%             tmpmat=BLKLABELS.materials(7,:);
%         else
%             tmpmat='AIR';
%         end
%         h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],cell2mat(BLKLABELSrot.BarName{kk,1}),l,tmpmat,'Uniform',[BLKLABELSrot.xy(kk,6),BLKLABELSrot.xy(kk,7), BLKLABELSrot.xy(kk,8)],BLKLABELSrot.xy(kk,4));
%     end

for kk=1:length(BLKLABELSrot.xy(:,1))
        switch BLKLABELSrot.xy(kk,3)
            case 1 %air
                h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],cell2mat(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(2,:),'Uniform',[BLKLABELSrot.xy(kk,6),BLKLABELSrot.xy(kk,7), BLKLABELSrot.xy(kk,8)],BLKLABELSrot.xy(kk,4));


            case 5 %rotor iron
                h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],'rotor',l,BLKLABELS.materials(5,:),'None',[0, 0, 1],BLKLABELSrot.xy(kk,4));
                
            case 6 %PM      
                if Br==0 || strcmp(geo.RotType,'SPM')
                   h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],cell2mat(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(2,:),'None',[0, 0, 1],BLKLABELSrot.xy(kk,4));
                else
                   MatName=BLKLABELS.materials(7,:);
                   h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],cell2mat(BLKLABELSrot.BarName{kk,1}),l,BLKLABELS.materials(7,:),'Uniform',[BLKLABELSrot.xy(kk,6),BLKLABELSrot.xy(kk,7), BLKLABELSrot.xy(kk,8)],BLKLABELSrot.xy(kk,4));
                end
            case 7 %shaft               
                h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],'shaft',l,BLKLABELS.materials(8,:),'None',[0, 0, 1],-1);
        end
    
    end
% end

% if strcmp(geo.RotType,'SPM')
%     for jj = kk+1:kk+geo.ps*2
%         h = MakeComponentMagnet(h,[CENTRIrot.xy(jj,1), CENTRIrot.xy(jj,2)],['Rotor_Air_Zone_',num2str(jj-kk)],L_assiale,CENTRI.materials(1,:),'None',[0, 0, 1],-1);
%     end
%     
%     jj = 1;
%     rotor Iron
%     jj=jj+1;
%     h = MakeComponentMagnet(h,[BLKLABELSrot.xy(jj,1), BLKLABELSrot.xy(jj,2)],'rotor',l,BLKLABELS.materials(5,:),'None',[0, 0, 1],BLKLABELSrot.xy(jj,4));
%     
%     Shaft
%     jj=jj+1;
%     h = MakeComponentMagnet(h,[BLKLABELSrot.xy(jj,1), BLKLABELSrot.xy(jj,2)],'shaft',l,BLKLABELS.materials(8,:),'None',[0, 0, 1],-1);
% else
%     rotor Iron
%     kk= kk+1;
%     h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],'rotor',l,BLKLABELS.materials(5,:),'None',[0, 0, 1],BLKLABELSrot.xy(kk,4));
%     
%     Shaft
%     kk= kk+1;
%     h = MakeComponentMagnet(h,[BLKLABELSrot.xy(kk,1), BLKLABELSrot.xy(kk,2)],'shaft',l,BLKLABELS.materials(8,:),'None',[0, 0, 1],-1);
% end

%%%%%%%%%%%%
%% Traferro
%%%%%%%%%%%%
for kk=1:size(BLKLABELStraf.names,1)
    h = MakeComponentMagnet(h,[BLKLABELStraf.xy(kk,1), BLKLABELStraf.xy(kk,2)],BLKLABELStraf.names{kk,1},l,BLKLABELS.materials(BLKLABELStraf.xy(kk,3),:),'None',[0, 0, 1],-1);
end

