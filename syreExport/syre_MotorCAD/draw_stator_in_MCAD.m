function draw_stator_in_MCAD(mcad,geo,~,~,~)

if geo.parallel_slot==0
    invoke(mcad,'SetVariable','SlotType',0);    %slot type (ParallelTooth)
else
    invoke(mcad,'SetVariable','SlotType',2);    %slot type (ParallelSlot)
end

Q = 6*geo.p*geo.q;
invoke(mcad,'SetVariable','Slot_number',Q);

tmp = 2*geo.R;
invoke(mcad,'SetVariable','Stator_Lam_Dia',tmp);

tmp = geo.l;
invoke(mcad,'SetVariable','Stator_Lam_Length',tmp);
invoke(mcad,'SetVariable','Rotor_Lam_Length',tmp);
invoke(mcad,'SetVariable','Magnet_Length',tmp);

invoke(mcad,'SetVariable','Motor_Length',tmp+80);

tmp = 2*geo.R+20;
invoke(mcad,'SetVariable','Housing_Dia',tmp);

tmp = 2*(geo.r+geo.g);
invoke(mcad,'SetVariable','Stator_Bore',tmp);

tmp = geo.g;
tmp = num2str(tmp);
tmp(tmp=='.') = ',';
invoke(mcad,'SetVariable','Airgap',tmp);

if geo.parallel_slot==0
    tmp=geo.wt;
    invoke(mcad,'SetVariable','Tooth_Width',tmp);      %ParallelTooth
else
    tmp=(geo.r+geo.g+geo.lt/15)*sin(pi/geo.p/geo.Qs)-geo.wt;
    invoke(mcad,'SetVariable','Slot_Width',tmp);       %ParallelSlot
end

tmp = geo.lt;
invoke(mcad,'SetVariable','Slot_Depth',tmp);
if geo.parallel_slot==0
    tmp=num2str(geo.SFR);
    invoke(mcad,'SetVariable','Slot_Corner_Radius',tmp); %ParallelTooth
end

tmp = num2str(geo.ttd);
tmp(tmp=='.') = ',';
invoke(mcad,'SetVariable','Tooth_Tip_Depth',tmp);

tmp = (geo.acs)*((geo.r+geo.g)*2*pi/Q);
tmp = num2str(tmp);
tmp(tmp=='.') = ',';
invoke(mcad,'SetVariable','Slot_Opening',tmp);

tmp = num2str(geo.tta);
tmp(tmp=='.') = ',';
invoke(mcad,'SetVariable','Tooth_Tip_Angle',tmp);
invoke(mcad,'SetVariable','Sleeve_Thickness','0');

%% Wedge
if geo.nmax>9999
    invoke(mcad,'SetVariable','Wedge_Model', '0');  %% with wedge
else
    invoke(mcad,'SetVariable','Wedge_Model', '1');  %% without wedge
end

end