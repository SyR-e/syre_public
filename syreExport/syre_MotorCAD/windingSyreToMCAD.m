function windingSyreToMCAD(mcad,pathname,filename,file_mot)
load([pathname filename])

geo.avvtot = geo.win.avv;
cyclew=1;
for k=2:1:(geo.p*2)
    if cyclew==1;
        geo.avvtot=[geo.avvtot (-geo.win.avv)];
        cyclew=0;
    else  geo.avvtot=[geo.avvtot geo.win.avv];
        cyclew=1;
    end
end

i=1; m=1; e=1; n=1; j=1; o=1;
for k=1:1:(geo.Qs*2*geo.p)
    value=geo.avvtot(1,k);
    if value==1;
        ph1go(i)=k; i=i+1;
    end
    
    if value==-1;
        ph1ret(m)=k; m=m+1;
    end
    
    if value==2;
        ph2go(n)=k; n=n+1;
    end
    
    if value==-2;
        ph2ret(o)=k; o=o+1;
    end
    
    if value==3;
        ph3go(j)=k; j=j+1;
    end
    
    if value==-3;
        ph3ret(e)=k; e=e+1;
    end
    
end

nbob=(e-1)*2;
geo.NbobInteger=round(geo.win.Nbob);

invoke(mcad,'LoadFromFile',[pathname file_mot]);

invoke(mcad,'SetVariable','MagWindingType', 1);
invoke(mcad,'SetVariable','MagPathType', 1);
invoke(mcad,'SetVariable','NumberOfCoils', nbob);
invoke(mcad,'SetVariable','Coil_Divider_Width',0);

a=1;
for i=0:2:(nbob-1)
    invoke(mcad,'SetVariable',['Phase_1_Go1[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_1_Go2[', num2str(i), ']'], ph1go(a));
    invoke(mcad,'SetVariable',['Phase_1_Return1[', num2str(i), ']'], ph1ret(a));
    invoke(mcad,'SetVariable',['Phase_1_Return2[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_1_Turns[', num2str(i), ']'], geo.NbobInteger);
    
    
    i=i+1;
    
    invoke(mcad,'SetVariable',['Phase_1_Go1[', num2str(i), ']'], ph1go(a));
    invoke(mcad,'SetVariable',['Phase_1_Go2[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_1_Return1[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_1_Return2[', num2str(i), ']'], ph1ret(a));
    invoke(mcad,'SetVariable',['Phase_1_Turns[', num2str(i), ']'], geo.NbobInteger);
    
    a=a+1;
end

a=1;
for i=0:2:(nbob-1)
    
    invoke(mcad,'SetVariable',['Phase_2_Go1[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_2_Go2[', num2str(i), ']'], ph2go(a));
    invoke(mcad,'SetVariable',['Phase_2_Return1[', num2str(i), ']'], ph2ret(a));
    invoke(mcad,'SetVariable',['Phase_2_Return2[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_2_Turns[', num2str(i), ']'], geo.NbobInteger);
    
    i=i+1;
    
    invoke(mcad,'SetVariable',['Phase_2_Go1[', num2str(i), ']'], ph2go(a));
    invoke(mcad,'SetVariable',['Phase_2_Go2[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_2_Return1[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_2_Return2[', num2str(i), ']'], ph2ret(a));
    invoke(mcad,'SetVariable',['Phase_2_Turns[', num2str(i), ']'], geo.NbobInteger);
    
    a=a+1;
end

a=1;
for i=0:2:(nbob-1)
    invoke(mcad,'SetVariable',['Phase_3_Go1[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_3_Go2[', num2str(i), ']'], ph3go(a));
    invoke(mcad,'SetVariable',['Phase_3_Return1[', num2str(i), ']'], ph3ret(a));
    invoke(mcad,'SetVariable',['Phase_3_Return2[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_3_Turns[', num2str(i), ']'], geo.NbobInteger);
    
    i=i+1;
    
    invoke(mcad,'SetVariable',['Phase_3_Go1[', num2str(i), ']'], ph3go(a));
    invoke(mcad,'SetVariable',['Phase_3_Go2[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_3_Return1[', num2str(i), ']'], 0);
    invoke(mcad,'SetVariable',['Phase_3_Return2[', num2str(i), ']'], ph3ret(a));
    invoke(mcad,'SetVariable',['Phase_3_Turns[', num2str(i), ']'], geo.NbobInteger);
    a=a+1;
end


if strcmp(geo.win.condType,'Square')
    %Hairping winding
    invoke(mcad,'SetVariable','Armature_CoilStyle', 1);
    
    %Lap Winding type (no custom)
    invoke(mcad,'SetVariable','MagneticWindingType', 0);
    
    %Conductor dimensions
    invoke(mcad,'SetVariable','Copper_Width', geo.win.wCond);
    invoke(mcad,'SetVariable','Copper_Height', geo.win.hCond);
    invoke(mcad,'SetVariable','Copper_Corner_Radius', geo.win.rCond);
    
    invoke(mcad,'SetVariable','Insulation_Thickness', geo.win.condIns);
    
    %Number of conductors in a slot
    invoke(mcad,'SetVariable','WindingLayers', geo.win.nCond);
    
    %Conductor in parallel
    tmp = geo.win.Ns/dataSet.SlotConductorNumber/geo.q;
    tmp = num2str(tmp);
    tmp(tmp=='.') = ',';
    invoke(mcad,'SetVariable','ParallelPaths', tmp);
    
    throw =  geo.q * 3 * geo.win.n3phase;
    tmp = throw;
    tmp = num2str(tmp);
    tmp(tmp=='.') = ',';
    invoke(mcad,'SetVariable','MagThrow', tmp);
    
    tmp = throw * (geo.r+geo.g+geo.lt/15)*sin(pi/geo.p/geo.Qs)/2;
    tmp = num2str(tmp);
    tmp(tmp=='.') = ',';
    invoke(mcad,'SetVariable','EWdg_Overhang_[F]', tmp);
    invoke(mcad,'SetVariable','EWdg_Overhang_[R]', tmp);
end
end