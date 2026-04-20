%% autoSyRe_TimingTest_v3.m

close all
clear
clc

[filename,pathname,~] = uigetfile([cd '\.mat'],'Load motor model');

mot = load([pathname filename]);

answer = inputdlg('Number of workers?','Setup',1,{'8'});
numWorker = eval(answer{1});

% numWorker = 16;

solver = 'XFEMM';

setup.filename  = filename;
setup.pathname  = pathname;
setup.mot       = mot;
setup.numWorker = numWorker;
setup.solver    = solver;


date=fix(clock);
year=date(1);
month=date(2);
day=date(3);
hour=date(4);
mins=date(5);
logname = ['LogFile_TimingTest_v3_', num2str(year) '_' num2str(month) '_' num2str(day) '_' num2str(hour) 'h_' num2str(mins) '.txt'];

logfile = [pathname logname];

logHeader(logfile,'SyR-e Timing Test v3')

fid = fopen(logfile,'a');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Timing test details:\r\n');
fprintf(fid,['Number of cores   : ' int2str(feature('numcores')) '\r\n']);
fprintf(fid,['Number of workers : ' int2str(numWorker) '\r\n']);
fprintf(fid,['Solver            : ' solver '\r\n']);
fprintf(fid,['Pathname          : ' strrep(pathname,'\','\\') '\r\n']);
fprintf(fid,['Motor model       : ' filename(1:end-4) '\r\n']);
fclose(fid);
edit(logfile)

%% test description

fid = fopen(logfile,'a');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Timing test v3 description:\r\n');
fprintf(fid,'Parallel pool enabled:\r\n');
fprintf(fid,'p1)  Flux Map    : flux and torque map (30 rot pos), 16x16 grid -->  30 x 256 FEMM simulations\r\n');
fprintf(fid,'p2)  Loss Map    : iron loss map (90 rot pos), 8x8 grid         -->  90 x  64 FEMM simulations\r\n');
fprintf(fid,'p3)  AC Loss Map : AC loss map, time-harmonic simulations 1+512 -->  1+1x512  FEMM simulations\r\n');
fprintf(fid,'p4)  Effy Map    : efficiency map, 128x128 grid                 -->   0 x  16384 FEMM simulations\r\n');
fprintf(fid,'s1)  Single point                                               -->  30 x 1  FEMM simulations\r\n');
fprintf(fid,'s2)  Single point                                               -->  90 x 1  FEMM simulations\r\n');
fprintf(fid,'s3)  Single point + iron loss                                   -->  90 x 1  FEMM simulations\r\n');
fprintf(fid,'s4)  Single point + NVH                                         -->  90 x 1  FEMM simulations\r\n');
fprintf(fid,'s5)  Demagnetization limit @ 1 temperature                      -->  ?? x 1  FEMM simulations\r\n');
fprintf(fid,'s6)  Structural simulation (10 repetition)                      -->  10 x 1  FEMM simulations\r\n');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'\r\n\r\n');
fclose(fid);
edit(logfile)

logMessage(logfile,'Parallel pool check...');

% fid=fopen(logfile,'a');
% fprintf(fid,[datestr(now) '-->']);
% fprintf(fid,['Parallel pool check...\r\n']);
% fclose(fid);
% edit(logfile)

% parallel pool creation
if ~isempty(isprop(gcp('nocreate'),'NumWorkers')) % no parallel pool enabled
    logMessage(logfile,'Parallel pool enabled. Shutting down...');
    
    delete(gcp)
    
    logMessage(logfile,'Parallel pool disabled');
end

parpool(numWorker);
logMessage(logfile,['Parallel pool enabled with ' int2str(numWorker) ' workers']);

if strcmp(solver,'FEMM')
    mot.dataSet.XFEMMsimulation = 0;
elseif strcmp(solver,'XFEMM')
    mot.dataSet.XFEMMsimulation = 1;
else
    error('wrong solver setup')
end


%% flux map test

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singm';
mot.dataSet.NumGrid          = 16;
mot.dataSet.CurrLoPP         = 3;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 30;
mot.dataSet.AngularSpanPP    = 60;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = interp1(mot.mat.LayerMag.temp.temp,mot.mat.LayerMag.temp.Br,mot.dataSet.tempPP);

logMessage(logfile,'Flux map: start simulations');

tic
eval_fluxMap(mot.dataSet);
compTime = toc;
close all

logMessage(logfile,['Flux map: evaluated in ' num2str(compTime) ' s']);
results.fluxMap.time    = compTime;
results.fluxMap.sim     = mot.dataSet.NumGrid^2;
results.fluxMap.FEAcall = mot.dataSet.NumOfRotPosPP*mot.dataSet.NumGrid^2;

logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')

%% iron loss map

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singmIron';
mot.dataSet.NumGrid          = 8;
mot.dataSet.CurrLoPP         = 3;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 90;
mot.dataSet.AngularSpanPP    = 180;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = interp1(mot.mat.LayerMag.temp.temp,mot.mat.LayerMag.temp.Br,mot.dataSet.tempPP);

logMessage(logfile,'Loss map: start simulations');

tic
eval_fluxMap(mot.dataSet);
compTime = toc;
close all

logMessage(logfile,['Loss map: evaluated in ' num2str(compTime) ' s']);
results.lossMap.time    = compTime;
results.lossMap.sim     = mot.dataSet.NumGrid^2;
results.lossMap.FEAcall = mot.dataSet.NumOfRotPosPP*mot.dataSet.NumGrid^2;

logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')

%% AC loss

logMessage(logfile,'Slot model creation...')
mot.dataSet.XFEMMcreation = 0;
skinEffect_draw(mot.dataSet);
logMessage(logfile,'Slot model created!')

mot.dataSet.SlotConductorFrequency = [1 5 10 20 50 100 150 200 225 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 4000 8000 10000 15000 20000];
mot.dataSet.SlotConductorTemperature = [-40 -20 0 20 30 40 50 60 80 100 120 140 160 180 200 220];

logMessage(logfile,'AC loss: start simulations...')
ACresults = skinEffect_eval(mot.dataSet);
close all
logMessage(logfile,['AC loss: DC simulation done in ' num2str(ACresults.timeDC) ' s'])
logMessage(logfile,['AC loss: AC simulations done in ' num2str(ACresults.timeAC) ' s'])

results.ACloss.time    = ACresults.timeAC+ACresults.timeDC;
results.ACloss.sim     = 1+numel(mot.dataSet.SlotConductorFrequency)*numel(mot.dataSet.SlotConductorTemperature);
results.ACloss.FEAcall = 1+numel(mot.dataSet.SlotConductorFrequency)*numel(mot.dataSet.SlotConductorTemperature);
results.ACloss.timeAC  = ACresults.timeAC;
results.ACloss.timeDC  = ACresults.timeDC;

logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')



%% Efficiency map

mot.motorModel.TnSetup.nmin           = 0;
mot.motorModel.TnSetup.nmax           = 18000;
mot.motorModel.TnSetup.nstep          = 128;
mot.motorModel.TnSetup.Tmin           = 0;
mot.motorModel.TnSetup.Tmax           = 550;
mot.motorModel.TnSetup.Tstep          = 128;
mot.motorModel.TnSetup.IronLossFlag   = 'Yes';
mot.motorModel.TnSetup.PMLossFlag     = 'Yes';
mot.motorModel.TnSetup.SkinEffectFlag = 'Yes';
mot.motorModel.TnSetup.temperature    = 120;
mot.motorModel.TnSetup.Control        = 'MTPA';

logMessage(logfile,'Effy map: start evaluation');

tic
MMM_MaxTw(mot.motorModel,[],0);
compTime = toc;
close all

logMessage(logfile,['Effy map: evaluated in ' num2str(compTime) ' s']);
results.effyMap.time    = compTime;
results.effyMap.sim     = mot.motorModel.TnSetup.Tstep*mot.motorModel.TnSetup.nstep;
results.effyMap.FEAcall = 0;

%% Parallel pool closing
logMessage(logfile,'Parallel pool disabling...')
delete(gcp)
logMessage(logfile,'Parallel pool disabled')

%% Single point simulation

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singt';
mot.dataSet.CurrLoPP         = 1;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 30;
mot.dataSet.AngularSpanPP    = 60;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = interp1(mot.mat.LayerMag.temp.temp,mot.mat.LayerMag.temp.Br,mot.dataSet.tempPP);
mot.dataSet.GammaPP          = 135;
mot.dataSet.axisType         = 'PM';

logMessage(logfile,['Single point (30pos): start simulations']);
tic
eval_operatingPoint(mot.dataSet,0)
compTime = toc;
close all

logMessage(logfile,['Single point (30pos) evaluated in ' num2str(compTime) ' s']);
results.singt1.time    = compTime;
results.singt1.sim     = 1;
results.singt1.FEAcall = mot.dataSet.NumOfRotPosPP;


logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')

%% Single point simulation

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singt';
mot.dataSet.CurrLoPP         = 1;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 90;
mot.dataSet.AngularSpanPP    = 180;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = interp1(mot.mat.LayerMag.temp.temp,mot.mat.LayerMag.temp.Br,mot.dataSet.tempPP);
mot.dataSet.GammaPP          = 135;
mot.dataSet.axisType         = 'PM';

logMessage(logfile,['Single point (90pos): start simulations']);
tic
eval_operatingPoint(mot.dataSet,0)
compTime = toc;
close all

logMessage(logfile,['Single point (90pos) evaluated in ' num2str(compTime) ' s']);
results.singt.time    = compTime;
results.singt.sim     = 1;
results.singt.FEAcall = mot.dataSet.NumOfRotPosPP;


logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')

%% Single point - iron loss

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singtIron';
mot.dataSet.CurrLoPP         = 1;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 90;
mot.dataSet.AngularSpanPP    = 180;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = interp1(mot.mat.LayerMag.temp.temp,mot.mat.LayerMag.temp.Br,mot.dataSet.tempPP);
mot.dataSet.GammaPP          = 135;
mot.dataSet.axisType         = 'PM';

logMessage(logfile,['Single point + iron loss: start simulations']);
tic
eval_operatingPoint(mot.dataSet,0)
compTime = toc;
close all

logMessage(logfile,['Single point + iron loss evaluated in ' num2str(compTime) ' s']);
results.singtIron.time    = compTime;
results.singtIron.sim     = 1;
results.singtIron.FEAcall = mot.dataSet.NumOfRotPosPP;


logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')

%% Single point - NVH

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'force';
mot.dataSet.CurrLoPP         = 1;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 90;
mot.dataSet.AngularSpanPP    = 180;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = interp1(mot.mat.LayerMag.temp.temp,mot.mat.LayerMag.temp.Br,mot.dataSet.tempPP);
mot.dataSet.GammaPP          = 135;
mot.dataSet.axisType         = 'PM';

logMessage(logfile,['Single point + NVH sources: start simulations']);
tic
eval_operatingPoint(mot.dataSet,0)
compTime = toc;
close all

logMessage(logfile,['Single point + NVH sources evaluated in ' num2str(compTime) ' s']);
results.singtNVH.time    = compTime;
results.singtNVH.sim     = 1;
results.singtNVH.FEAcall = mot.dataSet.NumOfRotPosPP;


logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')

%% Demagnetization curve

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'demagCurve';
mot.dataSet.CurrLoPP         = 1;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 90;
mot.dataSet.AngularSpanPP    = 180;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 80;
mot.dataSet.BrPP             = mot.mat.LayerMag.temp.Br;
mot.dataSet.GammaPP          = 135;
mot.dataSet.axisType         = 'PM';

logMessage(logfile,['Demagnetization limit: start simulations']);
tic
[demagLimit,~] = eval_idemag_vol(mot.dataSet);
compTime = toc;
close all

iterations = numel(demagLimit.iterations.temp_1.Iiter);

logMessage(logfile,['Demagnetization limit evaluated in ' num2str(compTime) ' s with ' int2str(iterations) ' iterations']);
results.demagnetization.time    = compTime;
results.demagnetization.sim     = iterations;
results.demagnetization.FEAcall = 2*iterations;


logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')



%% Structural simulation

nVect = 2000:2000:20000;

results.structural.time    = 0;
results.structural.sim     = numel(nVect);
results.structural.FEAcall = numel(nVect);
results.structural.timeMin = -1;
results.structural.timeMax = 1e10;

for ii=1:length(nVect)

    mot.dataSet.EvalType  = 'structural';
    mot.dataSet.EvalSpeed = nVect(ii);

    logMessage(logfile,['Structural (PDE) simulation #' int2str(ii) ' : start simulations']);

    tic
    eval_vonMisesStress(mot.dataSet);
    compTime = toc;
    close all

    logMessage(logfile,['Structural (PDE) simulation #' int2str(ii) ' : evaluated in ' num2str(compTime) ' s']);
    results.structural.time = results.structural.time+compTime;
    if results.structural.timeMin>compTime
        results.structural.timeMin = compTime;
    end
    if results.structural.timeMax<compTime
        results.structural.timeMax = compTime;
    end
end

logMessage(logfile,['Structural simulations: ' int2str(numel(nVect)) ' simulations done in ' num2str(results.structural.time) ' s'])

logMessage(logfile,'Deleting temporary files...');
clearTmpFolder(0);
logMessage(logfile,'Temporary files deleted')


%% end test

save([logfile(1:end-4) '.mat'],'results');
logMessage(logfile,['results saved in: ' strrep([logfile(1:end-4) '.mat'],'\','\\')])

logMessage(logfile,'Timing test done!!!');







