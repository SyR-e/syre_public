%% autoSyRe_TimingTest_v2.m

close all
clear
clc

[filename,pathname,~] = uigetfile([cd '\.mat'],'Load motor model');

mot = load([pathname filename]);

answer = inputdlg('Number of workers?','Setup',1,{'8'});
numWorker = eval(answer{1});


setup.filename  = filename;
setup.pathname  = pathname;
setup.mot       = mot;
setup.numWorker = numWorker;


date=fix(clock);
year=date(1);
month=date(2);
day=date(3);
hour=date(4);
mins=date(5);
logname = ['LogFile_TimingTest_v2_', num2str(year) '_' num2str(month) '_' num2str(day) '_' num2str(hour) 'h_' num2str(mins) '.txt'];

logfile = [pathname logname];

logHeader(logfile,'SyR-e Timing Test v2')

fid = fopen(logfile,'a');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Timing test details:\r\n');
fprintf(fid,['Number of cores   : ' int2str(feature('numcores')) '\r\n']);
fprintf(fid,['Number of workers : ' int2str(numWorker) '\r\n']);
fprintf(fid,['Pathname          : ' strrep(pathname,'\','\\') '\r\n']);
fprintf(fid,['Motor model       : ' filename(1:end-4) '\r\n']);
fclose(fid);
edit(logfile)

%% test description

fid = fopen(logfile,'a');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Timing test v2 description:\r\n');
fprintf(fid,'1)  FEAfix1  : 1 FEAfix simulation (6 rot pos), 2 op points  -->    12 FEMM simulations\r\n');
fprintf(fid,'2)  FEAfix4  : 4 FEAfix simulation (6 rot pos), 2 op points  -->    48 FEMM simulations\r\n');
fprintf(fid,'3)  FEAfix5  : 5 FEAfix simulation (6 rot pos), 2 op points  -->    60 FEMM simulations\r\n');
fprintf(fid,'4)  FEAfix8  : 8 FEAfix simulation (6 rot pos), 2 op points  -->    96 FEMM simulations\r\n');
fprintf(fid,'5)  FEAfix16 : 16 FEAfix simulation (6 rot pos), 2 op points -->   192 FEMM simulations\r\n');
fprintf(fid,'6)  FEAfix16 : 16 FEAfix simulation (6 rot pos), MTPA search --> ~3840 FEMM simulations\r\n');
fprintf(fid,'7)  Flux Map : flux and torque map (30 rot pos), 15x15 grid  -->  6750 FEMM simulations\r\n');
fprintf(fid,'8)  Loss Map : iron loss map (90 rot pos), 9x9 grid          -->  7290 FEMM simulations\r\n');
fprintf(fid,'9)  Demag.   : demagnetization area (2 rot pos)              -->     2 FEMM simulations\r\n');
fprintf(fid,'10) Struct.  : PDE structural analysis (Matlab PDE)          -->     1  PDE simulations\r\n');
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

%% FEAfix test: FEAfix4, FEAfix5, FEAfix8 and FEAfix16 (without and with gammaFix)

nVector  = [1 4 5 8 16 16];
gfVector = [0 0 0 0  0  1];

for ii=1:length(nVector)
    
    logMessage(logfile,['FEAfix' int2str(nVector(ii)) ': start simulations']);
    
    tic;    
    dataSet = mot.dataSet;
    geo     = mot.geo;
    per     = mot.per;
    dataSet = back_compatibility(dataSet,geo,per,0);
    [~, ~, geo, per, mat] = data0(dataSet);
    mot.dataSet = dataSet;
    mot.geo = geo;
    mot.per = per;
    mot.mat = mat;
    
    dataSet.syrmDesignFlag.gf       = gfVector(ii);
    dataSet.syrmDesignFlag.ichf     = 0;
    dataSet.syrmDesignFlag.scf      = 0;
    dataSet.syrmDesignFlag.demag0   = 0;
    dataSet.syrmDesignFlag.demagHWC = 0;
    dataSet.syrmDesignFlag.mech     = 0;
    dataSet.syrmDesignFlag.therm    = 0;

    dataSet.FEAfixN = nVector(ii);

    
    eval_xbDesignPlane(dataSet,1);
    
    compTime = toc;
    
    logMessage(logfile,['FEAfix' int2str(nVector(ii)) ': evaluated in ' num2str(compTime) ' s']);


    close all

    logMessage(logfile,'Deleting temporary files...');
    try
        rmdir([cd,'\tmp'],'s');
    catch
    end
    if exist([cd,'\tmp'],'dir') == 0
        mkdir([cd,'\tmp']);
    end
    logMessage(logfile,'Temporary files deleted')
end

%% flux map test

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singm';
mot.dataSet.NumGrid          = 15;
mot.dataSet.CurrLoPP         = 3;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 30;
mot.dataSet.AngularSpanPP    = 60;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 20;
mot.dataSet.BrPP             = mot.mat.LayerMag.Br;

logMessage(logfile,'Flux map: start simulations');

tic
eval_fluxMap(mot.dataSet);
compTime = toc;
close all

logMessage(logfile,['Flux map: evaluated in ' num2str(compTime) ' s']);


logMessage(logfile,'Deleting temporary files...');
try
    rmdir([cd,'\tmp'],'s');
catch
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end
logMessage(logfile,'Temporary files deleted')

%% iron loss map

i0 = mot.per.i0;
mot.dataSet.EvalType         = 'singmIron';
mot.dataSet.NumGrid          = 9;
mot.dataSet.CurrLoPP         = 3;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 90;
mot.dataSet.AngularSpanPP    = 180;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;
mot.dataSet.MapQuadrants     = 1;
mot.dataSet.tempPP           = 20;
mot.dataSet.BrPP             = mot.mat.LayerMag.Br;

logMessage(logfile,'Loss map: start simulations');

tic
eval_fluxMap(mot.dataSet);
compTime = toc;
close all

logMessage(logfile,['Loss map: evaluated in ' num2str(compTime) ' s']);


logMessage(logfile,'Deleting temporary files...');
try
    rmdir([cd,'\tmp'],'s');
catch
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end
logMessage(logfile,'Temporary files deleted')

%% Demagnetization area

ipuVect = 10:10:100;

for ii=1:length(ipuVect)

    i0 = mot.per.i0;
    mot.dataSet.EvalType         = 'demagArea';
    mot.dataSet.NumGrid          = 9;
    mot.dataSet.CurrLoPP         = ipuVect(ii);
    mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
    mot.dataSet.NumOfRotPosPP    = 90;
    mot.dataSet.AngularSpanPP    = 180;
    mot.dataSet.currentpathname  = pathname;
    mot.dataSet.currentfilename  = filename;
    mot.dataSet.MapQuadrants     = 1;
    mot.dataSet.tempPP           = 20;
    mot.dataSet.BrPP             = mot.mat.LayerMag.Br;
    mot.dataSet.GammaPP          = 180;
    mot.dataSet.axisType         = 'PM';

    logMessage(logfile,['Demagnetization area #' int2str(ii) ' : start simulations']);

    tic
    eval_demagArea(mot.dataSet);
    compTime = toc;
    close all

    logMessage(logfile,['Demagnetization area #' int2str(ii) ' : evaluated in ' num2str(compTime) ' s']);

end


logMessage(logfile,'Deleting temporary files...');
try
    rmdir([cd,'\tmp'],'s');
catch
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end
logMessage(logfile,'Temporary files deleted')

%% Structural simulation

nVect = 2000:2000:20000;

for ii=1:length(nVect)

    mot.dataSet.EvalType  = 'structural';
    mot.dataSet.EvalSpeed = nVect(ii);

    logMessage(logfile,['Structural (PDE) simulation #' int2str(ii) ' : start simulations']);

    tic
    eval_vonMisesStress(mot.dataSet);
    compTime = toc;
    close all

    logMessage(logfile,['Structural (PDE) simulation #' int2str(ii) ' : evaluated in ' num2str(compTime) ' s']);

end

logMessage(logfile,'Deleting temporary files...');
try
    rmdir([cd,'\tmp'],'s');
catch
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end
logMessage(logfile,'Temporary files deleted')


%% end test

logMessage(logfile,'Timing test done!!!');




