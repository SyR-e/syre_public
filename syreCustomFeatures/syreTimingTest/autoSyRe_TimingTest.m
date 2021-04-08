%% autoSyRe_TimingTest.m

close all
clear
clc

[filename,pathname,~] = uigetfile([cd '\.mat'],'Load motor model');

mot = load([pathname filename]);

answer = inputdlg('Number of workers?','Setup',1,{'4'});
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
logname = ['LogFile_TimingTest', num2str(year) '_' num2str(month) '_' num2str(day) '_' num2str(hour) 'h_' num2str(mins) '.txt'];

logfile = [pathname logname];

fid = fopen(logfile,'w');
fprintf(fid,'Test computational time SyR-e\r\n');
fprintf(fid,['Date              : ' num2str(day) '/' num2str(month) '/' num2str(year) '\r\n']);
fprintf(fid,['Computer name     : ' getenv('computername') '\r\n']);
fprintf(fid,['Number of cores   : ' int2str(feature('numcores')) '\r\n']);
fprintf(fid,['Number of workers : ' int2str(numWorker) '\r\n']);
fprintf(fid,['Pathname          : ' strrep(pathname,'\','\\') '\r\n']);
fprintf(fid,['Motor model       : ' filename(1:end-4) '\r\n']);
% fprintf(fid,['Tester            : Simone Ferrari\r\n']);
fprintf(fid,'\r\n\r\n');
fclose(fid);
edit(logfile)

fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Parallel pool check...\r\n']);
fclose(fid);
edit(logfile)

% parallel pool creation
if ~isempty(isprop(gcp('nocreate'),'NumWorkers')) % no parallel pool enabled
    fid=fopen(logfile,'a');
    fprintf(fid,[datestr(now) '-->']);
    fprintf(fid,['Parallel pool enabled. Shutting down...\r\n']);
    fclose(fid);
    edit(logfile)
    
    delete(gcp)
    
    fid=fopen(logfile,'a');
    fprintf(fid,[datestr(now) '-->']);
    fprintf(fid,['Parallel pool disabled\r\n']);
    fclose(fid);
    edit(logfile)
end

parpool(numWorker);
fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Parallel pool enabled with ' int2str(numWorker) ' workers\r\n']);
fclose(fid);
edit(logfile)

%% FEAfix test: FEAfix4, FEAfix5, FEAfix8 and full plane

nVector = [1 4 5 8 1000];

for ii=1:length(nVector)
    
    fid=fopen(logfile,'a');
    fprintf(fid,[datestr(now) '-->']);
    fprintf(fid,['FEAfix' int2str(nVector(ii)) ': start simulations\r\n']);
    fclose(fid);
    edit(logfile)
    
    tic;    
    dataSet = mot.dataSet;
    dataSet = back_compatibility(dataSet,geo,per,0);
    [~, ~, geo, per, mat] = data0(dataSet);
    
    dataSet.FEAfixN = nVector(ii);
    
    % Design equations
    switch dataSet.TypeOfRotor
        case 'SPM'
            map = syrmDesign_SPM(dataSet);
        case 'Vtype'
            map = syrmDesign_Vtype(dataSet);
        otherwise
            map = syrmDesign_SyR(dataSet);
    end
    
    % FEAfix
    if dataSet.FEAfixN==0
        map.kd   = ones(size(map.xx));
        map.kq   = ones(size(map.xx));
        map.km   = ones(size(map.xx));
        map.k0   = ones(size(map.xx));
        map.xRaw = [];
        map.bRaw = [];
    else
        [FEAfixOut]=FEAfix(dataSet,geo,map);
        map.kd   = FEAfixOut.kd;
        map.kq   = FEAfixOut.kq;
        map.km   = FEAfixOut.km;
        map.k0   = FEAfixOut.k0;
        map.xRaw = FEAfixOut.xRaw;
        map.bRaw = FEAfixOut.bRaw;
    end
    
    if strcmp(dataSet.TypeOfRotor,'SPM')
        map.fd = map.fd.*map.kd;
        map.fq = map.fq.*map.kq;
        map.T  = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
        map.PF = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));
    elseif strcmp(dataSet.TypeOfRotor,'Vtype')
        map.fd  = map.fM.*map.km+(map.fd-map.fM).*map.kd;
        map.fq  = map.fq.*map.kq;
        map.fM  = map.fM.*map.km;
        map.T   = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
        map.ich = map.ich.*map.km./map.k0;
    else
        map.fd = map.fd.*map.kd;
        map.fq = map.fq.*map.kq;
        map.T  = 3/2*geo.p*(map.fd.*map.iq-map.fq.*map.id);
        map.PF = abs(sin(atan(map.iq./map.id)-atan(map.fq./map.fd)));
    end
    
    compTime = toc;
    
    fid=fopen(logfile,'a');
    fprintf(fid,[datestr(now) '-->']);
    fprintf(fid,['FEAfix' int2str(nVector(ii)) ': evaluated in ' num2str(compTime) ' s\r\n']);
    fclose(fid);
    edit(logfile)
    
    fid=fopen(logfile,'a');
    fprintf(fid,[datestr(now) '-->']);
    fprintf(fid,['Deleting temporary files...\r\n']);
    fclose(fid);
    edit(logfile)
    
    try
        rmdir([cd,'\tmp'],'s');
    catch
    end
    if exist([cd,'\tmp'],'dir') == 0
        mkdir([cd,'\tmp']);
    end
    
    fid=fopen(logfile,'a');
    fprintf(fid,[datestr(now) '-->']);
    fprintf(fid,['Temporary files deleted\r\n']);
    fclose(fid);
    edit(logfile)
    
    
end

%% flux map test

i0 = mot.per.i0;
mot.dataSet.EvalType = 'singm';
mot.dataSet.NumGrid  = 15;
mot.dataSet.CurrLoPP = 3;
mot.dataSet.SimulatedCurrent = mot.dataSet.CurrLoPP*i0;
mot.dataSet.NumOfRotPosPP    = 30;
mot.dataSet.currentpathname  = pathname;
mot.dataSet.currentfilename  = filename;

fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Flux map: start simulations\r\n']);
fclose(fid);
edit(logfile)

tic
eval_fluxMap(mot.dataSet);
compTime = toc;
close all

fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Flux map: evaluated in ' num2str(compTime) ' s\r\n']);
fclose(fid);
edit(logfile)


fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Deleting temporary files...\r\n']);
fclose(fid);
edit(logfile)

try
    rmdir([cd,'\tmp'],'s');
catch
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end

fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Temporary files deleted\r\n']);
fclose(fid);
edit(logfile)

fid=fopen(logfile,'a');
fprintf(fid,[datestr(now) '-->']);
fprintf(fid,['Timing test done!!!\r\n']);
fclose(fid);
edit(logfile)



