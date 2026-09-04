%% autoSyre_PMsegmentation_map.m

close all
clear
clc


[fileName, filePath] = uigetfile('*.mat', 'Select source .mat file');
sourceFile = fullfile(filePath, fileName);

[pathname,~,~] = fileparts(sourceFile);
pathname = checkPathSyntax([pathname '\']);

load(sourceFile);

logfile = [pathname 'PMsegmentation_loss_sensitivity_logfile.txt'];

prompt = {...
    'Minimum axial segments:',...
    'Maximum axial segments:',...
    'Minimum tangential segments:',...
    'Maximum tangential segments:',...
    'Number of workers:'};
dlgTitle = 'PM segmentation inputs';
defaultAns = {'1','134','1','20','16'};
answer = inputdlg(prompt, dlgTitle, 1, defaultAns);


min_nL     = round(str2double(answer{1}));
max_nL     = round(str2double(answer{2}));
min_nT     = round(str2double(answer{3}));
max_nT     = round(str2double(answer{4}));
numWorkers = round(str2double(answer{5}));

nL = min_nT:1:max_nL;
nT = min_nT:1:max_nT;


[nL,nT] = meshgrid(nL,nT);

Ppm = nan(size(nL));


% logfile header
logHeader(logfile,'PM Segmentation - loss analysis')
fid = fopen(logfile,'a');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'PM loss sensitivity to PM segmentation\r\n');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Evaluation details:\r\n');
fprintf(fid,['Number of cores   : ' int2str(feature('numcores')) '\r\n']);
fprintf(fid,['Number of workers : ' int2str(numWorkers) '\r\n']);
fprintf(fid,['Pathname          : ' strrep(pathname,'\','\\') '\r\n']);
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,['Minimum axial segments      : ' int2str(min(nL(:))) '\r\n']);
fprintf(fid,['Maximum axial segments      : ' int2str(max(nL(:))) '\r\n']);
fprintf(fid,['Minimum tangential segments : ' int2str(min(nT(:))) '\r\n']);
fprintf(fid,['Maximum tangential segments : ' int2str(max(nT(:))) '\r\n']);
fprintf(fid,['Number of configuration     : ' int2str(numel(nL)) '\r\n']);
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fclose(fid);

% parallel pool check and enabling
logMessage(logfile,'Parallel pool check...')
ppState = parallelComputingCheck();
if ppState==-1
    logMessage(logfile,'Parallel Computing Toolbox not installed. Evaluation without parallel computing.')
    ppEnable = 0;
    numWorkers = 1;
elseif ppState==0
    logMessage(logfile,'Parallel pool not enabled')
    ppEnable = 1;
elseif ppState~=numWorkers
    logMessage(logfile,['Parallel pool enabled on ' int2str(ppState) ' workers. Shutting down...'])
    delete(gcp)
    logMessage(logfile,'Parallel pool disabled')
    ppEnable = 1;
else
    logMessage(logfile,['Parallel pool already enabled on ' int2str(ppState) ' workers']);
    ppEnable = 0;
end

if ppEnable&&numWorkers>1
    logMessage(logfile,'Enabling parallel pool...')
    parpool(numWorkers);
    logMessage(logfile,['Parallel pool enabled on ' int2str(numWorkers) ' workers'])
end

logMessage(logfile,'Start evaluation')

Ppm = nan(size(nL));
if numWorkers==1
    for ii=1:numel(nL)
        [~,Ppm(ii)] = recompute_PMloss(geo,per,mat,out,nL(ii),nT(ii));
        logMessage(logfile,['Computed point: nL=' sprintf('%04d',nL(ii)) ' / nT=' sprintf('%04d',nT(ii)) ' --> Ppm = ' sprintf('%7.2f',Ppm(ii)) ' W'],1)
    end
else
    parfor ii=1:numel(nL)
        geo_local = geo;
        per_local = per;
        mat_local = mat;
        out_local = out;
        [~,Ppm(ii)] = recompute_PMloss(geo_local,per_local,mat_local,out_local,nL(ii),nT(ii));
        logMessage(logfile,['Computed point: nL=' sprintf('%04d',nL(ii)) ' / nT=' sprintf('%04d',nT(ii)) ' --> Ppm = ' sprintf('%7.2f',Ppm(ii)) ' W'],0);
    end


end

logMessage(logfile,'Evaluation done!')

hfig(1) = figure();
figSetting();
set(hfig(1),'FileName',[pathname 'PMloss_segmentationSense3D.fig']);
xlabel('axial segments')
ylabel('tangential segments (per PM)')
zlabel('PM loss (W)')
colormap turbo
view(3)
stem3(nL,nT,Ppm,'k')
scatter3(nL(:),nT(:),Ppm(:),[],Ppm(:),'filled')
% surf(nL,nT,Ppm)


hfig(2) = figure();
figSetting();
set(hfig(2),'FileName',[pathname 'PMloss_segmentationSense2D.fig']);
xlabel('axial segments')
ylabel('tangential segments (per PM)')
colormap turbo
colorbar
clim([min(Ppm(:)) max(Ppm(:))]);
contourf(nL,nT,Ppm,'EdgeColor','none','LevelList',linspace(min(Ppm(:)),max(Ppm(:)),501),'DisplayName','$P_{PM}$ (W)');
contour(nL,nT,Ppm,'-k','ShowText','on','HandleVisibility','off')
plot(nL(:),nT(:),'r.','DisplayName','tested configurations')
legend('show','Location','northeast')

for ii=1:length(hfig)
    savePrintFigure(hfig(ii))
end



























