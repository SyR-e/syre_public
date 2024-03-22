% Copyright 2023
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

function [results] = autoSyre_checkMotor(dataIn)

close all
clc

if nargin()==0
    [filename,pathname] = uigetfile([cd '\*.mat'],'Select motor');
    load([pathname filename],'dataSet');
    dataIn = dataSet;
    clear dataSet;
end

prompt={
    'Peak phase current [Apk]'
    'DC link voltage [Vdc]'
    'PM temperatures [degC]'
    'Cu temperature [degC]'};
name='Inverter limits';
numlines=1;
answer={
    num2str(dataIn.RatedCurrent)
    '565'
    '[20 60 100]'
    '[20 100 120]'};

answer=inputdlg(prompt,name,numlines,answer);

results.Imax        = eval(answer{1})*ones(1,3);
results.Vdc         = eval(answer{2})*ones(1,3);
results.tempPM      = eval(answer{3});
results.tempCu      = eval(answer{4});
results.dataSet     = dataIn;
results.Idq         = nan(1,3);
results.Vdq         = nan(1,3);
results.Fdq         = nan(1,3);
results.IHWC0       = nan(1,3);
results.IHWCmax     = nan(1,3);
results.Fm          = nan(1,3);
results.T           = nan(1,3);
results.dTpp        = nan(1,3);
results.nbase       = nan(1,3);
results.nmax        = dataIn.OverSpeed*ones(1,3);
results.nUGO        = nan(1,3);
results.dPM_Imax    = nan(1,3);
results.dPM_HWC0    = nan(1,3);
results.dPM_HWCmax  = nan(1,3);
results.Bmin_Imax   = nan(1,3);
results.Bmin_HWC0   = nan(1,3);
results.Bmin_HWCmax = nan(1,3);
results.Ich         = nan(1,3);


mat = material_properties_layer(dataIn.FluxBarrierMaterial);

for ii=1:length(results.tempPM)
    if results.tempPM(ii)<min(mat.temp.temp)
        results.tempPM(ii) = NaN;
        results.tempCu(ii) = NaN;
    end
    if results.tempPM(ii)>max(mat.temp.temp)
        results.tempPM(ii) = NaN;
        results.tempCu(ii) = NaN;
    end
end

results.tempPM = results.tempPM(~isnan(results.tempPM));
results.tempCu = results.tempCu(~isnan(results.tempCu));

results.Br = interp1(mat.temp.temp,mat.temp.Br,results.tempPM);
results.Bd = interp1(mat.temp.temp,mat.temp.Bd,results.tempPM);

results.dataSet.EvalType = 'singt';

resFolder = [dataIn.currentpathname dataIn.currentfilename(1:end-4) '_results\preliminaryTest\'];
mkdir(resFolder)
logfile    = [resFolder 'logfile.txt'];
reportfile = [resFolder 'report.txt'];
timeIni = datetime('now');

% logfile initialization
fid = fopen(logfile,'w');
fprintf(fid,'Preliminary motor check - logfile\r\n');
fprintf(fid,['Date          : ' char(timeIni) '\r\n']);
fprintf(fid,['Motor name    : ' dataIn.currentfilename(1:end-4) '\r\n']);
fprintf(fid,['Pathname      : ' strrep(dataIn.currentpathname,'\','\\') '\r\n']);
fprintf(fid,['Computer name : ' getenv('computername') '\r\n']);
fprintf(fid,'\r\n\r\n');
fclose(fid);
edit(logfile)

% parallel computing check
if ~isempty(isprop(gcp('nocreate'),'NumWorkers')) % no parallel pool enabled
    logMessage(logfile,'Parallel computing enabled. Shutting down...');
    delete(gcp)
    logMessage(logfile,'Parallel computing disabled')
end

parpool(feature('numcores'));
logMessage(logfile,['Parallel computing enabled on ' int2str(feature('numcores')) ' workers'])


for ii=1:length(results.tempPM)
    logMessage(logfile,['Starting evaluation at ' int2str(results.tempPM(ii)) 'degC']);
    
    dataSet = results.dataSet;
    dataSet.tempPP = results.tempPM(ii);
    mat = material_properties_layer(dataSet.FluxBarrierMaterial);
    dataSet.BrPP = interp1(mat.temp.temp,mat.temp.Br,dataSet.tempPP);
    % MTPA search
    dataSet.SimulatedCurrent = results.Imax(ii)*ones(1,24);
    dataSet.CurrLoPP         = results.Imax(ii)/dataSet.RatedCurrent*ones(1,24);
    if strcmp(dataSet.axisType,'SR')
        dataSet.GammaPP = linspace(0,90,24);
    else
        dataSet.GammaPP = linspace(90,180,24);
    end
    dataSet.NumOfRotPosPP = 12;
    dataSet.AngularSpanPP = 60;
    logMessage(logfile,'Start MTPA search...')
    [~,senseOut] = eval_operatingPoint(dataSet);
    logMessage(logfile,'MTPA identified')
    [~,index] = max(senseOut.T);
    results.Idq(ii) = senseOut.id(index)+j*senseOut.iq(index);
    close all

    % single point
    dataSet.SimulatedCurrent = dataSet.SimulatedCurrent(1);
    dataSet.CurrLoPP         = dataSet.CurrLoPP(1);
    dataSet.GammaPP          = angle(results.Idq(ii))*180/pi;
    dataSet.NumOfRotPosPP    = 30;
    dataSet.AngularSpanPP    = 60;
    logMessage(logfile,'Start single point simulation')
    [out,~] = eval_operatingPoint(dataSet);
    logMessage(logfile,'Single point simulated')
    results.Fdq(ii)  = out.fd+j*out.fq;
    results.T(ii)    = out.T;
    results.dTpp(ii) = out.dTpp;
    close all

    % HWC computation
    dataSet.NumOfRotPosPP    = 12;
    dataSet.AngularSpanPP    = 60;
    dataSet.CurrLoPP         = [0 dataSet.CurrLoPP];
    dataSet.SimulatedCurrent = [0 dataSet.SimulatedCurrent];
    dataSet.GammaPP          = [0 dataSet.GammaPP];
    logMessage(logfile,'Start HWC current computation')
    [pkSCout] = eval_peakShortCircuitCurrent(dataSet);
    logMessage(logfile,'HWC current computation done')
    results.IHWC0(ii)   = abs(pkSCout.idq(1));
    results.IHWCmax(ii) = abs(pkSCout.idq(2));
    close all

    % demagnetization @ Imax
    dataSet.CurrLoPP         = dataSet.CurrLoPP(2);
    dataSet.SimulatedCurrent = dataSet.SimulatedCurrent(2);
    if strcmp(dataSet.axisType,'SR')
        dataSet.GammaPP = 90;
    else
        dataSet.GammaPP = 180;
    end
    logMessage(logfile,'Start demagnetization @ Imax')
    [DemagArea] = eval_demagArea(dataSet);
    logMessage(logfile,'Demagnetization @ Imax computed')
    results.dPM_Imax(ii)  = max(DemagArea.dPM);
    results.Bmin_Imax(ii) = min(DemagArea.Bmin);
    close all

    % demagnetization @ HWC from no-load
    dataSet.CurrLoPP         = results.IHWC0(ii)/dataSet.RatedCurrent;
    dataSet.SimulatedCurrent = results.IHWC0(ii);
    if strcmp(dataSet.axisType,'SR')
        dataSet.GammaPP = 90;
    else
        dataSet.GammaPP = 180;
    end
    logMessage(logfile,'Start demagnetization @ HWC current from no-load')
    [DemagArea] = eval_demagArea(dataSet);
    logMessage(logfile,'Demagnetization @ HWC current from no-load computed')
    results.dPM_HWC0(ii)  = max(DemagArea.dPM);
    results.Bmin_HWC0(ii) = min(DemagArea.Bmin);
    close all

    % demagnetization @ HWC from Imax
    dataSet.CurrLoPP         = results.IHWCmax(ii)/dataSet.RatedCurrent;
    dataSet.SimulatedCurrent = results.IHWCmax(ii);
    if strcmp(dataSet.axisType,'SR')
        dataSet.GammaPP = 90;
    else
        dataSet.GammaPP = 180;
    end
    logMessage(logfile,'Start demagnetization @ HWC current from Imax')
    [DemagArea] = eval_demagArea(dataSet);
    logMessage(logfile,'Demagnetization @ HWC current from Imax computed')
    results.dPM_HWCmax(ii)  = max(DemagArea.dPM);
    results.Bmin_HWCmax(ii) = min(DemagArea.Bmin);
    close all

    % characteristic current
    dataSet.NumOfRotPosPP    = 12;
    dataSet.AngularSpanPP    = 60;
    logMessage(logfile,'Start evaluation of characteristic current')
    [ichOut,~] = eval_ich(dataSet);
    logMessage(logfile,'Charateristic current computed')
    results.Ich(ii) = ichOut.ichVect;
    results.Fm(ii)  = ichOut.FmVect;
    close all

    logMessage(logfile,['Simulations at ' int2str(results.tempPM(ii)) 'degC done!'])
end

% computation of the phase resistance
logMessage(logfile,'Resistance computation...')
tmp = load([dataIn.currentpathname dataIn.currentfilename]);
for ii=1:length(results.tempCu)
    tmp.per.tempcu = results.tempCu(ii);
    [per] = calc_i0(tmp.geo,tmp.per);
    results.Rs(ii) = per.Rs;
    clear per
end
logMessage(logfile,'Phase resistance computed')

% post-processing data
logMessage(logfile,'Data post-processing...')
results.wbase   = calcLimitPulsation(real(results.Idq),imag(results.Idq),real(results.Fdq),imag(results.Fdq),results.Rs,results.Vdc(1)/sqrt(3));
results.nbase   = results.wbase/dataSet.NumOfPolePairs*30/pi;
results.nUGO    = results.nbase.*abs(results.Fdq)./results.Fm;
results.Vdq     = results.Rs.*results.Idq+j*results.wbase.*results.Fdq;
results.V0nbase = results.Fm.*results.wbase;
results.V0nmax  = results.Fm.*results.nmax*pi/30*dataSet.NumOfPolePairs;
results.IPF     = sin(angle(results.Idq)-angle(results.Fdq));
results.PF      = cos(angle(results.Vdq)-angle(results.Idq));
results.P       = results.T.*results.nbase*pi/30;
logMessage(logfile,'Data post-processing done')

save([resFolder 'preliminaryCheckResults.mat'],'results');

% Report generation
logMessage(logfile,'Report file generation...')

fid = fopen(reportfile,'w');
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Preliminary motor check\r\n');
fprintf(fid,['Date          : ' char(timeIni) '\r\n']);
fprintf(fid,['Motor name    : ' dataIn.currentfilename(1:end-4) '\r\n']);
fprintf(fid,['Pathname      : ' strrep(dataIn.currentpathname,'\','\\') '\r\n']);
fprintf(fid,['Computer name : ' getenv('computername') '\r\n']);
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
fprintf(fid,'Inverter limits:\r\n');
fprintf(fid,['Peak phase current : ' num2str(results.Imax(1),'% 5.2f') ' Apk (' num2str(results.Imax(1)/sqrt(2),'% 5.2f') ' Arms)\r\n']);
fprintf(fid,['DC link voltage    : ' num2str(results.Vdc(1),'% 5.2f') ' V\r\n']);
fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
for ii=1:length(results.tempPM)
    fprintf(fid,['Performance figures at ' int2str(results.tempPM(ii)) 'degC PM / ' int2str(results.tempCu(ii)) 'degC Cu\r\n']);
    fprintf(fid,['d-axis current                            : ' num2str(real(results.Idq(ii)),'% 5.2f') ' Apk\r\n']);
    fprintf(fid,['q-axis current                            : ' num2str(imag(results.Idq(ii)),'% 5.2f') ' Apk\r\n']);
    fprintf(fid,['Current angle                             : ' num2str(angle(results.Idq(ii))*180/pi,'% 5.2f') '\r\n']);
    fprintf(fid,['Flux linkage amplitude                    : ' num2str(abs(results.Fdq(ii)),'% 5.2f') ' Vs\r\n']);
    fprintf(fid,['Flux linkage angle                        : ' num2str(angle(results.Fdq(ii))*180/pi,'% 5.2f') '\r\n']);
    fprintf(fid,['d-axis flux linkage                       : ' num2str(real(results.Fdq(ii)),'% 5.2f') ' Vs\r\n']);
    fprintf(fid,['q-axis flux linkage                       : ' num2str(imag(results.Fdq(ii)),'% 5.2f') ' Vs\r\n']);
    fprintf(fid,['PM flux linkage                           : ' num2str(real(results.Fm(ii)),'% 5.2f') ' Vs (' int2str(results.Fm(ii)/abs(results.Fdq(ii))*100) '%% of fdq)\r\n']);
    fprintf(fid,['Peak phase voltage                        : ' num2str(abs(results.Vdq(ii)),'% 5.2f') ' Vpk\r\n']);
    fprintf(fid,['d-axis voltage                            : ' num2str(real(results.Vdq(ii)),'% 5.2f') ' Vpk\r\n']);
    fprintf(fid,['q-axis voltage                            : ' num2str(imag(results.Vdq(ii)),'% 5.2f') ' Vpk\r\n']);
    fprintf(fid,['Internal power factor                     : ' num2str(results.IPF(ii),'% 5.2f') '\r\n']);
    fprintf(fid,['Power factor                              : ' num2str(results.PF(ii),'% 5.2f') '\r\n']);
    fprintf(fid,['Phase resistance                          : ' num2str(results.Rs(ii),'% 5.2f') ' Ohm\r\n']);
    fprintf(fid,['Torque                                    : ' num2str(abs(results.T(ii)),'% 5.2f') ' Nm\r\n']);
    fprintf(fid,['Power                                     : ' num2str(abs(results.P(ii)/1000),'% 5.2f') ' kW\r\n']);
    fprintf(fid,['Base speed                                : ' num2str(real(results.nbase(ii)),'% 5.2f') ' rpm\r\n']);
    fprintf(fid,['UGO limit                                 : ' num2str(real(results.nUGO(ii)),'% 5.2f') ' rpm (' int2str(results.nUGO(ii)/results.nbase(ii)*100) '%% of nbase)\r\n']);
    fprintf(fid,['Maximum speed                             : ' num2str(real(results.nmax(ii)),'% 5.2f') ' rpm (' int2str(results.nmax(ii)/results.nbase(ii)*100) '%% of nbase)\r\n']);
    fprintf(fid,['No-load line voltage @ nbase              : ' num2str(abs(results.V0nbase(ii)*sqrt(3)),'% 5.2f') ' Vpk (' int2str(results.V0nbase(ii)*sqrt(3)/results.Vdc(ii)*100) '%% of Vdc)\r\n']);
    fprintf(fid,['No-load line voltage @ nmax               : ' num2str(abs(results.V0nmax(ii)*sqrt(3)),'% 5.2f') ' Vpk (' int2str(results.V0nmax(ii)*sqrt(3)/results.Vdc(ii)*100) '%% of Vdc)\r\n']);
    fprintf(fid,['Characteristic current                    : ' num2str(results.Ich(ii),'% 5.2f') ' Apk (' int2str(results.Ich(ii)/results.Imax(ii)*100) '%% of Imax)\r\n']);
    fprintf(fid,['HWC short-circuit current (from no-load)  : ' num2str(results.IHWC0(ii),'% 5.2f') ' Apk (' int2str(results.IHWC0(ii)/results.Imax(ii)*100) '%% of Imax)\r\n']);
    fprintf(fid,['HWC short-circuit current (from Imax)     : ' num2str(results.IHWCmax(ii),'% 5.2f') ' Apk (' int2str(results.IHWCmax(ii)/results.Imax(ii)*100) '%% of Imax)\r\n']);
    fprintf(fid,['Demagnetized PM volume @ Imax             : ' num2str(real(results.dPM_Imax(ii)*100),'% 5.2f') '%%\r\n']);
    fprintf(fid,['Minimum flux density @ demag Imax         : ' num2str(results.Bmin_Imax(ii),'% 5.2f') ' T (' int2str(results.Bmin_Imax(ii)/(results.Br(ii)-results.Bd(ii))*100) '%% of the available range)\r\n']);
    fprintf(fid,['Demagnetized PM volume @ HWC from no-load : ' num2str(real(results.dPM_HWC0(ii)*100),'% 5.2f') '%%\r\n']);
    fprintf(fid,['Minimum flux density @ demag HWC no-load  : ' num2str(results.Bmin_HWC0(ii),'% 5.2f') ' T (' int2str(results.Bmin_HWC0(ii)/(results.Br(ii)-results.Bd(ii))*100) '%% of the available range)\r\n']);
    fprintf(fid,['Demagnetized PM volume @ HWC from Imax    : ' num2str(real(results.dPM_HWCmax(ii)*100),'% 5.2f') '%%\r\n']);
    fprintf(fid,['Minimum flux density @ demag HWC Imax     : ' num2str(results.Bmin_HWCmax(ii),'% 5.2f') ' T (' int2str(results.Bmin_HWCmax(ii)/(results.Br(ii)-results.Bd(ii))*100) '%% of the available range)\r\n']);
    fprintf(fid,'-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-\r\n');
end

fclose(fid);

logMessage(logfile,'Report file created');
edit(reportfile)

% delete temporary files and parallel computing

logMessage(logfile,'Deleting temporary files...')
try
    rmdir([cd,'\tmp'],'s');
catch
end
if exist([cd,'\tmp'],'dir') == 0
    mkdir([cd,'\tmp']);
end
logMessage(logfile,'Temporary files deleted')

logMessage(logfile,'Disabling parallel computing...')
delete(gcp)
logMessage(logfile,'Parallel computing disabled')

logMessage(logfile,'Motor performance chech completed!')

if nargout()==0
    clear results
end
