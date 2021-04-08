function [cost,geo,mat,out,pathname] = MCADfitness (RQ,geo,per,mat,eval_type,filenameIn)

% [~,filename,ext] = fileparts(filenameIn);
[pathname,filename,ext] = fileparts(filenameIn);
filename = [filename ext]; % fem file name

pathname=[pathname '\'];

% [thisfilepath,pathname]=createTempDir();
% [thisfilepath,pathname]=createTempDir();

%load Syr-e and MCAD model
%load([pathname filename])
mcad=actxserver('MotorCAD.AppAutomation');
file_mot=[filename(1:(end-4)) '.mot'];
invoke(mcad,'LoadFromFile',[pathname file_mot]);

%MCAD mesh
invoke(mcad,'SetVariable','AirgapMeshPoints_layers',1440);
invoke(mcad,'SetVariable','AirgapMeshPoints_mesh',1440);

%EMag to Thermal
invoke(mcad,'SetVariable','MagneticThermalCoupling',1);

%AC losses
invoke(mcad,'SetVariable','ProximityLossModel',1);   %Full FEA model

[SOL] = simulate_xdegMCAD(geo,per,mat,eval_type,pathname,filename);

%save outputs
out.id = mean(SOL.id);  %const
out.iq = mean(SOL.iq);  %const
out.fd = mean(SOL.fd);  %waveform
out.fq = mean(SOL.fq);  %waveform
out.T = mean(SOL.T);    %waveform
out.dT = std(SOL.T);
out.dTpu = std(SOL.T)/out.T;
out.dTpp = max(SOL.T)-min(SOL.T);
out.IPF = SOL.IPF;
out.SOL = SOL;

%check Torque sign
if sign(out.T)~=sign(out.fd*out.iq-out.fq*out.id)
    out.T = -out.T;
    out.SOL.T = -out.SOL.T;
end

%save losses
out.Pfes_h = SOL.Pfes_h;
out.Pfes_c = SOL.Pfes_c;
out.Pfer_h = SOL.Pfer_h;
out.Pfer_c = SOL.Pfer_c;
out.Pfe    = out.Pfes_h+out.Pfes_c+out.Pfer_h+out.Pfer_c;
out.velDim = per.EvalSpeed;

%unused output
cost=[];

%save MCAD model
invoke(mcad,'SaveToFile',[pathname file_mot]);

%save MCAD results and quit
% invoke(mcad,'SaveResults','EMagnetic');
% invoke(mcad,'Quit');

end