function [SOL]=simulate_xdegMCAD(geo,per,mat,eval_type,pathname,filename)

filename=[filename(1:end-4) '.mat'];
load([pathname filename])
mcad=actxserver('MotorCAD.AppAutomation');
invoke(mcad,'SetVariable','PhaseAdvance',dataSet.GammaPP);           % phase advance
io=calc_io(geo,per);
invoke(mcad,'SetVariable','PeakCurrent',dataSet.CurrLoPP*io);   % peak current
% tmp=input('Insert DC voltage bus [V]: ','s');
% tmp(tmp=='.')=',';
invoke(mcad,'SetVariable','DCBusVoltage',565);   
% invoke(mcad,'SetVariable','Shaft_Speed_Ref',1000);   
% invoke(mcad,'SetVariable','Shaft_Speed_[RPM]',1000);   

if dataSet.EvalSpeed~=0
    invoke(mcad,'SetVariable','Shaft_Speed_Ref',dataSet.EvalSpeed);
else  invoke(mcad,'SetVariable','Shaft_Speed_[RPM]',1000); 
    %invoke(mcad,'SetVariable','Shaft_Speed_Ref',1000);
    disp('simulation runs with a default value of 1000 rpm - No input speed from Syr-e')
end
invoke(mcad,'SetVariable','ArmatureConductor_Temperature',per.tempcu);

%Simulation settings

invoke(mcad,'SetVariable','BackEMFCalculation','False');
invoke(mcad,'SetVariable','CoggingTorqueCalculation','False');
invoke(mcad,'SetVariable','TorqueSpeedCalculation','False');
invoke(mcad,'SetVariable','DemagnetizationCalc','False');
invoke(mcad,'SetVariable','TorqueCalculation','True');
nPoints=dataSet.NumOfRotPosPP*6;                                 %over 360 eltDeg
invoke(mcad,'SetVariable','TorquePointsPerCycle',int2str(nPoints));
magnetic_solver=1;                                               %multi-static magnetic solver
invoke(mcad,'SetVariable','MagneticSolver',magnetic_solver);     %multi-static magnetic solver
invoke(mcad,'SetVariable','ArmatureEWdgMLT_Multiplier',0);       %no end-windings effect
invoke(mcad,'SetVariable','MagThreads_Option',0);                %single or multiple threads (0 or 1)

% disp('Magnetic simulation in progress...')
success=invoke(mcad,'DoMagneticCalculation'); 
if success==0     
    disp('Magnetic calculation successfully completed')  
else
    disp('Magnetic calculation failed') 
end

%save losses
if magnetic_solver==0
[tmp,Pfes_h_BackIron]=invoke(mcad,'GetVariable','StatorBackIronLoss_Hys');
[tmp,Pfes_h_Tooth]=invoke(mcad,'GetVariable','StatorToothLoss_Hys');
[tmp,Pfes_exc_backiron]=invoke(mcad,'GetVariable','StatorBackIronLoss_Excess');
[tmp,Pfes_exc_tooth]=invoke(mcad,'GetVariable','StatorToothLoss_Excess');
SOL.Pfes_h=Pfes_h_BackIron+Pfes_h_Tooth+Pfes_exc_backiron+Pfes_exc_tooth;

[tmp,Pfes_c_BackIron]=invoke(mcad,'GetVariable','StatorBackIronLoss_Eddy');
[tmp,Pfes_c_Tooth]=invoke(mcad,'GetVariable','StatorToothLoss_Eddy');
SOL.Pfes_c=Pfes_c_BackIron+Pfes_c_Tooth;

[tmp,Pfer_h_BackIron]=invoke(mcad,'GetVariable','RotorBackIronLoss_Hys');
[tmp,Pfer_h_Tooth]=invoke(mcad,'GetVariable','RotorMagnetPoleLoss_Hys');
[tmp,Pfer_exc_backiron]=invoke(mcad,'GetVariable','RotorBackIronLoss_Excess');
[tmp,Pfer_exc_tooth]=invoke(mcad,'GetVariable','RotorMagnetPoleLoss_Excess');
SOL.Pfer_h=Pfer_h_BackIron+Pfer_h_Tooth+Pfer_exc_tooth+Pfer_exc_backiron;

[tmp,Pfer_c_BackIron]=invoke(mcad,'GetVariable','RotorBackIronLoss_Eddy');
[tmp,Pfer_c_Tooth]=invoke(mcad,'GetVariable','RotorMagnetPoleLoss_Eddy');
SOL.Pfer_c=Pfer_c_BackIron+Pfer_c_Tooth;
end

if magnetic_solver==1
[tmp,Pfes_h_BackIron]=invoke(mcad,'GetVariable','StatorBackIronLoss_Hys_Static');
[tmp,Pfes_h_Tooth]=invoke(mcad,'GetVariable','StatorToothLoss_Hys_Static');
[tmp,Pfes_exc_backiron]=invoke(mcad,'GetVariable','StatorBackIronLoss_Exc_Static');
[tmp,Pfes_exc_tooth]=invoke(mcad,'GetVariable','StatorToothLoss_Exc_Static');
SOL.Pfes_h=Pfes_h_BackIron+Pfes_h_Tooth+Pfes_exc_backiron+Pfes_exc_tooth;

[tmp,Pfes_c_BackIron]=invoke(mcad,'GetVariable','StatorBackIronLoss_Eddy_Static');
[tmp,Pfes_c_Tooth]=invoke(mcad,'GetVariable','StatorToothLoss_Eddy_Static');
SOL.Pfes_c=Pfes_c_BackIron+Pfes_c_Tooth;

SOL.Pfer_c=0;
SOL.Pfer_h=0;
end
%save current dq (dq axis Syr-e)
[tmp,SOL.id]=invoke(mcad,'GetVariable','CurrentLoad_Q');
[tmp,SOL.iq]=invoke(mcad,'GetVariable','CurrentLoad_D');  
SOL.iq=-sqrt(2)*SOL.iq;      SOL.id=sqrt(2)*SOL.id;

%save IPF
[tmp,SOL.IPF]=invoke(mcad,'GetVariable','WaveformPowerFactor');

%%salvo perdite magnetiche che serviranno per calcolo termico
% [tmp,AmatureLoss]=invoke(mcad,'GetVariable','ConductorLoss');
% [tmp,MagnetLoss]=invoke(mcad,'GetVariable','MagnetLoss');
% [tmp,StatorBackIronLoss]=invoke(mcad,'GetVariable','StatorBackIronLoss_Total');
% [tmp,StatorToothLoss]=invoke(mcad,'GetVariable','StatorToothLoss_Total');
% [tmp,RotorBackIronLoss]=invoke(mcad,'GetVariable','RotorBackIronLoss_Total');

% save([cd, '\results\MCAD\loss_values.mat'], 'AmatureLoss','MagnetLoss','StatorBackIronLoss','StatorToothLoss','RotorBackIronLoss');

% disp('Magnetic results saved in:')
% disp([cd '\loss.mat'])
% disp(' ')

% loss=[AmatureLoss,MagnetLoss,StatorBackIronLoss,StatorToothLoss,RotorBackIronLoss];
% c=categorical({'Conductor','Magnet','Stator Back Iron','Stator Tooth','Rotor Back Iron'});
% figure()
% bar(c,loss)
% title('Magnetic loss')
% ylabel('Loss [Watt]')
% savefig([cd, '\results\MCAD\loss_graph.fig']);

%save Torque
RotorPosition = linspace(0,360,nPoints); 
SLOT.T =zeros(nPoints,1); 
for loop=1:nPoints 
    [success,x,y]=invoke(mcad,'GetMagneticGraphPoint','TorqueVW',loop);
    if success == 0       
        RotorPosition(loop)=x; 
        SOL.T(loop)=y;     
    end
end

%save flux dq
for loop=1:nPoints
    [success,x,y]=invoke(mcad,'GetMagneticGraphPoint','FluxLinkageLoadTotalD',loop);
    if success == 0       
        RotorPosition(loop)=x; 
        SOL.fq(loop)=-y;     
    end
end
for loop1=1:nPoints
    [success,x,y]=invoke(mcad,'GetMagneticGraphPoint','FluxLinkageLoadTotalQ',loop1);
    if success == 0       
        RotorPosition(loop1)=x; 
        SOL.fd(loop1)=y;     
    end
end

end